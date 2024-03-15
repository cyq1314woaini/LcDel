# coding=utf-8
"""
查询候选deletion位点
输入：bam文件
输出：
步骤：
1、首先根据cigar字符串查找比对内的deletion，将其表示为四元组【chron, pos, length, d-cigar】
2、分别查找是补充比对的和没有补充比对但具有SA字段的将其表示为四元组【chrom, pos, length, s-cigar】
3、按照在参考基因组上的位置进行排序，采用一定的算法进行聚类
"""

import pysam
import multiprocessing as mp
import os
import numpy as np
import pandas as pd
import sys


def decode_flag(Flag):
    signal = {1 << 2: 0, 1 >> 1: 1, 1 << 4: 2, 1 << 11: 3, 1 << 4 | 1 << 11: 4}
    return signal[Flag] if (Flag in signal) else 0


def find_read_length(chrom, bam, outpath):
    #     计算比对到chrom染色体上的读数长度，这里需要过滤掉没有比对上以及二次比对的读数
    totalMapLength = 0
    bamfile = pysam.AlignmentFile(bam, "r")
    allreads = bamfile.fetch(contig=chrom)
    for read in allreads:
        if read.is_unmapped or read.is_secondary:
            continue
        if read.is_supplementary:
            continue
        else:
            totalMapLength += read.query_length

    # 将计算的结果写入到文件中
    if totalMapLength != 0:
        f = open(outpath, "w")
        f.write(str(totalMapLength) + "\n")
        f.close()
    print("比对到{}染色体的读数总长度为{}".format(chrom, totalMapLength))


def find_coverage(bam, max_work, output):
    bamfile = pysam.AlignmentFile(bam, "r")

    chromosomes = bamfile.references
    chromolength = {}
    for chrom in chromosomes:
        tempchrom = chrom
        # if "chr" in chrom:
        #     tempchrom = tempchrom[3:]
        chromolength[tempchrom] = bamfile.get_reference_length(chrom)

    coveragePath = os.path.join(output, "coverageTemp")
    if not os.path.exists(coveragePath):
        os.mkdir(coveragePath)

        pool = mp.Pool(processes=max_work)
        for chrom in chromolength:
            chromCoveragePath = os.path.join(coveragePath, "coverage_{}.txt".format(chrom))
            pool.apply_async(find_read_length, args=(chrom, bam, chromCoveragePath))

        pool.close()
        pool.join()

    print("=============================================")
    print("比对到各染色体上的读数总长度已经计算完毕~~~")
    os.system("cat " + coveragePath + "/coverage_*.txt >" + coveragePath + "/coverage_allchrom.txt")

    readlengthlist = []
    f = open(coveragePath + "/coverage_allchrom.txt", "r")
    while True:
        lengt = f.readline().rstrip()
        if not lengt:
            break
        readlengthlist.append(int(lengt))
    # 计算染色体的总长度以及比对到其染色体上的总长度
    reflength = sum([chromolength[c] for c in chromolength])
    alignreadlength = sum(readlengthlist)
    depth = float(alignreadlength) / reflength

    return depth


def findlength(bam_path):
    bamfile = pysam.AlignmentFile(bam_path, "r")
    references = bamfile.references
    chromosomelength = {}

    for chrom in references:
        tempchrom = chrom
        # if "chr" in tempchrom:
        #     tempchrom = tempchrom[3:]
        length = bamfile.get_reference_length(chrom)
        if length >= 40000000:
            chromosomelength[tempchrom] = length

    return chromosomelength


def find_cigar_pos(position, cigar, chr):
    resultsv = []
    numberstr = "0123456789"
    num = ""
    reflen = 0
    readlen = 0
    leftclip = 0
    rightclip = 0
    gapthreshold = 30
    min_size = 40
    for c in cigar:
        if c in numberstr:
            num += c
        if c in 'MNP=X':
            readlen += int(num)
            reflen += int(num)
            num = ""
            continue
        if c == 'D':
            if int(num) >= min_size:
                resultsv += [[chr, position + reflen, int(num), 'd-cigar']]
            reflen += int(num)
            num = ""
            continue
        if c == "I":
            readlen += int(num)
            num = ""
            continue
        if c in "SH":
            if readlen == 0:
                leftclip = int(num)
            else:
                rightclip = int(num)
            num = ""
            continue

    # 合并
    if len(resultsv) > 1:
        index = len(resultsv) - 1
        while index > 0:
            rightstart = resultsv[index][1]
            leftend = resultsv[index - 1][1] + resultsv[index - 1][2]
            gap = rightstart - leftend
            if 0 < gap <= gapthreshold:
                resultsv = resultsv[:index - 1] + [
                    [resultsv[index][0], resultsv[index - 1][1], resultsv[index - 1][2] + resultsv[index][2],
                     'd-cigar']] + resultsv[index + 1:]
            index = index - 1

    return [resultsv, reflen, [leftclip, readlen, rightclip]]


def c_pos(cigar, temprefstart):
    number = ""
    numstr = "1234567890"
    readstart = False
    readend = False
    refend = False
    readloc = 0
    refloc = temprefstart

    for c in cigar:
        if c in numstr:
            number += c
        else:
            number = int(number)
            if readstart == False and c in ['M', 'I', '=', 'X']:
                readstart = readloc

            if readstart != False and c in ['H', 'S']:
                readend = readloc
                refend = refloc
                break
            if c in ['M', 'I', 'S', '=', 'X']:
                readloc += number
            if c in ['M', 'D', 'N', '=', 'X']:
                refloc += number
            number = ""
    return temprefstart, refend, readstart, readend


def find_clip_pos(position, cigar, chrom, read, min_size, max_size):
    resultsv = []
    rawsalist = read.get_tag("SA").split(";")[:-1]
    code = decode_flag(read.flag)
    for sa in rawsalist:
        sainfo = sa.split(",")
        tempchrom, temprefstart, strand, cigar = sainfo[0], int(sainfo[1]), sainfo[2], sainfo[3]
        if tempchrom != chrom:
            continue
        if (strand == "-" and code % 2 == 0) or (strand == "+" and code % 2 != 0):
            refstart1, refend1, readstart1, readend1 = read.reference_start, read.reference_end, read.query_alignment_start, read.query_alignment_end
            refstart2, refend2, readstart2, readend2 = c_pos(cigar, temprefstart)

            a = readend1 - readstart2
            b = refend1 - refstart2
            if abs(a) <= 500:
                if abs(b - a) < 30:
                    continue
                if b - a < 0:
                    deletionlen = abs(b - a)
                    if min_size <= deletionlen <= max_size:
                        resultsv += [[chrom, min(refend1, refstart2), deletionlen, "d-clip"]]
    return resultsv


def findcigarinfo(cigar):
    numbers = '1234567890'
    num = ''
    reflen = 0
    readlen = 0
    leftclip = 0
    rightclip = 0
    for c in cigar:
        if c in numbers:
            num += c
            continue
        if c in "MNP=X":
            readlen += int(num)
            reflen += int(num)
            num = ""
            continue
        if c == "I":
            readlen += int(num)
            num = ""
            continue
        if c == "D":
            reflen += int(num)
            num = ""
            continue
        if c in "SH":
            if readlen == 0:
                leftclip = int(num)
            else:
                rightclip = int(num)
            num = ""
            continue

    return [reflen, [leftclip, readlen, rightclip]]


def finddeletionbysegment(segments, min_size, max_size):
    # if len([c for c in segments if int(c[1]) <= 16]) == 0:
    #     return []
    if len(segments) <= 1:
        return []
    deletionset = []
    primary = [c for c in segments if len(c) == 8][0]
    segments = [c for c in segments if c != primary]
    readseqence = primary[7]
    chrom = primary[2]
    priflag = (int(primary[1]) % 32) > 15
    samedirchr = []
    # 找出跟主比对 比对相同的方向和染色体
    for c in segments:
        ch = c[2]
        f = (int(c[1]) % 32) > 15
        readlen = c[5][1]
        if readlen < 300:
            continue
        if ch != chrom:
            continue
        elif f != priflag:
            continue
        else:
            samedirchr += [c]

    for c in samedirchr:
        if c[3] > primary[3]:
            leftread = primary
            rightread = c
        elif c[3] < primary[3]:
            leftread = c
            rightread = primary
        else:
            continue
        leftinfo = leftread[5]
        rightinfo = rightread[5]
        overlapmap = leftinfo[0] + leftinfo[1] - rightinfo[0]
        window_max = 1500
        if -200 < overlapmap < window_max:
            deletionsize = rightread[3] - leftread[4] + overlapmap
            if min_size <= deletionsize <= max_size:
                deletionset += [
                    [chrom, max(0, leftread[4] - max(0, overlapmap)), deletionsize, 'd-clip']
                ]

    return deletionset


def async_find_sv_cigar_softclip(chrom, chromlength, bam_path, outpath):
    # print("{}:{}:{}:{}".format(chrom, chromlength, bam_path, out_path))
    print("==============={} find deletion is start===============".format(chrom))
    bamfile = pysam.AlignmentFile(bam_path, "r")
    allreads = bamfile.fetch(contig=chrom)
    min_size = 50
    max_size = 100000
    cigarsv = []
    clipsv = []
    segmentreads = {}
    for read in allreads:
        if read.is_unmapped or read.is_secondary:
            continue

        position = read.reference_start + 1
        readname = read.query_name
        flag = read.flag
        mapq = read.mapping_quality

        cigars = read.cigarstring
        cigarinfo = find_cigar_pos(position, cigars, chrom)
        cigarsv += cigarinfo[0]

        if read.has_tag("SA"):
            clipsv += find_clip_pos(position, cigars, chrom, read, min_size, max_size)

    sv = cigarsv + clipsv
    # sv = np.array(sv)

    print("{}号染色体通过cigar查找到{}个deletion".format(chrom, len(cigarsv)))
    print("{}号染色体通过clip查找到{}个deletion".format(chrom, len(clipsv)))

    cigarsv = pd.DataFrame(cigarsv, columns=["chrom", "position", "length", "type"])
    cigarsv = cigarsv.sort_values(by="position", ascending=True, axis=0)
    clipsv = pd.DataFrame(clipsv, columns=["chrom", "position", "length", "type"])
    clipsv = clipsv.sort_values(by="position", ascending=True, axis=0)
    sv = pd.DataFrame(sv, columns=["chrom", "position", "length", "type"])
    sv = sv.sort_values(by="position", ascending=True, axis=0)

    cigarsvpath = os.path.join(outpath, "cigarsv_{}.csv".format(chrom))
    clipsvpath = os.path.join(outpath, "clipsv_{}.csv".format(chrom))
    svpath = os.path.join(outpath, "sv_{}.csv".format(chrom))
    cigarsv.to_csv(cigarsvpath, header=True, index=False)
    clipsv.to_csv(clipsvpath, header=True, index=False)
    sv.to_csv(svpath, header=True, index=False)
    print("==============={} find deletion is over~~~===============".format(chrom))


def print_error(value):
    print("error:", value)
    # sys.exit(0)
    os._exit(0)


def find_sv(bam_path, max_work, outpath, chromosomeslength):
    pool = mp.Pool(processes=max_work)
    svpath = os.path.join(outpath, "svtemp")
    if not os.path.exists(svpath):
        os.mkdir(svpath)
    for chrom in chromosomeslength:
        ppath = os.path.join(svpath, "{}_temp".format(chrom))
        if not os.path.exists(ppath):
            os.mkdir(ppath)
        pool.apply_async(async_find_sv_cigar_softclip, args=(chrom, chromosomeslength[chrom], bam_path, ppath),
                         error_callback=print_error)

    pool.close()
    pool.join()


def read_sv_from_file(deletionsv, chromosomeslength, outpath):
    for chrom in chromosomeslength:
        svpath = os.path.join(outpath, "svtemp/{}_temp/sv_{}.csv".format(chrom, chrom))
        sv = pd.read_csv(svpath)
        sv = np.array(sv)
        for i in sv:
            deletionsv[chrom] += [list(i)]


def sortbyposlength(c):
    return [c[1], c[2]]


def merge_deletion(candi, minsupport):
    """
    对大的deletion进行聚类，并返回该簇所代表的SV
    :param candi:
    :param minsupport:
    :return:
    """
    candi.sort(key=lambda x: x[2])
    index = 0
    rawsv = [[c] for c in candi]
    while index + 1 < len(rawsv):
        average_len1 = 0
        cluster1 = rawsv[index]
        for i in range(len(cluster1)):
            average_len1 += cluster1[i][2]
        average_len1 = int(average_len1 / len(cluster1))

        average_len2 = 0
        cluster2 = rawsv[index + 1]
        for i in range(len(cluster2)):
            average_len2 += cluster2[i][2]
        average_len2 = int(average_len2 / len(cluster2))

        rate = abs(average_len2 - average_len1) / (max(average_len1, average_len2) * 1.0)
        if rate <= 0.2:
            rawsv[index] += rawsv[index + 1]
            rawsv.pop(index + 1)
        else:
            index = index + 1

    newsv = []
    for cluster in rawsv:
        if len(cluster) >= minsupport:
            average_length = 0
            for deletion in cluster:
                average_length += deletion[2]
            average_length = int(average_length / len(cluster))
            index = 0
            minvalue = 1000000
            for i in range(len(cluster)):
                deletion = cluster[i]
                length = deletion[2]
                interval = abs(length - average_length)
                if interval < minvalue:
                    index = i
                    minvalue = interval
            nearaverage = cluster[index]
            tempsv = [nearaverage[0], nearaverage[1], nearaverage[2], nearaverage[3], len(cluster)]
            newsv += [[tempsv, cluster]]

    return newsv


def calculate_cluster_average_lens(cluster):
    average_lens = 0
    for deletion in cluster:
        average_lens += deletion[2]
    average_lens = int(average_lens / len(cluster))
    return average_lens


def find_deletion_from_cluster(deletioncluster):
    average_lens = 0
    for deletion in deletioncluster:
        lens = deletion[2]
        average_lens += lens
    average_lens = average_lens / len(deletioncluster)
    mindis = 10000
    index = 0
    for i in range(len(deletioncluster)):
        deletion = deletioncluster[i]
        lens = deletion[2]
        dis = abs(lens - average_lens)
        if dis < mindis:
            mindis = dis
            index = i

    nearaverage = deletioncluster[index]
    svtemp = [nearaverage[0], nearaverage[1], nearaverage[2], nearaverage[3], len(deletioncluster)]
    return [svtemp]


def hierarchical_cluster(deletion, minsupport):
    """
    分层聚类：将每个deletion看作一个簇，每次找到差异最小的两个簇进行合并，直到只剩下两个簇，
    然后再判断这两个簇是否可以合并
    :param deletionset:
    :param minsupport:
    :return:
    """
    deletion.sort(key=lambda x: x[2])
    rawdeletion = []
    for i in deletion:
        rawdeletion += [[i]]

    num = 0
    while len(rawdeletion) > 2:
        ratematrix = []
        for i in range(len(rawdeletion)):
            ratematrix.append([0] * len(rawdeletion))
        minrate = 1
        index1 = 0
        index2 = 0
        for i in range(len(rawdeletion)):
            for j in range(len(rawdeletion)):
                if i != j:
                    cluster1 = rawdeletion[i]
                    average_len1 = calculate_cluster_average_lens(cluster1)
                    cluster2 = rawdeletion[j]
                    average_len2 = calculate_cluster_average_lens(cluster2)
                    rate = abs(average_len2 - average_len1) / (max(average_len1, average_len2) * 1.0)
                    ratematrix[i][j] = rate
                    if rate <= minrate:
                        index1 = i
                        index2 = j
                        minrate = rate

        rawdeletion[index1] += rawdeletion[index2]
        rawdeletion.pop(index2)

    # print(len(rawdeletion))

    cluster1 = rawdeletion[0]
    cluster2 = rawdeletion[1]
    # print("cluster1:", cluster1)
    # print("cluster2:", cluster2)
    average_len1 = calculate_cluster_average_lens(cluster1)
    average_len2 = calculate_cluster_average_lens(cluster2)
    # print("average_len1:", average_len1)
    # print("average_len2:", average_len2)

    if abs(average_len1 - average_len2) < 20 + 0.01 * max(average_len1, average_len2):
        rawdeletion[0] += rawdeletion[1]
        rawdeletion.pop(1)


    fdeletion = []
    if len(rawdeletion) == 2:
        minsup = minsupport // 2
        for deletioncluster in rawdeletion:
            if len(deletioncluster) >= minsup:
                fdeletion += find_deletion_from_cluster(deletioncluster)

    elif len(rawdeletion) == 1:
        deletioncluster = rawdeletion[0]
        fdeletion += find_deletion_from_cluster(deletioncluster)

    return fdeletion


def merge_small_deletion(deletionset, minsupport):
    """
    分成两种情况，若大于minsupport，则对其进行聚类，如果大于1.5倍minsupport，则分层聚类；
    否则直接聚类
    若小于minsupport，则不进行聚类，直接返回[]即可
    :param deletionset:
    :param minsupport:
    :return: 代表该簇的deletion
    """
    deletionset.sort(key=lambda x: x[2])

    clustersv = []
    if len(deletionset) >= minsupport:
        if len(deletionset) >= 1.5 * minsupport:
            upper = deletionset[len(deletionset) * 3 // 4][2]
            lower = deletionset[len(deletionset) // 4][2]
            if upper > 1.75 * lower and upper - lower > 50:
                # 分层聚类
                clustersv += hierarchical_cluster(deletionset, minsupport)
            else:
                result = merge_deletion(deletionset, minsupport)
                for elem in result:
                    clustersv += [elem[0]]

        else:
            result = merge_deletion(deletionset, minsupport)
            for elem in result:
                clustersv += [elem[0]]
    else:
        # 该deletion簇中个数小于最小支持读数，则不进行聚类处理
        pass

    # print(clustersv)

    return clustersv


def positionandsvlen(a):
    return [a[1], a[2]]


def async_cluster_deletion(allsv, outpath, minsupport, maxdepth, chrom, chromlength):
    """
    分成大的deletion和小的deletion来分别进行聚类
    对于大的deletion（svlen > 2000）来说，首先设置一个长度为1500的窗口，然后滑动窗口将该窗口内的大deletion放入到一个
    列表中，判断该列表中deletion的个数是否大于min_support，如果大于min_support，则对该列表中的deletion
    进行聚类，聚类过程如下：
    1、将列表中deletion按照长度进行排序，并将每个deletion视为一个簇；
    2、对于两个deletion d1和d2来说，若它们的长度差异小于20%，则将它们合并为一个簇
    3、重复该合并过程，直到最后没有可以合并的簇为止
    对于小的deletion(svlen <= 3000)来说，首先按照其在参考基因组的位置和svlen进行排序，然后使用基于密度的
    方法进行聚类，在聚类之后会形成多个簇，判断其四分位点数，以此来判断其是否是杂合的
    :param allsv:
    :param outpath:
    :param minsupport:
    :param maxdepth:
    :param chrom:
    :param chromlength:
    :return:
    """

    path = os.path.join(outpath, "cluster_temp")
    if not os.path.exists(path):
        os.mkdir(path)
    chrompath = os.path.join(path, "{}_cluster".format(chrom))
    if not os.path.exists(chrompath):
        os.mkdir(chrompath)

    clusterpath = os.path.join(chrompath, "{}_larger_cluster.txt".format(chrom))
    cpfile = open(clusterpath, "w")
    largeSVpath = os.path.join(chrompath, "{}_largerSV.txt".format(chrom))
    lspfile = open(largeSVpath, "w")

    # largeSV
    largesv = [c for c in allsv if c[2] > 2000]
    window = 1500
    largesv.sort(key=sortbyposlength)
    start = 0
    largerdel = []
    candi = []
    for event in largesv:
        pos = event[1]
        if pos <= start + window:
            candi += [event]
            continue
        if len(candi) >= minsupport:
            largerdel += merge_deletion(candi, minsupport)

        candi = [event]
        start = pos

    if len(candi) >= minsupport:
        largerdel += merge_deletion(candi, minsupport)
        candi = []

    if len(largerdel) != 0:
        for deletion in largerdel:
            cpfile.write(str(deletion) + "\n")
            lspfile.write(str(deletion[0]) + "\n")
    cpfile.close()
    lspfile.close()


    # small deletion
    smallSV = [c for c in allsv if c[2] <= 3000]
    genomeposition = [0] * chromlength
    # 计算出每个SV区域位点的支持读数覆盖
    for c in smallSV:
        start = c[1]
        end = start + c[2]
        svregion = genomeposition[start - 1: end - 1]
        new = [value + 1 for value in svregion]
        genomeposition[start - 1: end - 1] = new

    svregion = []
    inpeak = False
    threshold = 3
    maxdep = 0
    localdep = []

    # 遍历每个位点，去查找带有peak的区域，并记录其区域
    for i in range(len(genomeposition)):
        # print("{}:{}".format(chrom, i))
        if inpeak:
            if genomeposition[i] > max(maxdep / 10.0, threshold):
                localdep += [genomeposition[i]]
                if genomeposition[i] > maxdep:
                    maxdep = genomeposition[i]
            else:
                inpeak = False
                end = i
                if maxdep <= maxdepth:
                    peakpos = localdep.index(maxdep)
                    peakkleftsize = 0
                    for k in range(peakpos):
                        if localdep[peakpos - k - 1] >= maxdep / 10.0:
                            peakkleftsize += 1
                        else:
                            break
                    svregion += [
                        (start + peakpos - peakkleftsize, end, maxdep)
                    ]
                start = 0
                end = 0
                maxdep = 0
        else:
            if genomeposition[i] >= threshold:
                inpeak = True
                localdep = [genomeposition[i]]
                start = i
                maxdep = genomeposition[i]

    # 用找到的peak区域做键，将该区域内包含的deletion放入到字典中
    svregion = [s for s in svregion if s[2] < maxdepth]
    alldeletioninfor = {}
    for s in svregion:
        alldeletioninfor[s] = []

    for c in smallSV:
        start = c[1]
        end = start + c[2]
        for d in svregion:
            sregion = d[0]
            eregion = d[1]
            if min(end, eregion) - max(start, sregion) > 0:
                alldeletioninfor[d] += [c]

    # 将small deletion输出查看
    smallpeakpath = os.path.join(chrompath, "{}_peak_region.txt".format(chrom))
    sppfile = open(smallpeakpath, "w")
    for c in alldeletioninfor:
        sppfile.write(str(c) + "\t" + str(alldeletioninfor[c]) + "\n")
    sppfile.close()

    # small deletion cluster
    smalldeletion = []
    for c in svregion:
        svinfo = alldeletioninfor[c]
        smalldeletion += merge_small_deletion(svinfo, minsupport)

    """
    此时已经找到了大于2000bp（大的deletion）和小于3000bp（小的deletion）的deletion，
    此时需要将其合并起来重叠的
    """
    newSV = []
    for lagers in largerdel:
        lager = lagers[0]
        ovlapflag = 0
        for small in smalldeletion:
            large_start = lager[1]
            large_svlen = lager[2]
            large_end = large_start + large_svlen
            small_start = small[1]
            small_svlen = small[2]
            small_end = small_start + small_svlen
            if min(large_end, small_end) - max(large_start, small_start) > 0 and 0.8 * small_svlen <= large_svlen <= small_svlen / 0.8:
                ovlapflag = 1
                break
        if ovlapflag == 0:
            newSV += [lager]

    newSV += smalldeletion
    newSV.sort(key=positionandsvlen)
    # 将最后结果写出来检查
    finaldeletionpath = os.path.join(chrompath, "{}_final_deletion.txt".format(chrom))
    fdpfile = open(finaldeletionpath, "w")
    for deletion in newSV:
        fdpfile.write(str(deletion) + "\n")
    fdpfile.close()

    print("{} cluster is over ~~~".format(chrom))


def cluster_deletions(deletionsv, outpath, min_support, highcov, maxwork, chromosomeslength):
    max_work = maxwork // 2
    pool = mp.Pool(processes=max_work)

    for chrom in chromosomeslength:
        pool.apply_async(async_cluster_deletion,
                         args=(deletionsv[chrom], outpath, min_support, highcov, chrom, chromosomeslength[chrom]),
                         error_callback=print_error)

    pool.close()
    pool.join()


def read_deletion_from_finalDeletion(deletiondict, chromosomeslength, outpath):
    for chrom in deletiondict:
        path = os.path.join(outpath, "cluster_temp/{}_cluster/{}_final_deletion.txt".format(chrom, chrom))
        rfile = open(path, "r")
        while True:
            line = rfile.readline().rstrip()
            if not line:
                break
            line = list(eval(line))
            deletiondict[chrom] += [line]
        rfile.close()


def filter_small_copy_deletion(deletiondict):
    chroms = [ch for ch in deletiondict]
    for chromm in chroms:
        deletionset = deletiondict[chromm]
        newDeletion = []
        # filter small
        for deletion in deletionset:
            svlen = deletion[2]
            if 40 <= svlen < 50:
                chrom = deletion[0]
                start = deletion[1]
                svlen = 50
                types = deletion[3]
                supp = deletion[4]
                newDeletion += [[chrom, start, svlen, types, supp]]
            elif svlen < 40:
                continue
            elif svlen >= 50:
                newDeletion += [deletion]
        # print("{}染色体原先有：{}".format(chrom, len(deletionset)))
        # filter
        findex = 0
        bindex = findex + 1
        while bindex < len(newDeletion):
            fdeletion = newDeletion[findex]
            bdeletion = newDeletion[bindex]
            fstart = fdeletion[1]
            fsvlen = fdeletion[2]
            fend = fstart + fsvlen
            bstart = bdeletion[1]
            bsvlen = bdeletion[2]
            if abs(bstart - fstart) <= 50 and bsvlen * 0.8 < fsvlen < bsvlen * 1.2:
                index = 0
                if fsvlen >= bsvlen:
                    index = bindex
                else:
                    index = findex
                newDeletion.pop(index)
            else:
                findex = bindex
                bindex = findex + 1
        deletiondict[chromm] = newDeletion
        # print("{}号染色体过滤后有：{}".format(chrom, len(newDeletion)))
        print("{} filter is over ~~~~".format(chromm))


def generate_vcf(inslist, contiglength, out_vcfpath):
    head = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">\n"""
    body = ''
    for contig in contiglength:
        body += "##contig=<ID=" + contig + ",length=" + str(contiglength[contig]) + ">\n"
    tail = """##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of SV:DEL=Deletion, INS=Insertion">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t.\n"""
    vcfinfo = head + body + tail
    for rec in inslist:
        contig = rec[0]
        geno = rec[-1]
        recinfo = 'SVLEN=' + str(rec[2]) + ';SVTYPE=' + 'DEL' + ';END=' + str(rec[3]) + ';' + '\tGT\t' + str(
            geno) + '\n'
        vcfinfo += (str(contig) + '\t' + str(
            rec[1]) + '\t' + '.' + '\t' + '.' + '\t' + 'A' + '\t' + '.' + '\t' + 'PASS' + '\t' + recinfo)
    f = open(out_vcfpath, 'w')
    f.write(vcfinfo)
    f.close()


def finalsort(a):
    return [a[0], a[1], a[2]]


def main(bampath, maxwork, outpath):
    """
       1、首先计算该文件的覆盖深度，在后面进行聚类后筛选需要使用；
       2、分别通过读数内比对和读数间比对进行deletion查找，并将每条染色体查找到的候选位点排序后写入到文件中
       :param bampath:
       :param maxwork:
       :param outpath:
       :return:
       """
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    print("------开始计算覆盖度------")
    depth = find_coverage(bampath, maxwork, outpath)
    # min_support = round(depth / 10) + 0

    highcov = max(50, round(depth * 2 + 20))
    min_support = 5
    # highcov = 50
    print("测序深度约为：{}".format(depth))
    print("最小支持读数为：{}".format(min_support))
    print("最大覆盖度为：{}".format(highcov))

    # 这里仅仅查找前22条染色体的长度即可

    chromosomeslength = findlength(bampath)
    print("------开始遍历cigar查找deletion------")
    find_sv(bampath, maxwork, outpath, chromosomeslength)

    # 此时已经从cigar和softclip查找到sv，下面开始聚类
    # 首先把写入文件的SV按照染色体读进内存，存入一个字典，字典的键就是染色体
    # 然后对其进行聚类
    deletionsv = {}
    for i in chromosomeslength:
        deletionsv[i] = []

    read_sv_from_file(deletionsv, chromosomeslength, outpath)
    cluster_deletions(deletionsv, outpath, min_support, highcov, maxwork, chromosomeslength)

    """
    将最后聚类形成的deletion读到内存中，然后过滤掉长度小于45的和重复的deletion
    """
    deletiondict = {}
    for i in chromosomeslength:
        if "X" in i or "Y" in i:
            continue
        deletiondict[i] = []
    read_deletion_from_finalDeletion(deletiondict, chromosomeslength, outpath)
    filter_small_copy_deletion(deletiondict)

    """
    将结果写入到vcf文件中
    """
    finalDeletion = []
    for chrom in deletiondict:
        finalDeletion += deletiondict[chrom]
        # print("{}:{}, {}".format(chrom, len(deletiondict[chrom]), type(chrom)))
    finalDeletion.sort(key=finalsort)

    data = []
    for deletion in finalDeletion:
        ch = deletion[0]
        start = deletion[1]
        svlen = deletion[2]
        end = start + svlen
        data += [[ch, start, -svlen, end, "."]]

    vcfpath = os.path.join(outpath, "deletion10x_0.2+0_p3.vcf")
    generate_vcf(data, chromosomeslength, vcfpath)


# clr
# bampath = "/public_data/SV/long_read/HG002_PB_70x_RG_HP10XtrioRTG_donwsample5.bam"
# bampath = "/public_data/SV/long_read/HG002_PB_70x_RG_HP10XtrioRTG_downsample12.bam"
# bampath = "/public_data/SV/long_read/HG002_PB_70x_RG_HP10XtrioRTG_downsample25.bam"
# bampath = "/public_data/SV/long_read/HG002_PB_70x_RG_HP10XtrioRTG.bam"
# bampath = "/public_data/SV/long_read/HG002_PB_70x_RG_HP10XtrioRTG_downsample50.bam"


# ccs
bampath = "/public_data/SV/long_read/HG002.Sequel.15kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.MD.bam"
# bampath = "/public_data/SV/long_read/HG002.Sequel.15kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio_downsample33.MD.bam"
# bampath = "/public_data/SV/long_read/HG002.Sequel.15kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio_downsample16.MD.bam"


# ont
# bampath = "/public_data/SV/long_read/HG002_GRCh37_ONT-UL_UCSC_20200508.phased.MD.bam"
# bampath = "/public_data/SV/long_read/HG002_GRCh37_ONT-UL_UCSC_20200508.phased_downsample40.MD.bam"
# bampath = "/public_data/SV/long_read/HG002_GRCh37_ONT-UL_UCSC_20200508.phased_downsample20.MD.bam"
# bampath = "/public_data/SV/long_read/HG002_GRCh37_ONT-UL_UCSC_20200508.phased_downsample10.MD.bam"

max_work = 20
outpath = "/home/yuyanan/AssemblySv/ccs_test/28x"

main(bampath, max_work, outpath)

