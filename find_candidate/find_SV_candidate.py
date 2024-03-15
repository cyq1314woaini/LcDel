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


def find_clip_pos(position, cigar, chrom, read):
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
                    resultsv += [[chrom, min(refend1, refstart2), abs(b - a), "d-clip"]]
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

    if len([c for c in segments if int(c[1]) <= 16]) == 0:
        return []
    if len(segments) <= 1:
        return []
    deletionset = []
    primary = [c for c in segments if int(c[1]) <= 16][0]
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
        if read.is_supplementary:
            # 记录下补充比对的信息
            cigarinfo = [0, 0, 0]
            cigarleft = read.cigar[0]
            cigarright = read.cigar[-1]
            if cigarleft[0] == 4 or cigarleft[0] == 5:
                cigarinfo[0] = cigarleft[1]
            if cigarright[0] == 4 or cigarright[0] == 5:
                cigarinfo[2] = cigarright[1]
            cigarinfo[1] = read.query_alignment_length
            refend = read.reference_end + 1
            readinfo = [readname, flag, chrom, int(position), refend, cigarinfo, mapq, '']
            if readname in segmentreads:
                segmentreads[readname] += [readinfo]
            else:
                segmentreads[readname] = [readinfo]

            print(read.has_tag("SA"))
            print(read)
        else:
            # 此时说明该比对既不是补充比对也不是次比对，因此可以遍历cigar字符串查找，以及判断其是否具有SA补充字段
            cigars = read.cigarstring
            cigarinfo = find_cigar_pos(position, cigars, chrom)
            cigarsv += cigarinfo[0]

            if read.has_tag("SA"):
                # clipsv += find_clip_pos(position, cigars, chrom, read)
                # 此时该读数具有SA补充字段，分别将其比对信息记录下来放入到segmentreads字典中
                refend = position + cigarinfo[1]
                readinfo = [readname, flag, chrom, position, refend, cigarinfo[2], mapq, read.query_sequence]
                # print(readinfo)
                if readname in segmentreads:
                    segmentreads[readname] += [readinfo]
                else:
                    segmentreads[readname] = [readinfo]

                # 紧接着去查找SA补充字段的读数比对信息，并将其记录下来
                sainfo = read.get_tag("SA").split(";")[:-1]
                for c in sainfo:
                    if c.split(",")[0] != chrom:
                        continue
                    if c.split(",")[2] == "+":
                        saflag = 2048
                    else:
                        saflag = 2064
                    sacigarinfo = findcigarinfo(c.split(",")[3])
                    sareadinfo = [readname, saflag, c.split(",")[0], int(c.split(",")[1]), int(c.split(",")[1]) +
                                  sacigarinfo[0], sacigarinfo[1], c.split(",")[4]]
                    segmentreads[readname] += [sareadinfo]

    # 接下来判断clip中是否有代表deletion的
    # clipsv = []
    for readgroup in segmentreads:
        if len(readgroup) < 2 or len(readgroup) > 10:
            continue
        clipsv += finddeletionbysegment(segmentreads[readgroup], min_size, max_size)

    sv = cigarsv + clipsv
    # sv = np.array(sv)
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


def find_sv(bam_path, max_work, out_path):
    # 这里仅仅查找前22条染色体的长度即可

    chromosomeslength = findlength(bam_path)
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


def main(bampath, maxwork, outpath):
    """
       1、首先计算该文件的覆盖深度，在后面进行聚类后筛选需要使用；
       2、分别通过读数内比对和读数间比对进行deletion查找，并将每条染色体查找到的候选位点排序后写入到文件中
       :param bampath:
       :param maxwork:
       :param outpath:
       :return:
       """
    # depth = find_coverage(bampath, maxwork, outpath)
    # min_support = round(depth / 10) + 2
    # highcov = max(50, round(depth * 2 + 20))
    # print("测序深度约为：{}".format(depth))
    # print("最小支持读数为：{}".format(min_support))
    # print("最大覆盖度为：{}".format(highcov))

    find_sv(bampath, maxwork, outpath)


bampath = "/public_data/SV/long_read/HG002_PB_70x_RG_HP10XtrioRTG_downsample50.bam"
max_work = 10
outpath = "/home/yuyanan/AssemblySv"

main(bampath, max_work, outpath)
