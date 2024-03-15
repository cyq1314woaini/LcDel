# coding=utf-8
"""
当bam文件输入时会首先计算该文件的覆盖度，主要步骤如下：
1、计算bam文件中的染色体在reference上的长度
2、计算比对到不同染色体上的读数总长度
3、用读数的总长度/染色体的总长度

输入：排好序的bam文件
输出：覆盖深度
"""

import pysam
import multiprocessing as mp
import os


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
    print("{}染色体的长度为{}".format(chrom, totalMapLength))


def find_coverage(bam, max_work, output):
    bamfile = pysam.AlignmentFile(bam, "r")

    chromosomes = bamfile.references
    chromolength = {}
    for chrom in chromosomes:
        tempchrom = chrom
        if "chr" in chrom:
            tempchrom = tempchrom[3:]
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


bam_path = "/public_data/SV/long_read/HG002_PB_70x_RG_HP10XtrioRTG_downsample25.bam"
max_work = 10
output = "/home/yuyanan/AssemblySv"
depth = find_coverage(bam_path, max_work, output)
print("覆盖深度为：{}".format(depth))
