# coding=utf-8

import pysam

"""
第二个点：使用组装的方法去检测SV
"""

# print("hello world~~~")
# print("hello world2~~~")
# print("hello world3~~~")
#
# for i in range(1000000):
#     print(i)


# a = [1, 2, 3, 4]
# print(sum(a))

import pysam


def decode_flag(Flag):
    signal = {1 << 2: 0, 1 >> 1: 1, 1 << 4: 2, 1 << 11: 3, 1 << 4 | 1 << 11: 4}
    return signal[Flag] if (Flag in signal) else 0


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


path = "/public_data/SV/long_read/HG002_PB_70x_RG_HP10XtrioRTG_downsample25.bam"
bamfile = pysam.AlignmentFile(path, "r")
allreads = bamfile.fetch(contig="1")

chrom = "1"

for read in allreads:
    if read.has_tag("SA"):
        code = decode_flag(read.flag)  # 当code为偶数时，代表原读数反向比对；code为奇数时代表原读数正向比对
        rawsalist = read.get_tag("SA").split(";")[:-1]
        for sa in rawsalist:
            sainfo = sa.split(",")
            tempchrom, temprefstart, strand, cigar = sainfo[0], int(sainfo[1]), sainfo[2], sainfo[3]
            if tempchrom != chrom:
                continue

            if (strand == '-' and code % 2 == 0) or (strand == '+' and code % 2 != 0):
                refstart1, refend1, readstart1, readend1 = read.reference_start, read.reference_end, read.query_alignment_start, read.query_alignment_end
                # print(read)
                # print(refstart1, refend1, readstart1, readend1)
                refstart2, refend2, readstart2, readend2 = c_pos(cigar, temprefstart)
                a = readend1 - readstart2
                b = refend1 - refstart2
                if abs(a) < 500:
                    if abs(b-a) < 30:
                        continue
                    if b - a < 0:
                        data = [min(refend1, refstart2), min(refend1, refstart2) + abs(b-a), abs(b-a)]
                        if abs(b-a) <= 5000:
                            print(data)
                            print(read)
                            break




