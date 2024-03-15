# coding=utf-8

"""
将黄金标准集划分为0-200， 200-500，500-1000，1000-2000，2000-
"""
import os

inpath = "/home/yuyanan/AssemblySv/clr_test/70x/split_interval/gold/hg002.gold.del.1-22.vcf"
outpath = "/home/yuyanan/AssemblySv/clr_test/70x/split_interval/gold"

outpath1 = os.path.join(outpath, "hg002.gold.del.1-22.0-200.vcf")
outpath2 = os.path.join(outpath, "hg002.gold.del.1-22.200-500.vcf")
outpath3 = os.path.join(outpath, "hg002.gold.del.1-22.500-1000.vcf")
outpath4 = os.path.join(outpath, "hg002.gold.del.1-22.1000-2000.vcf")
outpath5 = os.path.join(outpath, "hg002.gold.del.1-22.2000+.vcf")

infile = open(inpath, "r")
headStr = ""
while True:
    line = infile.readline().rstrip()
    if not line:
        break
    if line[0] == "#":
        headStr += line + "\n"
    else:
        break

outfile1 = open(outpath1, "w")
outfile2 = open(outpath2, "w")
outfile3 = open(outpath3, "w")
outfile4 = open(outpath4, "w")
outfile5 = open(outpath5, "w")

outfile1.write(headStr)
outfile2.write(headStr)
outfile3.write(headStr)
outfile4.write(headStr)
outfile5.write(headStr)

while True:
    lines = infile.readline().rstrip()
    if not lines:
        break
    if lines[0] != "#":
        line = lines.split("\t")
        info = line[7]
        infos = info.split(";")[20]
        svlen = abs(int(infos[6:]))
        if svlen <= 200:
            outfile1.write(lines + "\n")
        elif 200 < svlen <= 500:
            outfile2.write(lines + "\n")
        elif 500 < svlen <= 1000:
            outfile3.write(lines + "\n")
        elif 1000 < svlen <= 2000:
            outfile4.write(lines + "\n")
        elif svlen > 2000:
            outfile5.write(lines + "\n")

infile.close()
outfile1.close()
outfile2.close()
outfile3.close()
outfile4.close()
outfile5.close()