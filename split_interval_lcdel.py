# coding=utf-8
"""
尝试将不同长度的SV进行划分，划分为50-200，200-500，500-1000，1000-2000，2000-
在输入文件夹下生成五个vcf文件
"""
import os

inpath = "/home/yuyanan/AssemblySv/clr_test/10x/deletion.clr.10x_0.2+s2_p3.vcf"
outpath = "/home/yuyanan/AssemblySv/clr_test/10x/interval"

outpath1 = os.path.join(outpath, "deletion.clr.10x.r0.2_s2_p3.50-200.vcf")
outpath2 = os.path.join(outpath, "deletion.clr.10x.r0.2_s2_p3.200-500.vcf")
outpath3 = os.path.join(outpath, "deletion.clr.10x.r0.2_s2_p3.500-1000.vcf")
outpath4 = os.path.join(outpath, "deletion.clr.10x.r0.2_s2_p3.1000-2000.vcf")
outpath5 = os.path.join(outpath, "deletion.clr.10x.r0.2_s2_p3.2000+.vcf")

outgfile1 = open(outpath1, "w")
outgfile2 = open(outpath2, "w")
outgfile3 = open(outpath3, "w")
outgfile4 = open(outpath4, "w")
outgfile5 = open(outpath5, "w")


infile = open(inpath, "r")
headStr = ""
while True:
    line = infile.readline().rstrip()
    if not line:
        break
    if line[0] == "#":
        headStr += line + "\n"
    else:
        outgfile1.write(headStr)
        outgfile2.write(headStr)
        outgfile3.write(headStr)
        outgfile4.write(headStr)
        outgfile5.write(headStr)
        break
infile.close()

infile1 = open(inpath, "r")
while True:
    lines = infile1.readline().rstrip()
    if not lines:
        break
    if lines[0] != "#":
        line = lines.split("\t")
        info = line[7]
        lens = int(info.split(";")[0][6:])
        lens = abs(lens)
        if 50 <= lens <= 200:
            outgfile1.write(lines + "\n")
        elif 200 < lens <= 500:
            outgfile2.write(lines + "\n")
        elif 500 < lens <= 1000:
            outgfile3.write(lines + "\n")
        elif 1000 < lens <= 2000:
            outgfile4.write(lines + "\n")
        elif 2000 < lens:
            outgfile5.write(lines + "\n")

infile1.close()
outgfile1.close()
outgfile2.close()
outgfile3.close()
outgfile4.close()
outgfile5.close()


