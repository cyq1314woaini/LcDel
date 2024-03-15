# coding=utf-8
"""
将cutesv运行的NA19240数据集结果转化为一般代表deletion的vcf文件，用来评价
"""

NA19240gold_deletion = "E:/SV的论文/assemblySV/NA19240/NA19240goldchr1-22del.vcf"
invcf = "E:/SV的论文/assemblySV/NA19240/cutesv.na19240.5x.s2.vcf"
outvcf = "E:/SV的论文/assemblySV/NA19240/cutesv.na19240.5x.deletion.s2.vcf"


nafile = open(NA19240gold_deletion, "r")
outfile = open(outvcf, "w")
while True:
    line = nafile.readline().rstrip()
    if not line:
        break
    if line[0] == "#":
        outfile.write(line + "\n")
        continue
nafile.close()

infile = open(invcf, "r")
while True:
    line = infile.readline().rstrip()
    if not line:
        break
    if line[0] == "#":
        continue
    if line[:3] == "chr":
        lineinfo = line.split("\t")
        typeinfo = lineinfo[2]
        types = typeinfo.split(".")[1]
        if types == "DEL":
            info = lineinfo[7]
            infos = info.split(";")
            svlen = int(infos[2][6: ])
            if abs(svlen) >= 50:
                outfile.write(line[3:] + "\n")

infile.close()
outfile.close()
