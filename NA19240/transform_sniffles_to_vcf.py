# coding=utf-8
"""
将sniffles运行的NA19240数据集结果转化为一般代表deletion的vcf文件，用来评价
"""

NA19240gold_deletion = "E:/SV的论文/assemblySV/NA19240/NA19240goldchr1-22del.vcf"
invcf = "E:/SV的论文/assemblySV/NA19240/sniffles.na19240.5x.vcf"
outvcf = "E:/SV的论文/assemblySV/NA19240/sniffles.na19240.5x.deletion.s2.vcf"
supportt = 2

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
    line_info = line.split("\t")
    typeinfo = line_info[2]
    types = typeinfo.split(".")[1]
    if types == "DEL":
        info = line_info[7]
        infos = info.split(";")
        svlen = int(infos[2][6:])
        support = int(infos[4][8:])
        # print(svlen, support)
        if abs(svlen) >= 50 and support >= supportt:
            outfile.write(line[3:] + "\n")

outfile.close()
infile.close()

