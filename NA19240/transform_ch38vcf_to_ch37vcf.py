# coding=utf-8
"""
由于我们代码跑完NA19240数据集后，生成的vcf文件头文件是针对hg38的并且里面都是带chr，因此需要写代码进行处理
"""

NA19240gold_deletion = "E:/SV的论文/assemblySV/NA19240/NA19240goldchr1-22del.vcf"
invcfpath = "E:/SV的论文/assemblySV/NA19240/deletion.NA19240.40x_r0.2_s5_p3.vcf"
outvcfpath = "E:/SV的论文/assemblySV/NA19240/deletion.NA19240.chr37.40x_r0.2_s5_p3.vcf"

# 首先将标准文件的头文件给读入进去
gfile = open(NA19240gold_deletion, "r")
outfile = open(outvcfpath, 'w')
while True:
    line = gfile.readline().rstrip()
    if not line:
        break
    if line[0] == "#":
        outfile.write(line + "\n")

gfile.close()

infile = open(invcfpath, "r")
while True:
    line = infile.readline().rstrip()
    if not line:
        break
    if line[0] != "#":
        outfile.write(line[3:] + "\n")

infile.close()
outfile.close()