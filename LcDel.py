# coding=utf-8
"""
输入一共有四个参数
1、bampath文件的路径
2、max_work最大进程数
3、outpath输出文件
4、min_support最小支持读数阈值

python LcDel.py bampath max_work outpath min_support
"""
import sys


def usage():
    sys.exit(
        '===========================LcDel====================\n'
        'python LcDel.py bampath max_work outpath min_support \n'
        'bampath: input path to bam file\n'
        'max_work: maximum number of processes\n'
        'outpath: the path to the output file\n'
        'min_support: minimum support read threshold\n'
        "======================================================"
    )


result = sys.argv[1:]
# print(result)
# print(len(result))


if len(result) != 4:
    # 此时输入格式不符合要求
    print("The input format does not meet the requirements!!!")
    usage()
else:
    bampath = result[0]
    max_work = int(result[1])
    outpath = result[2]
    min_support = int(result[3])

    print("bampath:", bampath)
    print("max_work:", max_work)
    print("outpath:", outpath)
    print("min_support:", min_support)

    from find_SV_candidate5 import main

    main(bampath, max_work, outpath, min_support)
