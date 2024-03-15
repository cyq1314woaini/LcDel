# LcDel
LcDel: Deletion variation detection using clustering based on long read


## Installation
### Requirements
* python 3.9, pysam, multiprocessing, numpy, pandas, pysam, sys, os
### 1. Create a virtual environment
```
#create
conda create -n LcDel python=3.9
#activate
conda activate LcDel
#deactivate
conda deactivate
```

### 2. clone LcDel
* After creating and activating the LcDel virtual environment, download LcDel from github:
```　 
https://github.com/cyq1314woaini/LcDel.git
cd LcDel
```

### 3. Install 
```　
conda activate LcDel
conda install pysam, multiprocessing, numpy, pandas, pysam, sys, os
```

## Usage
```　 
python LcDel.py bampath max_work outpath min_support
bampath: input path to bam file
max_work: maximum number of processes
outpath: the path to the output file
min_support: minimum support read threshold
eg:python LcDel.py /bampath 10 /outpath 10/5/3/2/
``` 


