#!/bin/bash

DRMAA_LIBRARY_PATH=/data/manke/repository/scripts/DNA_methylation/drmaa/lib/libdrmaa.so
export PATH=$PATH:/package/bwa-0.7.4/bin
export R_LIBS_USER=/data/manke/repository/scripts/DNA_methylation/Rlibs.3.3.1
source /data/boehm/sikora/miniconda3/bin/activate NGSpy2.7

python /data/manke/repository/scripts/DNA_methylation/WGBS_pipe/v0.02/WGBS.pipe.ruffus.py --readIn /data/processing3/WGBS_pipe_test_IN/reads --ref GRCz10 --wdir /data/processing3/WGBS_pipe_example_OUT --trimReads  --fqcIn /data/processing3/WGBS_pipe_test_IN/in_fastqc --batchSize 24 --intList /data/processing3/WGBS_pipe_test_IN/danRer10.cpgIsland.ext.sorted.chr25.noCHR.bed --sampleInfo /data/processing3/WGBS_pipe_test_IN/example_sampleSheet.csv --DMRpg metilene 


source deactivate

echo 'done all'
