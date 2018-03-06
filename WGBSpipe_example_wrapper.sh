#!/bin/bash

DRMAA_LIBRARY_PATH=/data/manke/repository/scripts/DNA_methylation/drmaa/lib/libdrmaa.so
export PATH=$PATH:/package/bwa-0.7.4/bin
export R_LIBS_USER=/data/manke/repository/scripts/DNA_methylation/Rlibs.3.3.1
source /data/boehm/sikora/miniconda3/bin/activate NGSpy2.7

python /data/manke/repository/scripts/DNA_methylation/WGBS_pipe/v1.1.0/WGBS_pipe_lite.py --readIn /data/processing3/WGBS_pipe_test_IN/reads --ref GRCz10 --wdir /data/processing3/WGBS_pipe_lite_v1.1.0_example_OUT_auto --trimReads 'auto' --fqcIn /data/processing3/WGBS_pipe_test_IN/in_fastqc --batchSize 24 --intList /data/processing3/WGBS_pipe_test_IN/danRer10.cpgIsland.ext.sorted.chr25.noCHR.bed --sampleInfo /data/processing3/WGBS_pipe_test_IN/example_sampleSheet.csv  


source deactivate

echo 'done all'
