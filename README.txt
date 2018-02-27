############################################################################
20180226
####v1.0.0
-the pipeline is now available in the 'lite' version i.e. with hardcoded software choices
-read mapping is done with bwa-meth, methylation extraction with MethylDackel, stats with limma, DMR calling with metilene
-functionality of trimReads has been modified: defaults to None (no read trimming). Other choices are: 'auto' - illumina adapter trimming with auto-detected 5'end hard-trimming or 'user' with user-supplied arguments to cutadapt.

