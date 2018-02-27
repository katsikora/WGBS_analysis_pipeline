#this module uses available FastQC reports to determine the 5' hard trimming length for R1 and R2
#this is taken as the maximum of the rowwise differences for nucleotide frequency at position in the read
#standard Illumina adapters are removed
#a post-trim FastQC is run and the user is strongly encouraged to examine the results

import os
import argparse
import sys
import zipfile
import re
import subprocess
import pandas
import numpy
import commands
import math
import time
from collections import OrderedDict
from ruffus.drmaa_wrapper import run_job, error_drmaa_job
import logging
from ruffus.proxy_logger import *
    
#####################################  DEFINITIONS ############################################################
def calc_cutThd (zipL,fqin,fqout,shared_logger,logging_mutex):
    fqcdir=fqin
    fqcout=fqout
    rNcutL=[]
    with logging_mutex:
        for zipi in [ item for item in zipL ]: 
            zf=os.path.basename(zipi)
            if not os.path.exists(os.path.join(fqcout,re.sub('\.zip','',zf))):
                with zipfile.ZipFile(zipi, "r") as z:
                        z.extractall(path=fqcout)
            fqtxt=os.path.join(fqcout,re.sub('\.zip','',zf),'fastqc_data.txt')
            shared_logger.info('Currently processing :'+ fqtxt)
            os.chdir(os.path.join(fqcout,re.sub('\.zip','',zf)))
            subprocess.check_output(['csplit', '-z' , fqtxt , '/>>/','{*}'])
            with open(fqtxt,'r') as file:
                line=file.readline().strip()
            if '0.11.2' in line or '0.11.6' in line:
                NTconTab=pandas.read_table(os.path.join(os.getcwd(), 'xx09'), sep='\t',skiprows=1,header=0,names=['Index','G','A','T','C'],dtype={'Index':'object','G':'float64','A':'float64','T':'float64','C':'float64'},engine='c')
            else:
                NTconTab=pandas.read_table(os.path.join(os.getcwd(), 'xx09'), sep='\t',skiprows=1,header=0,names=['Index','G','A','T','C'],dtype={'Index':'object','G':'float64','A':'float64','T':'float64','C':'float64'},engine='c')
                shared_logger.warning('Check fastqc version')
            difftab=NTconTab.set_index('Index').diff(periods=-1)
            difftabA=difftab.abs()
            maxv=difftabA.idxmax(axis=0)
            maxv=maxv.values.astype(int)
            rNmax=list(difftabA.index)
            rNcut=rNmax[(maxv.max()-1)]
            rNcutL.append(str(rNcut)) ##
            shared_logger.info(NTconTab.head(n=10))
            shared_logger.info(difftab.head(n=10))
            shared_logger.info('Maximal absolute difference per nucleotide :')
            shared_logger.info(difftabA.max(axis=0))
            shared_logger.info('Index of diffmax :')
            shared_logger.info(difftabA.idxmax(axis=0))
            shared_logger.info('Index of the maximal difference :')
            shared_logger.info(maxv.max())
            shared_logger.info('Number of nucleotides for 5prime trimming :' + rNcut)
            os.getcwd()
    zipLre=[ re.sub('_fastqc.zip','.fastq.gz',x ) for x in zipL ]    
    cutThdRes=OrderedDict(zip(zipLre, rNcutL))
    ctr1=filter(lambda x:'_R1.fastq.gz' in x, cutThdRes.keys())
    ctr2=filter(lambda x:'_R2.fastq.gz' in x, cutThdRes.keys())
    cutThdRes_R1=[ cutThdRes[x] for x in ctr1 ]
    cutThdRes_R2=[ cutThdRes[x] for x in ctr2 ]
    cutThdL=zip(cutThdRes_R1,cutThdRes_R2)
    return cutThdL


def cut_reads_auto(INfile1,INfile2,OUTfile1,OUTfile2,cutThdR1,cutThdR2,cutpath,my_session,cutout,logobject):
    read_root=re.sub('_R1.fastq.gz','',os.path.basename(INfile1))
    bshcmd=os.path.join(cutpath,'cutadapt') + ' -a AGATCGGAAGAGC -A AGATCGGAAGAGC --minimum-length 30  -n 5 -u ' + cutThdR1 + ' -U ' + cutThdR2 + ' -o ' + OUTfile1 + ' -p ' + OUTfile2 + ' ' + INfile1 + ' ' + INfile2
    with open(os.path.join(cutout,"logs","%s.trim_reads.out" % read_root),'w+') as stdoutF, open(os.path.join(cutout,"logs","%s.trim_reads.err" % read_root),'w+') as stderrF:
        try:
        
            stdout_res, stderr_res  = run_job(cmd_str           = bshcmd,
                                      job_name          = 'cut_reads',
                                      logger            = logobject,
                                      drmaa_session     = my_session,
                                      run_locally       = False,
                                      working_directory = os.getcwd(),
                                      job_other_options = '-p bioinfo')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        except error_drmaa_job as err:
            logobject.error("Cut_reads error: %s" % err)
            raise
    return

def cut_reads_user(INfile1,INfile2,OUTfile1,OUTfile2,cutpath,my_session,cutout,logobject, args):
    read_root = os.path.basename(INfile1)[:-12]
    adapterSeq = "AGATCGGAAGAGC"
    if args.nextera:
        adapterSeq = "CTGTCTCTTATA"
    bshcmd = "{} -a {} -A {} -q {} -m 30 -j {} {} -o {} -p {} {} {}".format(os.path.join(cutpath, 'cutadapt'),
                                                                            adapterSeq,
                                                                            adapterSeq,
                                                                            args.trimThreshold,
                                                                            args.trimOtherArgs,
                                                                            args.nthreads,
                                                                            OUTfile1,
                                                                            OUTfile2,
                                                                            INfile1,
                                                                            INfile2)
    with open(os.path.join(cutout,"logs","%s.trim_reads.out" % read_root),'w+') as stdoutF, open(os.path.join(cutout,"logs","%s.trim_reads.err" % read_root),'w+') as stderrF:
        try:
            stdout_res, stderr_res = run_job(cmd_str = bshcmd,
                                             job_name = 'cut_reads',
                                             logger = logobject,
                                             drmaa_session = my_session,
                                             run_locally = False,
                                             working_directory = os.getcwd(),
                                             job_other_options = '-p bioinfo --mincpus={}'.format(nThreads))
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        except error_drmaa_job as err:
            logobject.error("Cut_reads error: %s" % err)
            raise
    return

    
def post_trim_fqc(INfile1,INfile2,fqcout,FQCpath,my_session,logobject):
    read_root=re.sub('_R1.fastq.gz','',os.path.basename(INfile1))
    bshcmd=os.path.join(FQCpath,'fastqc ')+' --outdir ' + fqcout + ' -t 8 '+ INfile1 + ' ' + INfile2
    with open(os.path.join(fqcout,"logs","%s.post_fqc.out" % read_root),'w+') as stdoutF, open(os.path.join(fqcout,"logs","%s.post_fqc.err" % read_root),'w+') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str           = bshcmd,
                                          job_name          = 'post_fqc',
                                          logger            = logobject,
                                          drmaa_session     = my_session,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo --mincpus=8')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except error_drmaa_job as err:
            logobject.error("Post_trim_fastqc error: %s" % err)
            raise
    return
