#####IMPORTS####################################################

import os
import re
import logging
from ruffus.drmaa_wrapper import run_job, error_drmaa_job
import subprocess
import string
import numpy as np
import pandas as pd
from vcf import Reader
import collections as cll

#####DEFINITIONS#################################################

def single_CpG_limma(ii,sampleInfo,outdir,Rpath,Rlib,pipev,my_session,logobject):
    Rstat_cmd=os.path.join(Rpath,'Rscript') +' --no-save --no-restore /data/manke/repository/scripts/DNA_methylation/WGBS_pipe/' + pipev + '/WGBSpipe.singleCpGstats.limma.R ' + outdir + ' ' + sampleInfo + ' ' + ii + ' ' + Rlib +' ;sleep 300'
    logobject.info(Rstat_cmd)
    with open(os.path.join(outdir,"logs","singleCpG_stats.out" ),'w') as stdoutF, open(os.path.join(outdir,"logs","singleCpG_stats.err"),'w') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str       = Rstat_cmd,
                                          job_name          = 'sCpG',
                                          logger            = logobject,
                                          drmaa_session     = my_session,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except Exception as err:
            logobject.error("Single CpG stats error: %s" % err)
            raise
        else:
            logobject.info('Single CpG stats calculation complete')
    return

def mCT_get_CpGxInt(ii,outList,refG,bedList,outdir,bedpath,my_session,logobject):
    os.chdir(outdir)
    pozF=os.path.join(outdir,re.sub('.fa*','.poz',os.path.basename(refG)))
    imdF=os.path.join(outdir,re.sub('.poz','.CpG.bed',os.path.basename(pozF)))
    if not os.path.exists(imdF):
        cmd0='grep "+" ' + pozF + ' | awk \'{print $1, $5, $5+1, $6, $8}\' - | tr " " "\\t" | sort -k 1,1 -k2,2n - > ' + imdF
        cmd1=[bedpath +' bedtools' + ' intersect -wa -a ' + imdF + ' -b ' + bli + ' > ' + oli for bli,oli in zip(bedList,outList) ]
        cmd_all=cmd1
        cmd_all[0:0]=[cmd0]
        cmd_all_str=';'.join(cmd_all)
    else:
        cmd1=[bedpath +' bedtools' + ' intersect -wa -a ' + imdF + ' -b ' + bli + ' > ' + oli for bli,oli in zip(bedList,outList) ]
        cmd_all=cmd1
        cmd_all_str=';'.join(cmd_all)
    logobject.info(cmd_all_str)
    with open(os.path.join(outdir,"logs","interval_CpG_prep.out" ),'w') as stdoutF, open(os.path.join(outdir,"logs","interval_CpG_prep.err"),'w') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str       = cmd_all_str,
                                          job_name          = 'int_prep',
                                          logger            = logobject,
                                          drmaa_session     = my_session,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except Exception as err:
            logobject.error("Interval CpG preparation error: %s" % err)
            raise
        else:
            logobject.info('Interval CpG preparation complete')
    return
    

def int_stats_limma(ii,bedList,auxList,sampleSheet,outdir,Rpath,Rlib,pipev,my_session,logobject):
    cmd_all=[os.path.join(Rpath,'Rscript') +' --no-save --no-restore /data/manke/repository/scripts/DNA_methylation/WGBS_pipe/' + pipev + '/WGBSpipe.interval_stats.limma.R ' + outdir + ' ' + li +' '+ aui +' ' + ii + ' ' + sampleSheet + ' ' + Rlib for li,aui in zip(bedList,auxList)]
    cmd_all_str=';'.join(cmd_all)
    logobject.info(cmd_all_str)
    with open(os.path.join(outdir,"logs","interval_stats.out" ),'w') as stdoutF, open(os.path.join(outdir,"logs","interval_stats.err"),'w') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str       = cmd_all_str,
                                          job_name          = 'agg_stats',
                                          logger            = logobject,
                                          drmaa_session     = my_session,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except Exception as err:
            logobject.error("Interval stats error: %s" % err)
            raise
        else:
            logobject.info('Interval stats calculation complete')
    return

