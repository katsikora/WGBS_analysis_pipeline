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

def single_CpG_limma(ii,sampleInfo,outdir,my_session,logobject):
    Rstat_cmd='/package/R-3.3.1/bin/Rscript --no-save --no-restore /data/manke/repository/scripts/DNA_methylation/WGBS_pipe/v0.02/WGBSpipe.singleCpGstats.limma.R ' + outdir + ' ' + sampleInfo + ' ' + ii
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

def int_stats_limma(ii,bedList,sampleSheet,outdir,my_session,logobject):
    cmd_all=['/package/R-3.3.1/bin/Rscript --no-save --no-restore /data/manke/repository/scripts/DNA_methylation/WGBS_pipe/v0.02/WGBSpipe.interval_stats.limma.R ' + outdir + ' ' + li + ' ' + ii + ' ' + sampleSheet for li in bedList]
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

def single_CpG_DSS(ii,sampleInfo,outdir,my_session,logobject):
    Rstat_cmd='/package/R-3.3.1/bin/Rscript --no-save --no-restore /data/manke/repository/scripts/DNA_methylation/WGBS_pipe/v0.02/WGBSpipe.singleCpGstats.R ' + outdir + ' ' + sampleInfo + ' ' + ii
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


def int_stats_DSS(ii,mpath,bedList,sampleSheet,outdir,my_session,logobject):
    cmd_all=['/package/R-3.3.1/bin/Rscript --no-save --no-restore /data/manke/repository/scripts/DNA_methylation/WGBS_pipe/v0.02/WGBSpipe.interval_stats.DSS.R ' + outdir + ' ' + li + ' ' + ii + ' ' + mpath + ' ' + sampleSheet for li in bedList]
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
