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

def DMR_metilene(ii,sampleInfo,outfile,nthreads,metipath,my_session,logobject):
    outdir=os.path.dirname(outfile)
    met_cmd=os.path.join(metipath,'metilene')+ ' -a ' + list(set(pd.read_table(sampleInfo)['Group']))[0] + ' -b ' + list(set(pd.read_table(sampleInfo)['Group']))[1] + ' -t ' + str(nthreads) + ' ' + ii + ' | sort -k 1,1 -k2,2n > ' + outfile + ';sleep 300'
    logobject.info(met_cmd)
    with open(os.path.join(outdir,"logs","DMR.metilene.out" ),'w') as stdoutF, open(os.path.join(outdir,"logs","DMR.metilene.err"),'w') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str       = met_cmd,
                                          job_name          = 'metilene',
                                          logger            = logobject,
                                          drmaa_session     = my_session,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo --mincpus='+str(nthreads))
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except Exception as err:
            logobject.error("Metilene error: %s" % err)
            raise
        else:
            logobject.info('Metilene DMR calling complete')
    return


def clean_up_metilene(metilene_out,CpG_stats_out,sampleInfo,outdir,my_session,logobject):
    cmd='/package/R-3.3.1/bin/Rscript --no-save --no-restore /data/manke/repository/scripts/DNA_methylation/WGBS_pipe/v1.0.0/WGBSpipe.metilene_stats.limma.R ' + outdir + ' ' + metilene_out + ' ' + CpG_stats_out +' ' + sampleInfo
    logobject.info(cmd)
    with open(os.path.join(outdir,"logs","metilene.cleanup.out" ),'w') as stdoutF, open(os.path.join(outdir,"logs","metilene.cleanup.err"),'w') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str       = cmd,
                                          job_name          = 'cleanup',
                                          logger            = logobject,
                                          drmaa_session     = my_session,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except Exception as err:
            logobject.error("Metilene cleanup error: %s" % err)
            raise
        else:
            logobject.info('Metilene cleanup complete')
    return
    
    
