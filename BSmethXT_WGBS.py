#####IMPORTS####################################################

import os
import re
import logging
from ruffus.drmaa_wrapper import run_job, error_drmaa_job
import subprocess
import string
import numpy as np
import pandas
from vcf import Reader
import collections as cll

#####DEFINITIONS#################################################

#Ruffus-compatible way of removing intermediary files
# truncate a file to zero bytes, and preserve its original modification time
def zeroFile(file):
    if os.path.exists(file):
        # save the current time of the file
        timeInfo = os.stat(file)
        try:
            f = open(file,'w')
        except IOError:
            pass
        else:
            f.truncate(0)
            f.close()
            # change the time of the file back to what it was
            os.utime(file,(timeInfo.st_atime, timeInfo.st_mtime))



def mCT_prep_poz(refG,mCTpath,tabpath,outdir,pozF,logobject):
    logobject.info('Preparing an index of CG positions')
    cmd_fapos=os.path.join(mCTpath,'methylCtools') + ' fapos ' + refG + ' - | ' + os.path.join(tabpath,'bgzip') + ' > ' + pozF
    prep_poz='gzip -dc ' + pozF + '| grep \'+\'  | awk \'{print $1, $2, $2+1, $3, $4, $5}\' - | tr " " "\t" > ' + re.sub('gz','P.txt',pozF) + ' ; gzip -dc ' + pozF + ' | grep \'-\'  | awk \'{print $1, $2, $2+1, $3, $4, $5}\' - | tr " " "\t" > ' + re.sub('gz','M.txt',pozF)
    cmd_all=';'.join([cmd_fapos,prep_poz])
    logobject.info(cmd_all)
    try:
        logobject.debug(subprocess.check_call(cmd_all,shell=True))
    except Exception as poze:
        logobject.error("CpG index preparation error: %s" % poze)
        raise IOError
    else:
        logobject.info('CpG index preparation complete')
    return    


def methXT_POM(INfile,QCdir,OUTpfx,refG,POMpath,mextDir,mbias_ignore,nthreads,my_session,logobject):
    read_root=re.sub('.bam','',os.path.basename(INfile)) 
    if len(mbias_ignore) < 3:
        m_ignore=','.join([mbias_ignore]*4)
        POM_cmd=os.path.join(POMpath,'MethylDackel') + ' extract ' + refG + ' ' + INfile + ' -o ' + OUTpfx + ' -q 10 -p 20 --nOT ' + m_ignore + ' --nOB  ' + m_ignore + ' --minDepth 10 --mergeContext --maxVariantFrac 0.25 --minOppositeDepth 5 -@ '+str(nthreads) + ';sleep 300'
    elif mbias_ignore=="auto":
        with open(os.path.join(QCdir,"logs",read_root)+'.mbias.err', 'r') as f:
            first_line = f.readline()
        mbias_auto=re.sub('Suggested inclusion options: ','',first_line).strip('\n')
        POM_cmd=os.path.join(POMpath,'MethylDackel') + ' extract ' + refG + ' ' + INfile + ' -o ' + OUTpfx + ' -q 10 -p 20 ' + mbias_auto + ' --minDepth 10 --mergeContext --maxVariantFrac 0.25 --minOppositeDepth 5 -@ '+str(nthreads)  + ';sleep 300'
    elif len(mbias_ignore) > 4:
        POM_cmd=os.path.join(POMpath,'MethylDackel') + ' extract ' + refG + ' ' + INfile + ' -o ' + OUTpfx + ' -q 10 -p 20 ' + mbias_ignore + ' --minDepth 10 --mergeContext --maxVariantFrac 0.25 --minOppositeDepth 5 -@ '+str(nthreads)  + ';sleep 300'
    logobject.info(POM_cmd)
    with open(os.path.join(mextDir,"logs","%s.MetDack_extract.out" % read_root),'w') as stdoutF, open(os.path.join(mextDir,"logs","%s.MetDack_extract.err" % read_root),'w') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str       = POM_cmd,
                                          job_name          = 'MD_extract',
                                          logger            = logobject,
                                          drmaa_session     = my_session,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo --nodes=1=1 --mincpus='+str(nthreads))
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except Exception as err:
            logobject.error("Methylation extraction error: %s" % err)
            raise
        else:
            logobject.info('Methylation extraction complete')
    return
   
def filt_POM(INfile,bedpath,mextDir,Rpath,Rlib,pipev,my_session,logobject,blackListF=None):
    read_root=re.sub('_CpG.bedGraph','',os.path.basename(INfile))
    Rfilt_cmd=os.path.join(Rpath,'Rscript') +' --no-save --no-restore /data/manke/repository/scripts/DNA_methylation/WGBS_pipe/' + pipev + '/WGBSpipe.POM.filt.R ' + mextDir + ' ' + INfile + ' ' + Rlib 
    if blackListF is None:
        mv_cmd='mv -v '+ re.sub('_CpG.bedGraph','.CpG.filt.bed',INfile) + ' ' + re.sub('_CpG.bedGraph','.CpG.filt2.bed',INfile) + ';sleep 300'
        cmd_all=';'.join([Rfilt_cmd,mv_cmd])
    else:
        SNP_filt=bedpath +' bedtools' + ' intersect -v -a' + re.sub('_CpG.bedGraph','.CpG.filt.bed',INfile) + ' -b ' + blackListF + ' > ' + re.sub('_CpG.bedGraph','.CpG.filt2.bed',INfile) + ';sleep 300'
        cmd_all=';'.join([Rfilt_cmd,SNP_filt])
    logobject.info(cmd_all)
    with open(os.path.join(mextDir,"logs","%s.MetDack_filt.out" % read_root),'w') as stdoutF, open(os.path.join(mextDir,"logs","%s.MetDack_filt.err" % read_root),'w') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str       = cmd_all,
                                          job_name          = 'MD_filt',
                                          logger            = logobject,
                                          drmaa_session     = my_session,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except Exception as err:
            logobject.error("Methylation filtering error: %s" % err)
            raise
        else:
            logobject.info('Methylation filtering complete')
    return


   
       
    
