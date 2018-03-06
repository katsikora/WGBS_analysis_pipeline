#####IMPORTS####################################################

import os
import re
import logging
from ruffus.drmaa_wrapper import run_job, error_drmaa_job
import gzip
import subprocess

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

def bMeth_map_reads(INfile1,INfile2,outBam,bmethpath,sampath,bamoutDir,refpathBS,nthreads,my_session,logobject):
    read_root=re.sub('_R1.fastq.gz','',os.path.basename(INfile1))
    with gzip.open(INfile1, 'r') as f:
        file_content = f.readline().strip()
    PL=re.sub('@','',file_content).split(":")[0]
    PU=re.sub('@','',file_content).split(":")[2]
    RG='@RG"\t"ID:1"\t"SM:'+read_root+'"\t"LB:'+read_root+'"\t"PL:'+PL+'"\t"PU:'+PU
    # There no point in using more than a few sort threads, it just uses excess memory
    sortThreads = nthreads
    if nthreads > 4:
        sortThreads = 4
    mapcmd=bmethpath + ' bwameth.py --threads ' + str(nthreads)  + ' --read-group '+ RG + ' --reference ' + refpathBS + ' ' + INfile1 + ' ' + INfile2 + ' | ' + os.path.join(sampath,'samtools') + ' sort -T ' + os.path.join('/data/extended',read_root) + ' -m 3G -@ ' + str(sortThreads)  + ' -o ' + outBam
    logobject.info(mapcmd)
    with open(os.path.join(bamoutDir,"logs","%s.readmap.out.log" % read_root),'w') as stdoutF, open(os.path.join(bamoutDir,"logs","%s.readmap.err.log" % read_root),'w') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str       = mapcmd,
                                          job_name          = 'BSmap',
                                          logger            = logobject,
                                          drmaa_session     = my_session,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo --mem-per-cpu=10000 --nodes=1=1 --mincpus='+str(nthreads))
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except Exception as err:
            logobject.error("Map_reads error: %s" % err)
            raise
        else:
            logobject.info('Mapping complete')
    return

    
def BS_index_bam(INfile,sampath,bamoutDir,my_session,logobject):
    cmd_bamInd=os.path.join(sampath,'samtools') + ' index ' + INfile  + ';sleep 300'
    read_root=re.sub('.bam','',os.path.basename(INfile))
    logobject.info(cmd_bamInd)
    with open(os.path.join(bamoutDir,"logs","%s.bamIndex.out.log" % read_root),'w') as stdoutF, open(os.path.join(bamoutDir,"logs","%s.bamIndex.err.log" % read_root),'w') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str       = cmd_bamInd,
                                          job_name          = 'bamIndex',
                                          logger            = logobject,
                                          drmaa_session     = my_session,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except Exception as err:
            logobject.error("Bam_index error: %s" % err)
            raise
        else:
            logobject.info('Bam indexing complete')
    return

def BS_rm_dupes(INfile,OUTfile,sbpath,bamoutDir,nthreads,my_session,logobject):
    read_root=re.sub('.bam.bai','',os.path.basename(INfile))
    INfile0=os.path.join(bamoutDir,read_root)+'.bam'
    rmD_cmd=sbpath + 'sambamba_v0.6.6 markdup --remove-duplicates -t ' + str(nthreads) + ' --tmpdir /data/extended ' + INfile0 + ' ' + OUTfile + ';sleep 300'
    logobject.info(rmD_cmd)
    with open(os.path.join(bamoutDir,"logs","%s.rmDupes.out.log" % read_root),'w') as stdoutF, open(os.path.join(bamoutDir,"logs","%s.rmDupes.err.log" % read_root),'w') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str       = rmD_cmd,
                                          job_name          = 'rmPCRdupes',
                                          logger            = logobject,
                                          drmaa_session     = my_session,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo --mem-per-cpu=3000 --nodes=1=1 --mincpus='+str(nthreads))
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except Exception as err:
            logobject.error("Remove_PCR_duplicates error: %s" % err)
            raise
        else:
            logobject.info('PCR duplicate removal complete')
            zeroFile(INfile)
    return


