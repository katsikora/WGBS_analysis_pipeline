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

def methCT_convert_reads(INfile1,INfile2,mCTpath,readout,my_session,logobject):
    read_root=re.sub('_R1.fastq.gz','',os.path.basename(INfile1))
    R12out=os.path.join(readout,read_root) + '_R12.conv.fastq.gz'
    logobject.info('CT converting the reads')
    cmd_fqconv= os.path.join(mCTpath,'methylCtools') + ' fqconv -1 ' + INfile1 + ' -2 ' + INfile2 + ' - | gzip -c > ' + R12out
    logobject.info(cmd_fqconv)
    with open(os.path.join(readout,"logs","%s.fqconv.out.log" % read_root),'w') as stdoutF, open(os.path.join(readout,"logs","%s.fqconv.err.log" % read_root),'w') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str           = cmd_fqconv,
                                          job_name          = 'fqconv',
                                          logger            = logobject,
                                          drmaa_session     = my_session,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except Exception as err:
            logobject.error("Convert_reads error: %s" % err)
            raise
        else:
            logobject.info('CT converting complete')
    return

def methCT_map_reads(INfile,bwapath,sampath,bamoutDir,refpathBS,nthreads,my_session,logobject):
    read_root=re.sub('_R12.conv.fastq.gz','',os.path.basename(INfile))
    outBam=os.path.join(bamoutDir,read_root)+'.conv.bam'
    #get read group ID string from read name
    with gzip.open(INfile, 'r') as f:
        file_content = f.readline().strip()
    PL=re.sub('@','',file_content).split(":")[0]
    PU=re.sub('@','',file_content).split(":")[2]
    RG='@RG"\t"ID:1"\t"SM:'+read_root+'"\t"LB:'+read_root+'"\t"PL:'+PL+'"\t"PU:'+PU
    mapcmd= os.path.join(bwapath,'bwa') + ' mem -M -p -t ' + str(nthreads) + ' -R ' + RG + ' ' + refpathBS + ' ' + INfile + ' | ' + os.path.join(sampath,'samtools') + ' view -Sb - > ' + outBam
    logobject.info(mapcmd)
    with open(os.path.join(bamoutDir,"logs","%s.readmap.out.log" % read_root),'w') as stdoutF, open(os.path.join(bamoutDir,"logs","%s.readmap.err.log" % read_root),'w') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str           = mapcmd,
                                          job_name          = 'BSmap',
                                          logger            = logobject,
                                          drmaa_session     = my_session,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo --mincpus='+str(nthreads))
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except Exception as err:
            logobject.error("Map_reads error: %s" % err)
            raise
        else:
            logobject.info('Mapping complete')
            zeroFile(INfile)
    return

def methCT_bcov_bam(INfile,mCTpath,sampath,bamoutDir,nthreads,my_session,logobject):
    outBam=re.sub('.conv.bam','.bconv.s.bam',INfile)
    read_root=re.sub('.conv.bam','',os.path.basename(INfile))
    cmd_bconv= os.path.join(mCTpath,'methylCtools') + ' bconv ' + INfile + ' - | ' + os.path.join(sampath,'samtools') + ' sort -T ' + os.path.join('/data/extended',read_root) + ' -m 6G -@ ' + str(nthreads)  + ' -o ' + outBam 
    logobject.info(cmd_bconv)
    with open(os.path.join(bamoutDir,"logs","%s.bconv.out.log" % read_root),'w') as stdoutF, open(os.path.join(bamoutDir,"logs","%s.bconv.err.log" % read_root),'w') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str       = cmd_bconv,
                                          job_name          = 'bconv',
                                          logger            = logobject,
                                          drmaa_session     = my_session,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo --mincpus='+str(nthreads))
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except Exception as err:
            logobject.error("Backconvert_bam error: %s" % err)
            raise
        else:
            logobject.info('Backconverting complete')
            zeroFile(INfile)
    return

def Bism_map_reads(INfile1,INfile2,bismpath,bamoutDir,refpathBS,nthreads,my_session,logobject):
    read_root=re.sub('_R1.fastq.gz','',os.path.basename(INfile1))
    mapcmd= os.path.join(bismpath,'bismark') + ' -p ' + str(nthreads) + ' --dovetail --temp_dir /data/extended --path_to_bowtie /package/bowtie2-2.2.8 --output_dir ' + bamoutDir + ' --basename ' + read_root + ' --genome_folder ' + BSrefpath + ' -1 ' + INfile1 + ' -2 ' + INfile2
    logobject.info(mapcmd)
    with open(os.path.join(bamoutDir,"logs","%s.readmap.out.log" % read_root),'w') as stdoutF, open(os.path.join(bamoutDir,"logs","%s.readmap.err.log" % read_root),'w') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str       = mapcmd,
                                          job_name          = 'BSmap',
                                          logger            = logobject,
                                          drmaa_session     = my_session,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo --mincpus='+str(nthreads))
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except Exception as err:
            logobject.error("Map_reads error: %s" % err)
            raise
        else:
            logobject.info('Mapping complete')
    return

def bMeth_map_reads(INfile1,INfile2,outBam,bmethpath,sampath,bamoutDir,refpathBS,nthreads,my_session,logobject):
    read_root=re.sub('_R1.fastq.gz','',os.path.basename(INfile1))
    # There no point in using more than a few sort threads, it just uses excess memory
    sortThreads = nthreads
    if nthreads > 4:
        sortThreads = 4
    mapcmd=bmethpath + ' bwameth.py --threads ' + str(nthreads)  + ' --reference ' + refpathBS + ' ' + INfile1 + ' ' + INfile2 + ' | ' + os.path.join(sampath,'samtools') + ' sort -T ' + os.path.join('/data/extended',read_root) + ' -m 3G -@ ' + str(sortThreads)  + ' -o ' + outBam
    logobject.info(mapcmd)
    with open(os.path.join(bamoutDir,"logs","%s.readmap.out.log" % read_root),'w') as stdoutF, open(os.path.join(bamoutDir,"logs","%s.readmap.err.log" % read_root),'w') as stderrF:
        try:
            stdout_res, stderr_res = run_job(cmd_str = mapcmd,
                                             job_name = 'BSmap',
                                             logger = logobject,
                                             drmaa_session = my_session,
                                             run_locally = False,
                                             working_directory = os.getcwd(),
                                             job_other_options = '-p bioinfo --mincpus='+str(nthreads))
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except Exception as err:
            logobject.error("Map_reads error: %s" % err)
            raise
        else:
            logobject.info('Mapping complete')
    return

def BS_sort_bam(INfile,sampath,bamoutDir,nthreads,my_session,logobject):
    bamOut=re.sub('.bam','.sorted.bam',INfile)
    read_root=re.sub('.bam','',os.path.basename(INfile))
    sort_cmd=os.path.join(sampath,'samtools') + ' sort -T ' + os.path.join('/data/extended',read_root) + ' -o ' + bamOut + ' ' + INfile
    logobject.info(sort_cmd)
    with open(os.path.join(bamoutDir,"logs","%s.bamsort.out.log" % read_root),'w') as stdoutF, open(os.path.join(bamoutDir,"logs","%s.bamsort.err.log" % read_root),'w') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str       = sort_cmd,
                                          job_name          = 'bamSort',
                                          logger            = logobject,
                                          drmaa_session     = my_session,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo --mincpus='+str(nthreads))
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except Exception as err:
            logobject.error("Bam_sort error: %s" % err)
            raise
        else:
            logobject.info('Bam sorting complete')
            zeroFile(INfile)
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

def BS_rm_dupes(INfile,OUTfile,Picpath,bamoutDir,my_session,logobject):
    read_root=re.sub('.bam','',os.path.basename(INfile))
    rmD_cmd='java -Xmx20g -jar ' + os.path.join(Picpath,'picard.jar') + ' MarkDuplicates I=' + INfile + ' O=' + OUTfile + ' M=' + os.path.join(os.path.dirname(INfile)+re.sub('.bam','.matrix.dup.txt',os.path.basename(INfile))) + ' REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true; ln -f -s ' + re.sub('PCRrm.bam','PCRrm.bai',OUTfile) + ' ' + re.sub('PCRrm.bam','PCRrm.bam.bai',OUTfile)
    logobject.info(rmD_cmd)
    with open(os.path.join(bamoutDir,"logs","%s.rmDupes.out.log" % read_root),'w') as stdoutF, open(os.path.join(bamoutDir,"logs","%s.rmDupes.err.log" % read_root),'w') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str       = rmD_cmd,
                                          job_name          = 'rmPCRdupes',
                                          logger            = logobject,
                                          drmaa_session     = my_session,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo')
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

def BS_add_RGi(INfile,OUTfile,Picpath,bamoutDir,my_session,logobject):
    read_root=re.sub('.bam','',os.path.basename(INfile))
    read_str=subprocess.check_output('samtools view ' + INfile + '| head -n 1',shell=True)
    logobject.debug(read_str)
    PL=read_str.split(":")[0]
    PU=read_str.split(":")[2]
    cmd_addRGI='java -Xmx20g -Djava.io.tmpdir=/data/extended -jar ' +  os.path.join(Picpath,'picard.jar') + ' AddOrReplaceReadGroups I=' + INfile + ' O=' + OUTfile + ' SORT_ORDER=coordinate CREATE_INDEX=true RGPL=' + PL + ' RGSM=' + read_root + ' RGLB=' + read_root + ' RGPU=' + PU + ' VALIDATION_STRINGENCY=SILENT; ln -f -s '+ re.sub('RGi.bam','RGi.bai',OUTfile) + ' ' + re.sub('RGi.bam','RGi.bam.bai',OUTfile)
    logobject.info(cmd_addRGI)
    with open(os.path.join(bamoutDir,"logs","%s.rmDupes.out.log" % read_root),'w') as stdoutF, open(os.path.join(bamoutDir,"logs","%s.addRGi.err.log" % read_root),'w') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str       = cmd_addRGI,
                                          job_name          = 'addRGinfo',
                                          logger            = logobject,
                                          drmaa_session     = my_session,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except Exception as err:
            logobject.error("Add_RG_info error: %s" % err)
            raise
        else:
            logobject.info('Adding read group info complete')
            zeroFile(INfile)
    return
