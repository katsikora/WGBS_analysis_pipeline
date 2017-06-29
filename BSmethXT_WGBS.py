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


def methXT_mCT(INfile,QCdir,OUTfile,refG,mCTpath,tabpath,mextDir,mbias_ignore,my_session,logobject):
    read_root=re.sub('.PCRrm.bam','',os.path.basename(INfile))
    pozF=re.sub('.fa','.poz.gz',refG)
    if len(mbias_ignore) < 3:
        mCT_cmd=os.path.join(mCTpath,'methylCtools') + ' bcall --mapQ 10 --snv --trimPE --skipend ' + mbias_ignore + ' --genomic --zero ' + pozF + ' ' + INfile + ' - | ' +  os.path.join(tabpath,'bgzip') + ' > ' + OUTfile 
    elif mbias_ignore=="auto":
        OT=pandas.read_table(os.path.join(QCdir,"logs",read_root)+'.mbias.err',sep=' ',header=None)[4]
        OT=OT.str.split(pat=',').tolist()
        OB=pandas.read_table(os.path.join(QCdir,"logs",read_root)+'.mbias.err',sep=' ',header=None)[6]
        OB=OB.str.split(pat=',').tolist()
        auto_ignore=max([OT[0][0],OT[0][2],OB[0][0],OB[0][2]])
        auto_ignore_corr=min(10,auto_ignore)
        mCT_cmd=os.path.join(mCTpath,'methylCtools') + ' bcall --mapQ 10 --snv --trimPE --skipend ' + auto_ignore_corr + ' --genomic --zero ' + pozF + ' ' + INfile + ' - | ' +  os.path.join(tabpath,'bgzip') + ' > ' + OUTfile 
    logobject.info(mCT_cmd)
    with open(os.path.join(mextDir,"logs","%s.mCT_extract.out" % read_root),'w') as stdoutF, open(os.path.join(mextDir,"logs","%s.mCT_extract.err" % read_root),'w') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str       = mCT_cmd,
                                          job_name          = 'mCT_extract',
                                          logger            = logobject,
                                          drmaa_session     = my_session,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except Exception as err:
            logobject.error("Methylation extraction error: %s" % err)
            raise
        else:
            logobject.info('Methylation extraction complete')
    return

def mCT_prep_poz(refG,mCTpath,tabpath,outdir,pozF,logobject):
    logobject.info('Preparing an index of CG positions')
    cmd_fapos=[os.path.join(mCTpath,'methylCtools') + ' fapos ' + refpath + ' - | ' + os.path.join(tabpath,'bgzip') + ' > ' + pozF]
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

def filt_mCT(INfile,refG,mextDir,my_session,logobject):
    read_root=re.sub('.CG.call.gz','',os.path.basename(INfile))
    pozF=re.sub('.fa','.poz.gz',refG)
    gz_cmd='gzip -dc '+ INfile + ' > '+ re.sub('.gz','',INfile)
    filt_cmd='/package/R-3.3.1/bin/Rscript --no-save --no-restore /data/manke/repository/scripts/DNA_methylation/WGBS_pipe/v0.02/WGBSpipe.mCT.filt.R ' + mextDir + ' ' + re.sub('.gz','',INfile) + ' ' + pozF
    clean_cmd='rm -v ' + re.sub('.gz','',INfile)
    cmd_all=';'.join([gz_cmd,filt_cmd,clean_cmd])
    logobject.info(cmd_all)
    with open(os.path.join(mextDir,"logs","%s.CpG_filt.out" % read_root),'w') as stdoutF, open(os.path.join(mextDir,"logs","%s.CpG_filt.err" % read_root),'w') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str       = cmd_all,
                                          job_name          = 'CpG_filt',
                                          logger            = logobject,
                                          drmaa_session     = my_session,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except Exception as err:
            logobject.error("CpG filtering error: %s" % err)
            raise
        else:
            logobject.info('CpG filtering complete')
    return

def methXT_POM(INfile,QCdir,OUTpfx,refG,POMpath,mextDir,mbias_ignore,nthreads,my_session,logobject):
    read_root=re.sub('.bam','',os.path.basename(INfile)) 
    if len(mbias_ignore) < 3:
        m_ignore=','.join([mbias_ignore]*4)
        POM_cmd=os.path.join(POMpath,'MethylDackel') + ' extract ' + refG + ' ' + INfile + ' -o ' + OUTpfx + ' -q 10 -p 20 --nOT ' + m_ignore + ' --nOB  ' + m_ignore + ' --minDepth 10 --mergeContext --maxVariantFrac 0.25 --minOppositeDepth 5 -@ '+str(nthreads) 
    elif mbias_ignore=="auto":
        with open(os.path.join(QCdir,"logs",read_root)+'.mbias.err', 'r') as f:
            first_line = f.readline()
        mbias_auto=re.sub('Suggested inclusion options: ','',first_line).strip('\n')
        POM_cmd=os.path.join(POMpath,'MethylDackel') + ' extract ' + refG + ' ' + INfile + ' -o ' + OUTpfx + ' -q 10 -p 20 ' + mbias_auto + ' --minDepth 10 --mergeContext --maxVariantFrac 0.25 --minOppositeDepth 5 -@ '+str(nthreads)
    elif len(mbias_ignore) > 4:
        POM_cmd=os.path.join(POMpath,'MethylDackel') + ' extract ' + refG + ' ' + INfile + ' -o ' + OUTpfx + ' -q 10 -p 20 ' + mbias_ignore + ' --minDepth 10 --mergeContext --maxVariantFrac 0.25 --minOppositeDepth 5 -@ '+str(nthreads)
    logobject.info(POM_cmd)
    with open(os.path.join(mextDir,"logs","%s.MetDack_extract.out" % read_root),'w') as stdoutF, open(os.path.join(mextDir,"logs","%s.MetDack_extract.err" % read_root),'w') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str       = POM_cmd,
                                          job_name          = 'MD_extract',
                                          logger            = logobject,
                                          drmaa_session     = my_session,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo --mincpus='+str(nthreads))
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except Exception as err:
            logobject.error("Methylation extraction error: %s" % err)
            raise
        else:
            logobject.info('Methylation extraction complete')
    return
   
def filt_POM(INfile,bedpath,mextDir,my_session,logobject,blackListF=None):
    read_root=re.sub('_CpG.bedGraph','',os.path.basename(INfile))
    Rfilt_cmd='/package/R-3.3.1/bin/Rscript --no-save --no-restore /data/manke/repository/scripts/DNA_methylation/WGBS_pipe/v0.02/WGBSpipe.POM.filt.R ' + mextDir + ' ' + INfile 
    if blackListF is None:
        mv_cmd='mv -v '+ re.sub('_CpG.bedGraph','.CpG.filt.bed',INfile) + ' ' + re.sub('_CpG.bedGraph','.CpG.filt2.bed',INfile)
        cmd_all=';'.join([Rfilt_cmd,mv_cmd])
    else:
        SNP_filt=os.path.join(bedpath,'bedtools') + ' intersect -v -a' + re.sub('_CpG.bedGraph','.CpG.filt.bed',INfile) + ' -b ' + blackListF + ' > ' + re.sub('_CpG.bedGraph','.CpG.filt2.bed',INfile)
        cmd_all=';'.join([Rfilt_cmd,SNP_filt])    logobject.info(cmd_all)
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


def methXT_BisSNP(INfile,QCdir,outSNP,outCG,refG,BisSNPpath,mextDir,mbias_ignore,my_session,nt,logobject):
    read_root=re.sub('.PCRrm.bam','',os.path.basename(INfile))
    if len(mbias_ignore) < 3:
        mXT_cmd='java -Xmx20g -jar ' + BisSNPpath + ' -T BisulfiteGenotyper -C CG,1 -sm BM -nt '+ str(nt) +'  -mmq 10 -mbq 20  -trim5 ' + mbias_ignore +' -trim3 '+ mbias_ignore + ' -I ' + INfile + ' -R ' + refG + ' -vfn1 ' + outCG + ' -vfn2 ' + outSNP
    elif mbias_ignore=="auto":
        OT=pandas.read_table(os.path.join(QCdir,"logs",read_root)+'.mbias.err',sep=' ',header=None)[4]
        OT=OT.str.split(pat=',').tolist()
        OB=pandas.read_table(os.path.join(QCdir,"logs",read_root)+'.mbias.err',sep=' ',header=None)[6]
        OB=OB.str.split(pat=',').tolist()
        auto_ignore=max([OT[0][0],OT[0][2],OB[0][0],OB[0][2]])
        auto_ignore_corr=min(10,auto_ignore)
        mXT_cmd='java -Xmx20g -jar ' + BisSNPpath + ' -T BisulfiteGenotyper -C CG,1 -sm BM -nt '+ str(nt) +'  -mmq 10 -mbq 20  -trim5 ' + auto_ignore_corr +' -trim3 '+ auto_ignore_corr + ' -I ' + INfile + ' -R ' + refG + ' -vfn1 ' + outCG + ' -vfn2 ' + outSNP
    elif len(mbias_ignore) > 4:
        mXT_cmd='java -Xmx20g -jar ' + BisSNPpath + ' -T BisulfiteGenotyper -C CG,1 -sm BM -nt '+ str(nt) +'  -mmq 10 -mbq 20 ' + mbias_ignore +' -I ' + INfile + ' -R ' + refG + ' -vfn1 ' + outCG + ' -vfn2 ' + outSNP
    logobject.info(mXT_cmd)
    with open(os.path.join(mextDir,"logs","%s.BisSNP_extract.out" % read_root),'w') as stdoutF, open(os.path.join(mextDir,"logs","%s.BisSNP_extract.err" % read_root),'w') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str       = mXT_cmd,
                                          job_name          = 'bSNP_extract',
                                          logger            = logobject,
                                          drmaa_session     = my_session,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo --mincpus='+str(nt))
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except Exception as err:
            logobject.error("Methylation extraction error: %s" % err)
            raise
        else:
            logobject.info('Methylation extraction complete')
    return

def vcf_sort(INfile1,INfile2,OUTfile1,OUTfile2,refG,mextDir,my_session,logobject):
    read_root=re.sub('.SNP.vcf','',os.path.basename(INfile1))
    cmd_sortSNP='perl /data/boehm/sikora/tools/sortByRefAndCor.pl --k 1 --c 2 --tmp /data/extended/ ' + INfile1 + ' ' + re.sub('.fa','.fa.fai',refG) + ' > ' + OUTfile1
    cmd_sortCpG='perl /data/boehm/sikora/tools/sortByRefAndCor.pl --k 1 --c 2 --tmp /data/extended/ ' + INfile2 + ' ' + re.sub('.fa','.fa.fai',refG) + ' > ' + OUTfile2
    cmd_all=';'.join([cmd_sortSNP,cmd_sortCpG])
    logobject.info(cmd_all)
    with open(os.path.join(mextDir,"logs","%s.vcf_sort.out" % read_root),'w') as stdoutF, open(os.path.join(mextDir,"logs","%s.vcf_sort.err" % read_root),'w') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str       = cmd_all,
                                          job_name          = 'vcf_sort',
                                          logger            = logobject,
                                          drmaa_session     = my_session,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except Exception as err:
            logobject.error("Vcf sorting error: %s" % err)
            raise
        else:
            logobject.info('Vcf sorting complete')
            zeroFile(INfile1)
            zeroFile(INfile2)
    return

def BisSNP_SNP_filt(INfile,OUTfile,BisSNPpath,refG,mextDir,my_session,logobject):
    read_root=re.sub('.SNP.ro.vcf','',os.path.basename(INfile))
    cmd_SNP='java -Xmx20G -jar ' + BisSNPpath + ' -T VCFpostprocess -R '+ refG + ' -oldVcf '  + INfile + ' -newVcf ' + OUTfile + ' -snpVcf '+ INfile + ' -o ' + re.sub('.vcf','.summary.txt',OUTfile)
    logobject.info(cmd_SNP)
    with open(os.path.join(mextDir,"logs","%s.SNP_filt.out" % read_root),'w') as stdoutF, open(os.path.join(mextDir,"logs","%s.SNP_filt.err" % read_root),'w') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str       = cmd_SNP,
                                          job_name          = 'SNP_filt',
                                          logger            = logobject,
                                          drmaa_session     = my_session,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except Exception as err:
            logobject.error("SNP filtering error: %s" % err)
            raise
        else:
            logobject.info('SNP filtering complete')
    return

def BisSNP_rmSNP(INfile1,INfile2,OUTfile,BisSNPpath,refG,mextDir,my_session,logobject):
    read_root=re.sub('.SNP.filt.vcf','',os.path.basename(INfile1))
    cmd_rmSNP='java -Xmx20G -jar ' + BisSNPpath + ' -T VCFpostprocess -R '+ refG + ' -oldVcf '  + INfile2 + ' -newVcf ' + OUTfile + ' -snpVcf '+ INfile1 + ' -o ' + re.sub('.vcf','.summary.txt',OUTfile)
    logobject.info(cmd_rmSNP)
    with open(os.path.join(mextDir,"logs","%s.SNP_rm.out" % read_root),'w') as stdoutF, open(os.path.join(mextDir,"logs","%s.SNP_rm.err" % read_root),'w') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str       = cmd_rmSNP,
                                          job_name          = 'SNP_rm',
                                          logger            = logobject,
                                          drmaa_session     = my_session,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except Exception as err:
            logobject.error("SNP removal error: %s" % err)
            raise
        else:
            logobject.info('SNP removal complete')
    return

def BisSNP_vcf2bed(INfile,OUTfile,mextDir,my_session,logobject):
    cmd=['grep -v "#" '+ INfile + ' | wc -l ']
    rowNum=int(subprocess.check_output(cmd,shell=True))
    #create empty numpy array of desired size
    dtstr=['string','int','string','string','float','int','int']
    for i in range(7,17):
        dtstr.insert(i,'int')
    colN=['CHROM','POS','STRAND','REF','QUAL','MQ0','DP','DP_FREF','DP_RREF','DP_FALT','DP_RALT','C_Cstrand','T_Cstrand','Illegal_Cstrand','G_Gstrand','A_Gstrand','Illegal_Gstrand']
    dummyA=np.empty(shape=(rowNum,17),dtype='object')
    CpG_df=pd.DataFrame(dummyA,index=range(0,rowNum),columns=colN)
    #read in vcf
    vcf_reader = Reader(open(INfile, 'r'))
    for record, i in zip(vcf_reader, range(0,rowNum)):
        CHROM=pd.Series(record.CHROM,name='CHROM')
        POS=pd.Series(record.POS,name='POS')
        #skip  ALT for CpG calls; skip FILTER
        STRAND=pd.Series(record.INFO['CS'],name='STRAND')
        REF=pd.Series(record.REF,name='REF')
        QUAL=pd.Series(record.QUAL,name='QUAL')
        MQ0=pd.Series(record.INFO['MQ0'],name='MQ0')
        #Context=pd.Series(record.INFO['REF'],name='Context')
        #Sample=pd.Series(record.samples[0].sample,name='Sample') #don't need for single sample vcf
        GT,BQ,BRC6,CM,CP,CU,DP,DP4,GP,GQ,SS=record.samples[0].data
        #skip GT as non-homozygous CGs have been filtered out
        DP=pd.Series(DP,name='DP')
        DP4_dict = cll.OrderedDict(zip(['DP_FREF','DP_RREF','DP_FALT','DP_RALT'], DP4))
        DP4_df=pd.DataFrame.from_dict(DP4_dict,orient='index').transpose()
        BRC6_dict = cll.OrderedDict(zip(['C_Cstrand','T_Cstrand','Illegal_Cstrand','G_Gstrand','A_Gstrand','Illegal_Gstrand'], BRC6))
        BRC6_df=pd.DataFrame.from_dict(BRC6_dict,orient='index').transpose()
        comb_df=pd.concat([CHROM,POS,STRAND,REF,QUAL,MQ0,DP,DP4_df,BRC6_df],axis=1,join='inner')
        CpG_df.iloc[i,]=comb_df.iloc[0,]
    CpG_df.to_csv(OUTfile,sep='\t',na_rep='NA',index=False)
    logobject.info('CpG vcf to txt conversion complete')
    return

def BisSNP_CpGfilt(INfile,mextDir,my_session,logobject):
    read_root=re.sub('.CpG.filt.bed','',os.path.basename(INfile))
    filt_cmd='/package/R-3.3.1/bin/Rscript --no-save --no-restore /data/manke/repository/scripts/DNA_methylation/WGBS_pipe/v0.02/WGBSpipe.BisSNP.filt.R ' + mextDir + ' ' + INfile 
    logobject.info(filt_cmd)
    with open(os.path.join(mextDir,"logs","%s.CpG_filt.out" % read_root),'w') as stdoutF, open(os.path.join(mextDir,"logs","%s.CpG_filt.err" % read_root),'w') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str       = filt_cmd,
                                          job_name          = 'CpG_filt',
                                          logger            = logobject,
                                          drmaa_session     = my_session,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except Exception as err:
            logobject.error("CpG filtering error: %s" % err)
            raise
        else:
            logobject.info('CpG filtering complete')
    return
    
       
    
