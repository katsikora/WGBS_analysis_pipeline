#####IMPORTS####################################################

import os
import re
import logging
from ruffus.drmaa_wrapper import run_job, error_drmaa_job
import subprocess
import string
import shutil

#####DEFINITIONS#################################################

def mCT_get_ranCG(INfile,OUTfile,refG,pozDir,mCTpath,logobject):
    os.chdir(pozDir)
    if os.path.exists(OUTfile):
        return
    pozF=os.path.join(pozDir,re.sub('.fa*','.poz.gz',os.path.basename(refG)))
    cmd_from0='source activate NGSpy2.7; ' + os.path.join(mCTpath,'methylCtools') + ' fapos ' + refG + ' ' + re.sub('.gz','',pozF) + ';cat '+ re.sub('.gz','',pozF) +' | grep "+" - | shuf | head -n 1000000 | awk \'{print $1, $5, $5+1, $6, $8}\' - | tr " " "\\t" | sort -k 1,1 -k2,2n - > ' + OUTfile +' ;source deactivate'
    cmd_from_poz='zcat ' +  pozF + ' | grep "+" - | shuf | head -n 1000000 | awk \'{print $1, $5, $5+1, $6, $8}\' - | tr " " "\\t" | sort -k 1,1 -k2,2n - > ' + OUTfile 
    if os.path.exists(pozF):
        logobject.info(cmd_from_poz)
        logobject.debug(subprocess.check_output(cmd_from_poz,shell=True))
    else:
        logobject.info(cmd_from0)
        logobject.debug(subprocess.check_output(cmd_from0,shell=True))
    logobject.info('Random 1mln CpG indexing complete')
    return

def BS_doc_XT(INfile,OUTfileList,bedFiles,auxdir,refG,GATKpath,metDir,my_session,logobject):
    read_root=re.sub('.bam','',os.path.basename(INfile))
    OUTlist2=OUTfileList[2:]
    allExist = True
    for f in OUTfileList:
        if not os.path.exists(f):
            allExist = False
    if allExist:
        return
    WG_mean_cmd='java -Xmx30g -Djava.io.tmpdir=/data/extended -jar '+ os.path.join(GATKpath,'GenomeAnalysisTK.jar')+' -R '+ refG + ' -T DepthOfCoverage -o ' + str(OUTfileList[0]) + ' -I ' + INfile + ' -ct 0 -ct 1 -ct 2 -ct 5 -ct 10 -ct 15 -ct 20 -ct 30 -ct 50  -omitBaseOutput -mmq 10 --partitionType sample '
    CG_cmd='java -Xmx50g -Djava.io.tmpdir=/data/extended -jar '+ os.path.join(GATKpath,'GenomeAnalysisTK.jar')+' -R '+ refG + ' -T DepthOfCoverage -o ' + str(OUTfileList[1]) + ' -I ' + INfile + ' -ct 0 -ct 1 -ct 2 -ct 5 -ct 10 -ct 15 -ct 20 -ct 30 -ct 50  -omitBaseOutput -omitIntervals -mmq 10 --partitionType sample -L ' + os.path.join(auxdir,re.sub('.fa','.poz.ran1M.sorted.bed',os.path.basename(refG)))
    cmd_Xi=['java -Xmx30g -Djava.io.tmpdir=/data/extended -jar '+ os.path.join(GATKpath,'GenomeAnalysisTK.jar')+' -R '+ refG + ' -T DepthOfCoverage -o ' + oi + ' -I ' + INfile + ' -ct 0 -ct 1 -ct 2 -ct 5 -ct 10 -ct 15 -ct 20 -ct 30 -ct 50  -omitBaseOutput -omitIntervals -mmq 10 --partitionType sample -L ' + bi for oi,bi in zip(OUTlist2,bedFiles)]
    cmd_all=cmd_Xi
    cmd_all[0:0]=[CG_cmd]
    cmd_all[0:0]=[WG_mean_cmd]
    cmd_all_str=';'.join(cmd_all)
    logobject.info(cmd_all_str)
    with open(os.path.join(metDir,"logs","%s.depth_cov.out.log" % read_root),'w') as stdoutF, open(os.path.join(metDir,"logs","%s.depth_cov.err.log" % read_root),'w') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str       = cmd_all_str,
                                          job_name          = 'depth_cov',
                                          logger            = logobject,
                                          drmaa_session     = my_session,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo --mem-per-cpu=30000')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except Exception as err:
            logobject.error("Depth of coverage error: %s" % err)
            raise
        else:
            logobject.info('Depth of coverage calculation complete')
    return

def BS_doc(INfile,OUTfileList,refG,auxdir,GATKpath,metDir,my_session,logobject):
    allExist = True
    for f in OUTfileList:
        if not os.path.exists(f):
            allExist = False
    if allExist:
        return
    read_root=re.sub('.bam','',os.path.basename(INfile))
    WG_mean_cmd='java -Xmx30g -Djava.io.tmpdir=/data/extended -jar '+ os.path.join(GATKpath,'GenomeAnalysisTK.jar')+' -R '+ refG + ' -T DepthOfCoverage -o ' + str(OUTfileList[0]) + ' -I ' + INfile + ' -ct 0 -ct 1 -ct 2 -ct 5 -ct 10 -ct 15 -ct 20 -ct 30 -ct 50  -omitBaseOutput -mmq 10 --partitionType sample '
    CG_cmd='java -Xmx30g -Djava.io.tmpdir=/data/extended -jar '+ os.path.join(GATKpath,'GenomeAnalysisTK.jar')+' -R '+ refG + ' -T DepthOfCoverage -o ' + str(OUTfileList[1]) + ' -I ' + INfile + ' -ct 0 -ct 1 -ct 2 -ct 5 -ct 10 -ct 15 -ct 20 -ct 30 -ct 50  -omitBaseOutput -omitIntervals -mmq 10 --partitionType sample -L ' + os.path.join(auxdir,re.sub('.fa','.poz.ran1M.sorted.bed',os.path.basename(refG)))
    cmd_all=[WG_mean_cmd,CG_cmd]
    cmd_all_str=';'.join(cmd_all)
    logobject.info(cmd_all_str)
    with open(os.path.join(metDir,"logs","%s.depth_cov.out.log" % read_root),'w') as stdoutF, open(os.path.join(metDir,"logs","%s.depth_cov.err.log" % read_root),'w') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str       = cmd_all_str,
                                          job_name          = 'depth_cov',
                                          logger            = logobject,
                                          drmaa_session     = my_session,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo --mem-per-cpu=30000')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except Exception as err:
            logobject.error("Depth of coverage error: %s" % err)
            raise
        else:
            logobject.info('Depth of coverage calculation complete')
    return

def BS_conv_rate(ii1sub,oo,metDir,my_session,logobject):
    read_root=os.path.basename(ii1sub)
    CR_cmd='/data/manke/repository/scripts/DNA_methylation/DEEP_scripts/conversionRate_KS.sh '+ ii1sub + ' ' + oo
    logobject.info(CR_cmd)
    with open(os.path.join(metDir,"logs","%s.conv_rate.out.log" % read_root),'w') as stdoutF, open(os.path.join(metDir,"logs","%s.conv_rate.err.log" % read_root),'w') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str       = CR_cmd,
                                          job_name          = 'conv_rate',
                                          logger            = logobject,
                                          drmaa_session     = my_session,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except Exception as err:
            logobject.error("Conversion rate error: %s" % err)
            raise
        else:
            logobject.info('Conversion rate calculation complete')
    return

def BS_Mbias(INfile,OUTfile,POMpath,refG,metDir,nthreads,my_session,logobject):
    read_root=re.sub('.bam','',os.path.basename(INfile))
    Mb_cmd=os.path.join(POMpath,'MethylDackel') + ' mbias --txt ' + refG + ' ' + INfile +' ' + OUTfile +' -@ '+str(nthreads) +' > ' + OUTfile + '.txt'
    logobject.info(Mb_cmd)
    with open(os.path.join(metDir,"logs","%s.mbias.out" % read_root),'w') as stdoutF, open(os.path.join(metDir,"logs","%s.mbias.err" % read_root),'w') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str       = Mb_cmd,
                                          job_name          = 'mbias',
                                          logger            = logobject,
                                          drmaa_session     = my_session,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo --mem-per-cpu=3000 --nodes=1=1 --mincpus='+str(nthreads))
            
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))


        # relay all the stdout, stderr, drmaa output to diagnose failures
        except Exception as err:
            logobject.error("Methylation bias error: %s" % err)
            raise
        else:
            logobject.info('Methylation bias calculation complete')
    return
    
def BS_flagstat(INfile,OUTfile,sampath,outdir,my_session,logobject):
    if os.path.exists(OUTfile):
        return
    read_root=re.sub('.bam','',os.path.basename(INfile))
    cmd=os.path.join(sampath,'samtools') + ' flagstat ' + INfile +' > ' + OUTfile 
    logobject.info(cmd)
    with open(os.path.join(outdir,"logs","%s.flagstat.out" % read_root),'w') as stdoutF, open(os.path.join(outdir,"logs","%s.flagstat.err" % read_root),'w') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str       = cmd,
                                          job_name          = 'fstat',
                                          logger            = logobject,
                                          drmaa_session     = my_session,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo ' )
            
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))


        # relay all the stdout, stderr, drmaa output to diagnose failures
        except Exception as err:
            logobject.error("Flagstat error: %s" % err)
            raise
        else:
            logobject.info('Flagstat calculation complete')
    return
 
def BS_QC_rep(Rpath,QCdir,logobject):
    shutil.copyfile("{}/WGBS_QC_report_template.Rmd".format(os.path.dirname(os.path.realpath(__file__))), os.path.join(QCdir, "WGBS_QC_report_template.Rmd"))
    cmd = os.path.join(Rpath, 'Rscript -e "rmarkdown::render(\''+os.path.join(QCdir, "WGBS_QC_report_template.Rmd")+'\', params=list(QCdir=\'"' + QCdir +'"\' ), output_file =\'"'+ os.path.join(QCdir,'QC_report.pdf"\'')+')"')
    logobject.info(cmd)
    logobject.debug(subprocess.check_output(cmd,shell=True))
    logobject.info('Generating QC report complete')
    
