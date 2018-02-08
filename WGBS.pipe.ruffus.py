import os
import sys
import time
import argparse

#to-do: reduce the imports, take too long
import zipfile
import re
import subprocess
import pandas
import numpy
import commands
import math
import time
from collections import OrderedDict
import logging
from operator import is_not
from functools import partial
import statistics
import pysam

sys.path.insert(0, "/data/manke/repository/scripts/DNA_methylation/ruffus")
from ruffus import *
from ruffus.proxy_logger import *
from ruffus.combinatorics import *
from ruffus.drmaa_wrapper import run_job, run_job_using_drmaa, error_drmaa_job

from PIL import Image
import string

parser = argparse.ArgumentParser(prog="methPipe", version="0.0.2", description="runs complete CpG methylation pipeline", formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-ri", "--readIn", dest="readdir", action="store", default=False, help="input read folder")
parser.add_argument("-bi", "--bamIn", dest="bamdir", action="store", default=False, help="input bam folder")
parser.add_argument("-mi", "--methTabIn", dest="methtabdir", action="store", default=False, help="input methylation table folder")
parser.add_argument("-qi", "--fqcIn", dest="fqcdir", action="store", default=False, help="folder with fastqc.zip results for raw reads")
parser.add_argument("-g", "--ref", dest="refpath", action="store", default=False, help="path to indexed reference genome")
parser.add_argument("-cg", "--cref", dest="convrefpath", action="store", default=False, help="path to converted and indexed reference genome")
parser.add_argument("-il", "--intList", dest="intList", action="append", default=None, help="target interval file(s)")
parser.add_argument("-bl", "--blackList", dest="blackList", action="store", default=None, help="SNP black list")
parser.add_argument("-si", "--sampleInfo", dest="sampleInfo", action="store", default=None, help="sample sheet")
parser.add_argument("-wd", "--wdir", dest="wdir", action="store", default=False, help="output folder")
parser.add_argument("-vb", "--verbose", dest="verbose", action="store_true", help="more detailed log")
parser.add_argument("-bs", "--batchSize", dest="bsize", action="store",metavar="INT",type=int,default=10, help="number of samples to process in parallel")
parser.add_argument("-nt", "--numThr", dest="nthreads", action="store",metavar="INT",type=int,default=8, help="number of threads to use per sample")
parser.add_argument("-tR", "--trimReads", dest="trimReads", action="store_true", help="adapter-trim and hardclip the reads")
parser.add_argument("--nextera", action="store_true", help="Trim Nextera adapters, rather than the default TruSeq adapters. Requires --trimReads")
parser.add_argument("--trimThreshold", default=10, type=int, help="If --trimReads is specified, this sets the phred score threshold.")
parser.add_argument("--trimOtherArgs", default="", help="Other arguments you would like passed to cutadapt.")
parser.add_argument("-cr", "--convRef", dest="convRef", action="store_true", help="BS-convert reference genome")
parser.add_argument('--aligner', choices=['methylCtools', 'Bismark', 'bwaMeth'],dest="aligner",action="store",default="bwaMeth",help="mapping software to use")
parser.add_argument('--extractor', choices=['methylCtools', 'MethylDackel', 'BisSNP'],dest="methEXT",action="store",default="MethylDackel",help="methylation extraction software to use")
parser.add_argument('--stats', choices=['limma', 'DSS'],dest="stats",action="store",default="limma",help="stats package to use")
parser.add_argument('--mbias', dest="mbias_ignore",action="store",default="auto",help="number of nucleotides with mbias to ignore during methylation extraction")
parser.add_argument('--DMRpg', choices=['methylSeekR', 'metilene'],dest="DMRpg",action="store",default="None",help="stats package to use")
parser.add_argument("--touchOnly", dest="touchOnly", action="store_true", help="only touch files")
parser.add_argument("--target_tasks", dest="target_tasks", action="store",default=[], help="target tasks")
parser.add_argument("--forcedtorun_tasks", dest="forcedtorun_tasks", action="store",default=[], help="forced to run tasks")

args = parser.parse_args()

#setup central working directory
wdir=args.wdir
if not os.path.exists(wdir):
    os.makedirs(wdir)
    
#setup logging
logger = logging.getLogger(__name__)
fhandler = logging.FileHandler(filename=os.path.join(wdir,'pipeline.log'), mode='a')
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fhandler.setFormatter(formatter)
logger.addHandler(fhandler)
logger.setLevel(logging.DEBUG)

logger.debug(subprocess.check_output('echo $DRMAA_LIBRARY_PATH',shell=True))

import drmaa

#initiate 1 drmaa session for the whole pipeline
mySession=drmaa.Session()
mySession.initialize()


logger.info('Working directory is : ' + wdir)
logger.info(args)

#identify pipeline input files
if args.readdir:
    readdir=args.readdir
    os.chdir(readdir)

    libset2=[ os.path.join(readdir,x) for x in os.listdir(readdir) ]
    libset2_R1=filter(lambda x:'_R1.fastq.gz' in x, libset2)
    libset2_R1.sort()
    libset2_R2=filter(lambda x:'_R2.fastq.gz' in x, libset2)
    libset2_R2.sort()
    read_root=[ re.sub('_R1.fastq.gz','',R1f) for R1f in libset2_R1 ]
    INfiles=list(zip(libset2_R1,libset2_R2))
    
elif args.bamdir:
    bamdir=args.bamdir
    os.chdir(bamdir)
    INfiles0=[ os.path.join(bamdir,x) for x in os.listdir(bamdir) ]
    infiles=[re.match( '.+\.bam$', x) for x in INfiles0] ####use .bam$ only later on
    ifl=filter(partial(is_not,None),infiles)
    INfiles=[ x.group() for x in ifl ]
    INfiles.sort()
elif args.methtabdir:
    methtabdir=args.methtabdir
    os.chdir(methtabdir)
    INfiles0=[os.path.join(methtabdir,x) for x in os.listdir(methtabdir)]
    if args.methEXT=='methylCtools':
        infiles=[re.match( '.+\.CG\.call\.gz', x) for x in INfiles0]  
    elif args.methEXT=='BisSNP': 
        infiles=[re.match( '.+\.CpG\.filt\.bed', x) for x in INfiles0]
    elif args.methEXT=='MethylDackel': 
        infiles=[re.match( '.+_CpG\.bedGraph', x) for x in INfiles0]
    ifl=filter(partial(is_not,None),infiles)
    INfiles=[ x.group() for x in ifl ]
    INfiles.sort()
    
logger.debug(INfiles)

######################## paths to unconverted reference genome #######################################
if os.path.isfile(args.refpath):
    refG=args.refpath
    refGpath=os.path.dirname(refG)
    if args.aligner=="methylCtools" :
        logger.warning('Make sure the reference chromosome names contain no spaces!')
elif os.path.isfile(os.path.join('/data/repository/organisms',(args.refpath+'_ensembl'),'genome_fasta/genome.fa')):
    refG=os.path.join('/data/repository/organisms',(args.refpath+'_ensembl'),'genome_fasta/genome.fa')
    refGpath=os.path.dirname(refG)
    if args.aligner=="methylCtools" :
        logger.warning('Make sure the reference chromosome names contain no spaces!')
else:
    logger.info('Reference genome not recognized. Please check spelling and/or path.')
        
######################## paths to converted reference genome #######################################
ref_conv_dict={'bwaMeth':'BWAmethIndex','Bismark':'BismarkIndex/Bisulfite_Genome','methylCtools':'NA'}    

if os.path.isfile(str(args.convrefpath)):
    crefG=args.convrefpath
    crefGpath=os.path.dirname(crefG)
elif os.path.isfile(os.path.join('/data/repository/organisms',(str(args.refpath)+'_ensembl'),ref_conv_dict[args.aligner],'genome.fa')):
    
    crefG=os.path.join('/data/repository/organisms',(args.refpath+'_ensembl'),ref_conv_dict[args.aligner],'genome.fa')
    crefGpath=os.path.dirname(crefG)
elif not args.convrefpath:
    print('Converted reference genome not specified, will be generated by the pipeline.')
    args.convRef=True
else:
    logger.error('Converted reference genome not recognized. Please check spelling and/or path.')    
    
logger.debug('Reference genome is '+refG)
logger.debug('Converted reference genome is '+crefG)
logger.debug('Convert genome? ' + str(args.convRef))

##################PATHS TO EXECUTABLES###############################################################
FQCpath='/package/FastQC-0.11.3'
cutpath='/package/anaconda3/envs/cutadapt-1.15/bin'
mCTpath='/data/manke/repository/scripts/DNA_methylation/methylCtools'
tabpath='/package/tabix-1.2.1'
bwapath='/package/anaconda3/envs/bwa-0.7.17/bin'
bismpath='/data/manke/repository/scripts/DNA_methylation/Bismark-master'
BTpath='/package/bowtie2-2.2.8/'
#bmethpath='module load bwameth && /package/anaconda3/envs/bwameth-0.2.0/bin'
bmethpath='module load bwameth; '
sampath='/package/samtools-1.3/bin'
Picpath='/package/picard-tools-1.136'
GATKpath='/package/GenomeAnalysisTK-3.5'
POMpath='/package/MethylDackel-0.3.0/bin'
bedpath='/package/bedtools2-2.25.0/bin'
BisSNPpath='/data/manke/repository/scripts/DNA_methylation/BisSNP-0.82.2.jar'
metipath='/data/manke/repository/scripts/DNA_methylation/metilene_v0.2-6'
Rpath='/package/R-3.3.1/bin'

######################################################################################################

#add folder with custom modules to path            
sys.path.append(os.path.dirname(os.path.realpath(__file__)))

#############################PIPELINE#################################################################################
#TRIM READS############################################################################################################
#if args.trimReads:
#    fqcout=os.path.join(wdir,'fastqc_cut')
#    cutout=os.path.join(wdir,'reads_cut')
#    #setup shared logging 
#    log_args={}
#    log_args["file_name"] = os.path.join(fqcout,"LOG")
#    log_args["formatter"] = "%(asctime)s - %(name)s - %(levelname)6s - %(message)s"
#    log_args["delay"]     = True
#    log_args["level"]     = logging.DEBUG
#    (shared_logger,logging_mutex) = make_shared_logger_and_proxy (setup_std_shared_logger,
#                                               "shared_logger",
#                                               log_args)  
#    if args.fqcdir:    
#        fqcdir=args.fqcdir
#        fqcL=[ os.path.join(fqcdir,z ) for z in os.listdir(fqcdir) ]
#        zipR1=filter(lambda x:'_R1_fastqc.zip' in x, fqcL)
#        zipR1.sort()
#        zipR2=filter(lambda x:'_R2_fastqc.zip' in x, fqcL)
#        zipR2.sort()
#        zipL=list(zip(zipR1,zipR2))
#        zipL2=[ list(z) for z in zipL ]
#        os.chdir(fqcdir)
#        #from BSprecut import calc_cutThd
#        @mkdir(fqcout)
#        @transform(zipL2,suffix('_R1_fastqc.zip'),'.R12.ct.txt',output_dir=fqcout)
#        def get_cutThd(input_files,output_file):
#            ii1 = input_files[0]
#            ii2 = input_files[1]
#            from BSprecut import calc_cutThd
#            with open(output_file, "w") as oo:
#                cutThdRes=calc_cutThd([ii1,ii2],args.fqcdir,fqcout,shared_logger,logging_mutex)
#                oo.write('\n'.join('%s\t%s\n' % x for x in cutThdRes))
#                    
#        
#    os.chdir(readdir)
    
if args.trimReads:
    @mkdir(cutout,os.path.join(cutout,'logs'))
    @transform(input=INfiles,filter=suffix('_R1.fastq.gz'),output=['_R1.fastq.gz','_R2.fastq.gz'],output_dir=cutout) 
    def trim_reads(input_files,output_files):
        ii1 = input_files[0]
        ii2 = input_files[1]
        oo1 = output_files[0]
        oo2 = output_files[1]
        from BSprecut import cut_reads
        cut_reads(ii1,ii2,oo1,oo2,cutpath,mySession,cutout,logger, args)
    @mkdir(fqcout,os.path.join(fqcout,'logs'))
    @transform(input=trim_reads,filter=suffix('_R1.fastq.gz'),output=['_R1.zip','_R1.html','_R2.zip','_R2.html'],output_dir=fqcout)
    def postTrim_fastqc(input_files,output_files):
        ii1 = input_files[0]
        ii2 = input_files[1]
        from BSprecut import post_trim_fqc
        post_trim_fqc(ii1,ii2,fqcout,FQCpath,mySession,logger)   
        
#CONVERT REFERENCE GENOME##########################################################################################
if args.convRef:
    crefGpath=os.path.join(wdir,'conv_ref_'+args.aligner)
    if not os.path.exists(crefGpath):
        os.makedirs(crefGpath)
    
    if ( args.aligner=='methylCtools' and not args.bamdir and not args.methtabdir ):
        @active_if(args.convRef)
        @transform(refG,suffix('.fa'),'.conv.fa',output_dir=crefGpath)
        def conv_ref(input_file,output_files):
            os.chdir(refGpath)
            ii = input_file
            from BSrefprep import methCT_refprep
            methCT_refprep(mCTpath,tabpath,bwapath,ii,crefGpath,logger)
            os.chdir(wdir)
    elif ( args.aligner=='Bismark' and not args.bamdir and not args.methtabdir ): 
        @active_if(args.convRef)
        @transform(refG,suffix('.fa'),'_mfa.CT_conversion.fa',output_dir=os.path.join(crefGpath,'Bisulfite_Genome','CT_conversion'))
        def conv_ref(input_file,output_files):
            os.chdir(refGpath)
            ii = input_file
            from BSrefprep import Bism_refprep
            Bism_refprep(bismpath,BTpath,ii,crefGpath,logger)
            os.chdir(wdir)
    elif ( args.aligner=='bwaMeth' and not args.bamdir and not args.methtabdir ):
        @active_if(args.convRef)
        @transform(refG,suffix('.fa'),'.fa.bwameth.c2t',output_dir=crefGpath)
        def conv_ref(input_file,output_files):
            os.chdir(refGpath)
            ii = input_file
            from BSrefprep import bmeth_refprep
            bmeth_refprep(bmethpath,ii,crefGpath,logger)
            os.chdir(wdir)
            
            
#CONVERT AND MAP READS############################################################################################
bamout=os.path.join(wdir,'bams')
bamoutO=bamout+'_'+args.aligner
###################################################################################################################

if ( args.aligner=='methylCtools' and not args.bamdir and not args.methtabdir ):
    readout=os.path.join(wdir,'conv_reads')
    if args.trimReads:
        @mkdir(readout,os.path.join(readout,'logs'))
        @transform(trim_reads,suffix('_R1.fastq.gz'),'_R12.conv.fastq.gz',output_dir=readout)
        def convert_reads(input_files,output_file):
            ii1 = input_files[0]
            ii2 = input_files[1]
            oo = output_file
            from BSmapWGBS import methCT_convert_reads
            methCT_convert_reads(ii1,ii2,mCTpath,readout,mySession,logger)
    else:
        @mkdir(readout,os.path.join(readout,'logs'))
        @transform(INfiles,suffix('_R1.fastq.gz'),'_R12.conv.fastq.gz',output_dir=readout)
        def convert_reads(input_files,output_file):
            ii1 = input_files[0]
            ii2 = input_files[1]
            oo = output_file
            from BSmapWGBS import methCT_convert_reads
            methCT_convert_reads(ii1,ii2,mCTpath,readout,mySession,logger)
    if args.convRef:  
        @follows("conv_ref")
        @mkdir(bamoutO,os.path.join(bamoutO,'logs'))
        @transform(convert_reads,suffix('_R12.conv.fastq.gz'),".conv.bam",output_dir=bamoutO)
        def map_reads(input_file,output_file): 
            ii = input_file
            oo = output_file
            crefG=os.path.join(crefGpath,re.sub('.fa','.conv.fa',os.path.basename(refG)))
            from BSmapWGBS import methCT_map_reads
            methCT_map_reads(ii,bwapath,sampath,bamoutO,crefG,args.nthreads,mySession,logger)
                               
    else:
        @mkdir(bamoutO,os.path.join(bamoutO,'logs'))
        @transform(convert_reads,suffix('_R12.conv.fastq.gz'),'.conv.bam',output_dir=bamoutO)
        def map_reads(input_file,output_file):
            ii = input_file
            oo = output_file
            crefG2=os.path.join(crefGpath,re.sub('.fa','.fa.bwameth.c2t',os.path.basename(crefG)))
            from BSmapWGBS import methCT_map_reads
            methCT_map_reads(ii,bwapath,sampath,bamoutO,crefG2,args.nthreads,mySession,logger)
                               
    @transform(map_reads,suffix('.conv.bam'),'.bconv.s.bam',output_dir=bamoutO) 
    def backconvert_bam(input_file,output_file):
        ii = input_file
        oo = output_file
        from BSmapWGBS import methCT_bcov_bam
        methCT_bcov_bam(ii,mCTpath,sampath,bamoutO,args.nthreads,mySession,logger)
                               
elif ( args.aligner=='Bismark' and not args.bamdir and not args.methtabdir ):
    if args.convRef:
        if args.trimReads:
            @follows("conv_ref")
            @mkdir(bamoutO,os.path.join(bamoutO,'logs'))
            @transform(trim_reads,suffix('_R1.fastq.gz'),"_pe.bam",output_dir=bamoutO)
            def map_reads(input_files,output_file): 
                ii1 = input_files[0]
                ii2 = input_files[1]
                oo = output_file
                from BSmapWGBS import Bism_map_reads
                Bism_map_reads(ii1,ii2,bismpath,bamoutO,crefGpath,args.nthreads,mySession,logger)
        else:   
            @follows("conv_ref")
            @mkdir(bamoutO,os.path.join(bamoutO,'logs'))
            @transform(INfiles,suffix('_R1.fastq.gz'),"_pe.bam",output_dir=bamoutO)
            def map_reads(input_files,output_file):
                ii1 = input_files[0]
                ii2 = input_files[1]
                oo = output_file
                from BSmapWGBS import Bism_map_reads
                Bism_map_reads(ii1,ii2,bismpath,bamoutO,crefGpath,args.nthreads,mySession,logger)
    else:
        if args.trimReads:
            @mkdir(bamoutO,os.path.join(bamoutO,'logs'))
            @transform(trim_reads,suffix('_R1.fastq.gz'),'_pe.bam',output_dir=bamoutO)
            def map_reads(input_files,output_file):
                ii1 = input_files[0]
                ii2 = input_files[1]
                oo = output_file
                from BSmapWGBS import Bism_map_reads
                Bism_map_reads(ii1,ii2,bismpath,bamoutO,crefGpath,args.nthreads,mySession,logger)
        else: 
            @mkdir(bamoutO,os.path.join(bamoutO,'logs'))
            @transform(INfiles,suffix('_R1.fastq.gz'),'_pe.bam',output_dir=bamoutO)
            def map_reads(input_files,output_file):
                ii1 = input_files[0]
                ii2 = input_files[1]
                oo = output_file
                from BSmapWGBS import Bism_map_reads
                Bism_map_reads(ii1,ii2,bismpath,bamoutO,crefGpath,args.nthreads,mySession,logger)
elif ( args.aligner=='bwaMeth' and not args.bamdir and not args.methtabdir ):
    if args.convRef:
        if args.trimReads:
            @follows("conv_ref")
            @mkdir(bamoutO,os.path.join(bamoutO,'logs'))
            @transform(trim_reads,suffix('_R1.fastq.gz'),'.sorted.bam',output_dir=bamoutO)
            def map_reads(input_files,output_file):
                ii1 = input_files[0]
                ii2 = input_files[1]
                oo = output_file
                crefG=os.path.join(crefGpath,os.path.basename(refG))               
                from BSmapWGBS import bMeth_map_reads
                bMeth_map_reads(ii1,ii2,oo,bmethpath,sampath,bamoutO,crefG,args.nthreads,mySession,logger)
        else:
            @follows("conv_ref")
            @mkdir(bamoutO,os.path.join(bamoutO,'logs'))
            @transform(INfiles,suffix('_R1.fastq.gz'),'.sorted.bam',output_dir=bamoutO)
            def map_reads(input_files,output_file):
                ii1 = input_files[0]
                ii2 = input_files[1]
                oo = output_file
                crefG=os.path.join(crefGpath,os.path.basename(refG))
                from BSmapWGBS import bMeth_map_reads
                bMeth_map_reads(ii1,ii2,oo,bmethpath,sampath,bamoutO,crefG,args.nthreads,mySession,logger)
    else:
        if args.trimReads:
            @mkdir(bamoutO,os.path.join(bamoutO,'logs'))
            @transform(trim_reads,suffix('_R1.fastq.gz'),'.sorted.bam',output_dir=bamoutO)
            def map_reads(input_files,output_file):
                ii1 = input_files[0]
                ii2 = input_files[1]
                oo = output_file
                from BSmapWGBS import bMeth_map_reads
                bMeth_map_reads(ii1,ii2,oo,bmethpath,sampath,bamoutO,crefG,args.nthreads,mySession,logger)
        else:
            @mkdir(bamoutO,os.path.join(bamoutO,'logs'))
            @transform(INfiles,suffix('_R1.fastq.gz'),'.sorted.bam',output_dir=bamoutO)
            def map_reads(input_files,output_file):
                ii1 = input_files[0]
                ii2 = input_files[1]
                oo = output_file
                from BSmapWGBS import bMeth_map_reads
                bMeth_map_reads(ii1,ii2,oo,bmethpath,sampath,bamoutO,crefG,args.nthreads,mySession,logger)  

################## BAM POSTPROCESSING #########################################################
#bam postprocessing; bwa-meth is sorted and indexed, Bismark and methylCtools require these steps

if ( args.aligner=='methylCtools' and not args.bamdir and not args.methtabdir ):
    @transform(backconvert_bam,suffix('.bconv.s.bam'),'.bconv.s.bam.bai')
    def index_bam(input_file,output_file):
        ii = input_file
        oo = output_file
        from BSmapWGBS import BS_index_bam
        BS_index_bam(ii,sampath,bamoutO,mySession,logger)
    @transform(index_bam,suffix('.bconv.s.bam.bai'),'.PCRrm.bam')
    def PCRdup_rm(input_file,output_file):
        ii=re.sub('.bai','',input_file)
        oo=output_file
        from BSmapWGBS import BS_rm_dupes
        BS_rm_dupes(ii,oo,Picpath,bamoutO,mySession,logger)
        
elif ( args.aligner=='Bismark' and not args.bamdir and not args.methtabdir ):    
    @transform(map_reads,suffix("_pe.bam"),'.sorted.bam')
    def sort_bam(input_file,output_file):
        ii = input_file
        oo = output_file
        from BSmapWGBS import BS_sort_bam
        BS_sort_bam(ii,sampath,bamoutO,args.nthreads,mySession,logger)
        
    @transform(sort_bam,suffix('.sorted.bam'),'.sorted.bam.bai') 
    def index_bam(input_file,output_file):
        ii = input_file
        oo = output_file  
        from BSmapWGBS import BS_index_bam
        BS_index_bam(ii,sampath,bamoutO,mySession,logger)
    @transform(index_bam,suffix('.sorted.bam.bai'),'.RGi.bam') #for compatibility with GATK; depth of coverage
    def add_ReadGroupInfo(input_file,output_file):
        ii=re.sub('.bai','',input_file)
        oo=output_file
        from BSmapWGBS import BS_add_RGi
        BS_add_RGi(ii,oo,Picpath,bamoutO,mySession,logger)    
    @transform(add_ReadGroupInfo,suffix('.RGi.bam'),'.PCRrm.bam')
    def PCRdup_rm(input_file,output_file):
        ii=rinput_file
        oo=output_file
        from BSmapWGBS import BS_rm_dupes
        BS_rm_dupes(ii,oo,Picpath,bamoutO,mySession,logger)  
        
    #remember to symlink sorted.bai to sorted.bam.bai; is this samtools or Picard tools?
elif ( args.aligner=='bwaMeth' and not args.bamdir and not args.methtabdir ):
    @transform(map_reads,suffix('.sorted.bam'),'.sorted.bam.bai') 
    def index_bam(input_file,output_file):
        ii = input_file
        oo = output_file  
        from BSmapWGBS import BS_index_bam
        BS_index_bam(ii,sampath,bamoutO,mySession,logger)
    @transform(index_bam,suffix('.sorted.bam.bai'),'.RGi.bam') #for compatibility with GATK; depth of coverage
    def add_ReadGroupInfo(input_file,output_file):
        ii=re.sub('.bai','',input_file)
        oo=output_file
        from BSmapWGBS import BS_add_RGi
        BS_add_RGi(ii,oo,Picpath,bamoutO,mySession,logger)
        
    @transform(add_ReadGroupInfo,suffix('.RGi.bam'),'.PCRrm.bam')
    def PCRdup_rm(input_file,output_file):
        ii=input_file
        oo=output_file
        from BSmapWGBS import BS_rm_dupes
        BS_rm_dupes(ii,oo,Picpath,bamoutO,mySession,logger)
        
        
######################################################################################################################  
#RUN VARIOUS COVERAGE AND QUALITY METRICS#######################################################################
metout=os.path.join(wdir,'QC_metrics')
auxdir=os.path.join(wdir,'aux_files')

#prepare bed file with random CpGs
if not args.methtabdir :
    @mkdir(auxdir,os.path.join(auxdir,'logs'))
    @transform(refG,suffix('.fa'),'.poz.ran1M.sorted.bed',output_dir=auxdir)
    def get_ran_CG(input_file,output_file):
        ii=input_file
        oo=output_file
        from BSmetricsWGBS import mCT_get_ranCG
        mCT_get_ranCG(ii,oo,refG,auxdir,mCTpath,logger)


#GATK depth of coverage
if ( args.intList and not args.methtabdir and not args.bamdir ): 
    int_dest=[re.sub('.bed','.mean.doc.sample_summary',os.path.basename(x)) for x in args.intList]
    int_dest[:0]=['mean.CG.doc.sample_summary']
    int_dest[:0]=['mean.genome.doc.sample_summary']
    @mkdir(metout,os.path.join(metout,'logs'))
    @follows(get_ran_CG)
    @transform(PCRdup_rm,suffix('PCRrm.bam'),int_dest,output_dir=metout)#
    def depth_of_cov(input_file,output_files):
        ii=input_file
        oos=output_files 
        oos2=[w.replace('.sample_summary', '') for w in oos]
        from BSmetricsWGBS import BS_doc_XT
        BS_doc_XT(ii,oos2,args.intList,auxdir,refG,GATKpath,metout,mySession,logger)
        
        
elif ( args.intList and args.bamdir and not args.methtabdir ): 
    int_dest=[re.sub('.bed','.mean.doc.sample_summary',os.path.basename(x)) for x in args.intList]
    int_dest[:0]=['mean.CG.doc.sample_summary']
    int_dest[:0]=['mean.genome.doc.sample_summary']
    @mkdir(metout,os.path.join(metout,'logs'))
    @follows(get_ran_CG)
    @transform(INfiles,suffix('bam'),int_dest,output_dir=metout)
    def depth_of_cov(input_file,output_files):
        ii=input_file
        oos=output_files 
        oos2=[w.replace('.sample_summary', '') for w in oos]
        from BSmetricsWGBS import BS_doc_XT
        BS_doc_XT(ii,oos2,args.intList,auxdir,refG,GATKpath,metout,mySession,logger)        
                
elif ( not args.intList and not args.methtabdir and not args.bamdir ):
    @mkdir(metout,os.path.join(metout,'logs'))
    @follows(get_ran_CG)
    @transform(PCRdup_rm,suffix('.PCRrm.bam'),['.mean.genome.doc.sample_summary','.mean.CG.doc.sample_summary'],output_dir=metout)#,
    def depth_of_cov(input_file,output_files):
        ii=input_file
        oos=output_files
        oos2=[w.replace('.sample_summary', '') for w in oos]
        from BSmetricsWGBS import BS_doc
        BS_doc(ii,oos2,refG,auxdir,GATKpath,metout,mySession,logger)
        
elif ( args.bamdir and not args.intList and not args.methtabdir ):
    @mkdir(metout,os.path.join(metout,'logs'))
    @follows(get_ran_CG)
    @transform(INfiles,suffix('bam'),['.mean.genome.doc.sample_summary','.mean.CG.doc.sample_summary'],output_dir=metout)#,
    def depth_of_cov(input_file,output_files):
        ii=input_file
        oos=output_files
        oos2=[w.replace('.sample_summary', '') for w in oos]
        from BSmetricsWGBS import BS_doc
        BS_doc(ii,oos2,refG,auxdir,GATKpath,metout,mySession,logger)        
    
#conversion rate: currently from fastq, implement phiX/lambda control!
if ( args.trimReads and not args.methtabdir and not args.bamdir ):
    @mkdir(metout,os.path.join(metout,'logs'))
    @transform(trim_reads,suffix('_R1.fastq.gz'),'.conv.rate.txt',output_dir=metout) #it's ok to have both reads summarized in 1 file
    def conv_rate(input_files,output_file):
        ii1=input_files[0]
        ii1sub=re.sub('_R1.fastq.gz','',ii1)
        oo=output_file
        from BSmetricsWGBS import BS_conv_rate
        BS_conv_rate(ii1sub,oo,metout,mySession,logger)
elif ( not args.trimReads and not args.methtabdir and not args.bamdir ): 
    @mkdir(metout,os.path.join(metout,'logs'))
    @transform(INfiles,suffix('_R1.fastq.gz'),'.conv.rate.txt',output_dir=metout) #it's ok to have both reads summarized in 1 file
    def conv_rate(input_files,output_file):
        ii1=input_files[0][0]
        ii1sub=re.sub('_R1.fastq.gz','',ii1)
        oo=output_file
        from BSmetricsWGBS import BS_conv_rate
        BS_conv_rate(ii1sub,oo,metout,mySession,logger)
    
#methylation bias
if ( not args.methtabdir and not args.bamdir ):
    @mkdir(metout,os.path.join(metout,'logs'))
    @transform(PCRdup_rm,suffix('.PCRrm.bam'),'.Mbias.txt',output_dir=metout)
    def calc_Mbias(input_file,output_file):
        ii=input_file
        oo=output_file
        oos=re.sub('.txt','',oo)
        from BSmetricsWGBS import BS_Mbias
        BS_Mbias(ii,oos,POMpath,refG,metout,args.nthreads,mySession,logger)
        
if ( args.bamdir and not args.methtabdir ):
    @mkdir(metout,os.path.join(metout,'logs'))
    @transform(INfiles,suffix('.bam'),'.Mbias.txt',output_dir=metout)
    def calc_Mbias(input_file,output_file):
        ii=input_file
        oo=output_file
        oos=re.sub('.txt','',oo)
        from BSmetricsWGBS import BS_Mbias
        BS_Mbias(ii,oos,POMpath,refG,metout,args.nthreads,mySession,logger)  
        
#flagstat -> mapping rate

if ( not args.methtabdir and not args.bamdir ):
    @mkdir(metout,os.path.join(metout,'logs'))
    @transform(PCRdup_rm,suffix('.PCRrm.bam'),'.flagstat',output_dir=metout)
    def get_flagstat(input_file,output_file):
        ii=input_file
        oo=output_file
        from BSmetricsWGBS import BS_flagstat
        BS_flagstat(ii,oo,sampath,metout,mySession,logger)
        
if ( args.bamdir and not args.methtabdir ):
    @mkdir(metout,os.path.join(metout,'logs'))
    @transform(INfiles,suffix('.bam'),'.flagstat',output_dir=metout)
    def get_flagstat(input_file,output_file):
        ii=input_file
        oo=output_file
        from BSmetricsWGBS import BS_flagstat
        BS_flagstat(ii,oo,sampath,metout,mySession,logger)  

################################################################################################################## 
###final QC report: only if started from reads or bam files
if ( args.readdir and not args.methtabdir and not args.bamdir ):
    if not os.path.exists(metout):
        os.makedirs(metout)
    os.chdir(metout)
    @merge(input=[depth_of_cov,conv_rate,calc_Mbias,get_flagstat],output='QC_report.pdf')
    def produce_report(input_file,output_file):
        ii=input_file
        oo=output_file
        from BSmetricsWGBS import BS_QC_rep
        BS_QC_rep(Rpath,metout,logger) 
        
if ( args.bamdir and not args.methtabdir and not args.readdir ): 
    if not os.path.exists(metout):
        os.makedirs(metout)
    os.chdir(metout)
    @merge(input=[depth_of_cov,calc_Mbias,get_flagstat],output='QC_report.pdf')
    def produce_report(input_file,output_file):
        ii=input_file
        oo=output_file
        from BSmetricsWGBS import BS_QC_rep
        BS_QC_rep(Rpath,metout,logger)
        
#####EXTRACT METHYLATION COUNTS  ########################################################################
mextout=os.path.join(wdir,'methXT_')+args.methEXT
#auxdir=os.path.join(wdir,'aux_files') #defined previously


if ( args.methEXT=='methylCtools' and not args.methtabdir ):
    @transform(refG,suffix('.fa'),['.poz.P.txt','.poz.M.txt'],output_dir=auxdir)
    def get_poz(input_file,output_files):
        ii=input_file
        pozF=os.path.join(auxdir,re.sub('.fa','.poz.gz',os.path.basename(refG)))
        oo1=output_files[0]
        oo2=output_files[1]
        from BSmethXT_WGBS import mCT_prep_poz
        mCT_prep_poz(refG,mCTpath,tabpath,auxdir,pozF,logger)
    if ( not args.bamdir ):    
        @follows('calc_Mbias',mkdir(mextout),mkdir(os.path.join(mextout,'logs')))
        @transform(PCRdup_rm,suffix('.PCRrm.bam'),'.CG.call.gz',output_dir=mextout)
        def methyl_extract(input_file,output_file):
            ii = input_file
            oo = output_file
            from BSmethXT_WGBS import methXT_mCT
            methXT_mCT(ii,metout,oo,refG,mCTpath,tabpath,mextout,args.mbias_ignore,mySession,logger)
    elif ( args.bamdir ):    
        @follows('calc_Mbias',mkdir(mextout),mkdir(os.path.join(mextout,'logs')))
        @transform(INfiles,suffix('.PCRrm.bam'),'.CG.call.gz',output_dir=mextout)
        def methyl_extract(input_file,output_file):
            ii = input_file
            oo = output_file
            from BSmethXT_WGBS import methXT_mCT
            methXT_mCT(ii,metout,oo,refG,mCTpath,tabpath,mextout,mySession,logger)        
    @follows('get_poz',mkdir(mextout),mkdir(os.path.join(mextout,'logs')))    
    @transform(methyl_extract,suffix('.CG.call.gz'),'.CpG.filt2.bed',output_dir=mextout) 
    def CpG_filt(input_file,output_file):
        ii = input_file
        oo = output_file 
        from BSmethXT_WGBS import filt_mCT
        filt_mCT(ii,refG,mextout,mySession,logger)
        
elif ( args.methEXT=='methylCtools' and  args.methtabdir ): 
    @transform(INfiles,suffix('.CG.call.gz'),'.CpG.filt2.bed',output_dir=mextout) 
    def CpG_filt(input_file,output_file):
        ii = input_file
        oo = output_file 
        from BSmethXT_WGBS import filt_mCT
        filt_mCT(ii,refG,mextout,mySession,logger)
        
elif ( args.methEXT=='MethylDackel' and not args.methtabdir ): 
    if ( not args.bamdir ):
        @follows('calc_Mbias',mkdir(mextout),mkdir(os.path.join(mextout,'logs')))
        @transform(PCRdup_rm,suffix('.PCRrm.bam'),'_CpG.bedGraph',output_dir=mextout)
        def methyl_extract(input_file,output_file):
            ii = input_file
            oo = output_file
            oos=re.sub('_CpG.bedGraph','',oo)
            from BSmethXT_WGBS import methXT_POM
            methXT_POM(ii,metout,oos,refG,POMpath,mextout,args.mbias_ignore,args.nthreads,mySession,logger)
    elif (  args.bamdir ):
        @follows('calc_Mbias',mkdir(mextout),mkdir(os.path.join(mextout,'logs')))
        @transform(INfiles,suffix('.bam'),'_CpG.bedGraph',output_dir=mextout)
        def methyl_extract(input_file,output_file):
            ii = input_file
            oo = output_file
            oos=re.sub('_CpG.bedGraph','',oo)
            from BSmethXT_WGBS import methXT_POM
            methXT_POM(ii,metout,oos,refG,POMpath,mextout,args.mbias_ignore,args.nthreads,mySession,logger)        
    @transform(methyl_extract,suffix('_CpG.bedGraph'),'.CpG.filt2.bed',output_dir=mextout) 
    def CpG_filt(input_file,output_file):
        ii = input_file
        oo = output_file
        from BSmethXT_WGBS import filt_POM
        filt_POM(ii,bedpath,mextout,mySession,logger,args.blackList)
        
elif ( args.methEXT=='MethylDackel' and  args.methtabdir ): 
    @transform(INfiles,suffix('_CpG.bedGraph'),'.CpG.filt2.bed',output_dir=mextout) 
    def CpG_filt(input_file,output_file):
        ii = input_file
        oo = output_file
        from BSmethXT_WGBS import filt_POM
        filt_POM(ii,bedpath,mextout,mySession,args.blackList,logger)
        
elif ( args.methEXT=='BisSNP' and not args.methtabdir ):
    if ( not args.bamdir ):
        @follows('calc_Mbias',mkdir(mextout),mkdir(os.path.join(mextout,'logs')))
        @transform(PCRdup_rm,metout,suffix('.PCRrm.bam'),['.SNP.vcf','.CpG.vcf'],output_dir=mextout)
        def methyl_extract(input_file,output_file):
            ii = input_file
            ooSNP = output_file[0]
            ooCG = output_file[1]
            from BSmethXT_WGBS import methXT_BisSNP
            methXT_BisSNP(ii,ooSNP,ooCG,refG,BisSNPpath,mextout,args.mbias_ignore,mySession,args.nthreads,logger)
    elif (  args.bamdir ):
        @follows('calc_Mbias',mkdir(mextout),mkdir(os.path.join(mextout,'logs')))
        @transform(INfiles,suffix('.PCRrm.bam'),['.SNP.vcf','.CpG.vcf'],output_dir=mextout)
        def methyl_extract(input_file,output_file):
            ii = input_file
            ooSNP = output_file[0]
            ooCG = output_file[1]
            from BSmethXT_WGBS import methXT_BisSNP
            methXT_BisSNP(ii,metout,ooSNP,ooCG,refG,BisSNPpath,mextout,args.mbias_ignore,mySession,args.nthreads,logger)        
    @transform(methyl_extract,suffix('.SNP.vcf'),['.SNP.ro.vcf','.CpG.ro.vcf'],output_dir=mextout)
    def sort_vcf(input_files,output_files):
        ii1=input_files[0]
        ii2=input_files[1]
        oo1 = output_files[0]
        oo2 = output_files[1]
        from BSmethXT_WGBS import vcf_sort
        vcf_sort(ii1,ii2,oo1,oo2,refG,mextout,mySession,logger)
        
    @transform(sort_vcf,suffix('.SNP.ro.vcf'),'.SNP.filt.vcf',output_dir=mextout) 
    def SNP_filt(input_files,output_file):
        ii = input_files[0]
        oo = output_file
        from BSmethXT_WGBS import BisSNP_SNP_filt
        BisSNP_SNP_filt(ii,oo,BisSNPpath,refG,mextout,mySession,logger)
        
    @follows('sort_vcf')    
    @transform(SNP_filt,suffix('.SNP.filt.vcf'),'.CpG.filt.vcf',output_dir=mextout)
    def CpG_m_SNP(input_file,output_file):
        ii1 = input_file
        ii2 = re.sub('.SNP.filt.vcf','.CpG.ro.vcf',ii1)
        oo = output_file 
        from BSmethXT_WGBS import BisSNP_rmSNP
        BisSNP_rmSNP(ii1,ii2,oo,BisSNPpath,refG,mextout,mySession,logger)
        
    @transform(CpG_m_SNP,suffix('.CpG.filt.vcf'),'.CpG.filt.bed',output_dir=mextout)
    def vcf2bed(input_file,output_file):
        ii=input_file
        oo=output_file
        from BSmethXT_WGBS import BisSNP_vcf2bed
        BisSNP_vcf2bed(ii,oo,mextout,mySession,logger)
        
    @transform(vcf2bed,suffix('.CpG.filt.bed'),'.CpG.filt2.bed',output_dir=mextout)  
    def CpG_filt(input_file,output_file):
        ii=input_file
        oo=output_file
        from BSmethXT_WGBS import BisSNP_CpGfilt
        BisSNP_CpGfilt(ii,mextout,mySession,logger)
        
elif ( args.methEXT=='BisSNP' and  args.methtabdir ):  
    @transform(INfiles,suffix('.CpG.filt.bed'),'.CpG.filt2.bed',output_dir=mextout)  
    def CpG_filt(input_file,output_file):
        ii=input_file
        oo=output_file
        from BSmethXT_WGBS import BisSNP_CpGfilt
        BisSNP_CpGfilt(ii,mextout,mySession,logger)
        
        
#########################################################################################################     
####RUN SINGLE CYTOSINE STATISTICS IN R VIA LIMMA#########################################################
if (args.sampleInfo and args.stats=='limma'):
    CpGstat_out=os.path.join(wdir,'singleCpG_stats_limma_')+args.methEXT
    if not os.path.exists(CpGstat_out) or not os.path.exists(os.path.join(CpGstat_out,'logs')):
        os.makedirs(CpGstat_out)
        os.makedirs(os.path.join(CpGstat_out,'logs'))
    os.chdir(CpGstat_out)
    @merge(input=CpG_filt,output=[os.path.join(CpGstat_out,'singleCpG.RData'),os.path.join(CpGstat_out,'limdat.LG.RData'),os.path.join(CpGstat_out,'metilene.IN.txt')])
    def CpG_stats(input_files,output_file):
        ii = input_files[0]
        oo = output_file
        from BSstats_WGBS import single_CpG_limma
        single_CpG_limma(os.path.dirname(ii),args.sampleInfo,CpGstat_out,mySession,logger)
elif (args.sampleInfo and args.stats=='DSS'):
    CpGstat_out=os.path.join(wdir,'singleCpG_stats_DSS_')+args.methEXT
    if not os.path.exists(CpGstat_out) or not os.path.exists(os.path.join(CpGstat_out,'logs')):
        os.makedirs(CpGstat_out)
        os.makedirs(os.path.join(CpGstat_out,'logs'))
    os.chdir(CpGstat_out)
    @merge(input=CpG_filt,output=[os.path.join(CpGstat_out,'singleCpG.RData'),os.path.join(CpGstat_out,'DSS.res.RData')])
    def CpG_stats(input_files,output_file):
        ii = input_files[0]
        oo = output_file
        from BSstats_WGBS import single_CpG_DSS
        single_CpG_DSS(os.path.dirname(ii),args.sampleInfo,CpGstat_out,mySession,logger)        

########################################################################################################
####RUN INTERVAL AGGREGATE STATISTICS IN R ####################################################
if ( args.intList and args.stats=='limma' and args.sampleInfo ):
    int_dest=[re.sub('.bed','.aggCpG.RData',os.path.basename(x)) for x in args.intList]
    intStat_out=os.path.join(wdir,args.methEXT)+'_aggregate_stats_limma'
    @mkdir(intStat_out,os.path.join(intStat_out,'logs'))
    @transform(CpG_stats,suffix('.RData'),int_dest,output_dir=intStat_out)
    def intAgg_stats(input_files,output_files):
        ii = os.path.join(CpGstat_out,input_files[1])
        oo = output_files
        from BSstats_WGBS import int_stats_limma
        int_stats_limma(ii,args.intList,args.sampleInfo,intStat_out,mySession,logger)

elif ( args.intList and args.stats=='DSS' and args.sampleInfo ):
    int_dest=[re.sub('.bed','.aggCpG.RData',os.path.basename(x)) for x in args.intList]
    intStat_out=os.path.join(wdir,args.methEXT)+'_aggregate_stats_DSS'
    @mkdir(intStat_out,os.path.join(intStat_out,'logs'))
    @transform(CpG_stats,suffix('.RData'),int_dest,output_dir=intStat_out)
    def intAgg_stats(input_files,output_files):
        ii = os.path.join(CpGstat_out,input_files[1])
        oo = output_files
        from BSstats_WGBS import int_stats_DSS
        int_stats_DSS(ii,mextout,args.intList,args.sampleInfo,intStat_out,mySession,logger)

#####################################################################################################
####RUN DMR calling ################################################################################
if args.DMRpg=='methylSeekR' and args.sampleInfo : 
    msRout=os.path.join(wdir,'methylSeekR_out')
    @mkdir(msRout,os.path.join(msRout,'logs'))
    @transform(CpG_stats,suffix('.RData'),'methylSeekR.RData',output_dir=msRout)
    def run_mseekR (input_files,output_files):
        ii = open(input_files[1])
        oo = open(output_file, "w")
        print('This is a placeholder. MethylSeekR is currently not implemented.')
elif (args.DMRpg=='metilene' and args.stats=='limma' and args.sampleInfo ):
    DMRout=os.path.join(wdir,'metilene_out')
    @mkdir(DMRout,os.path.join(DMRout,'logs'))
    @transform(CpG_stats,suffix('.RData'),'.metilene.bed',output_dir=DMRout)
    def run_metilene (input_files,output_file):
        ii = input_files[2]
        oo = output_file
        from BS_DMR_WGBS import DMR_metilene
        DMR_metilene(ii,args.sampleInfo,oo,args.nthreads,metipath,mySession,logger)
    os.chdir(CpGstat_out)
    @merge(input=[run_metilene,CpG_stats],output=os.path.join(DMRout,"singleCpG.metilene.limma.bed"))
    def cleanup_metilene (input_files,output_file):
        ii1 = input_files[0]
        ii2 = input_files[1][1]
        oo = output_file
        from BS_DMR_WGBS import clean_up_metilene
        clean_up_metilene(ii1,ii2,args.sampleInfo,DMRout,mySession,logger)
###################################################################################################

#####main

if __name__ == '__main__':
    with open(os.path.join(wdir,"pipelineGraph.png"),'w') as pipeGraph:
        pipeline_printout_graph(stream=pipeGraph,output_format='png',pipeline_name='WGBS',target_tasks=args.target_tasks)
    with open (os.path.join(wdir,"pipelinePrint.txt"),'w') as pipePrint:
        pipeline_printout(verbose_abbreviated_path=0,output_stream=pipePrint,target_tasks=args.target_tasks)    

    pipeline_run(touch_files_only=args.touchOnly,multiprocess=args.bsize,target_tasks=args.target_tasks,forcedtorun_tasks=args.forcedtorun_tasks,logger=logger)
    mySession.exit()
