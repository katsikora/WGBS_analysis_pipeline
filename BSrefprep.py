#####IMPORTS####################################################

import os
import subprocess
import re


#######DEFINITIONS########
def methCT_refprep(mCTpath,tabpath,bwapath,refpath,outdir,logobject):
    os.chdir(os.path.dirname(refpath))
    outFile=os.path.join(outdir,re.sub('.fa','.conv.fa',os.path.basename(refpath)))
    logobject.info('Preparing reference genome under : ' + os.getcwd())
    logobject.info('CT converting the reference')
    cmd_faconv=[ os.path.join(mCTpath,'methylCtools') + ' faconv ' + refpath + ' ' + outFile]
    logobject.debug(subprocess.check_output(cmd_faconv,shell=True))
    logobject.info('CT converting complete')
    logobject.info('BWA indexing the reference')
    cmd_bwaInd=[ os.path.join(bwapath,'bwa') + ' index ' + outFile]
    logobject.debug(subprocess.check_output(cmd_bwaInd,shell=True))
    logobject.info('BWA indexing complete')
    return


def Bism_refprep(bismpath,BTpath,refpath,outdir,logobject):
    os.chdir(os.path.dirname(refpath))
    logobject.info('Indexing reference genome under : ' + os.getcwd())
    cmd=[os.path.join(bismpath,'bismark_genome_preparation') + ' --path_to_bowtie '+ BTpath + ' ' + outdir ]
    logobject.debug(subprocess.check_output(cmd,shell=True))
    logobject.info('Reference indexing complete')
    return

def bmeth_refprep(bmethpath,refpath,outdir,logobject):
    os.chdir(os.path.dirname(refpath))
    logobject.info('Indexing reference genome under : ' + os.getcwd())
    cmd=['ln -s '+ refpath + ' '+ os.path.join(outdir,os.path.basename(refpath))+'; python ' + os.path.join(bmethpath,'bwameth.py') + ' index ' + os.path.join(outdir,os.path.basename(refpath)) ]
    logobject.debug(subprocess.check_output(cmd,shell=True,))
    logobject.info('Reference indexing complete')
    return
