#####IMPORTS####################################################

import os
import subprocess
import re


#######DEFINITIONS########
def bmeth_refprep(bmethpath,refpath,outdir,logobject):
    os.chdir(os.path.dirname(refpath))
    logobject.info('Indexing reference genome under : ' + os.getcwd())
    cmd=['ln -s '+ refpath + ' '+ os.path.join(outdir,os.path.basename(refpath))+'; python ' + os.path.join(bmethpath,'bwameth.py') + ' index ' + os.path.join(outdir,os.path.basename(refpath)) ]
    logobject.debug(subprocess.check_output(cmd,shell=True,))
    logobject.info('Reference indexing complete')
    return
