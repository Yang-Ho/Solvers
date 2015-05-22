import os
import time
import glob
import importlib
import sys

thisFile = "all_python"

global all_info
global all_valu

all_info = {}
all_valu = {}

workDir = os.getcwd()
thisDir = os.getcwd() 

print ".. importing files from this SANDBOX directory \n{}".format(thisDir)

#sys.path.extend(["./xLib","../../../xBed/xLib"])

sandboxPath = os.path.dirname(thisDir) 
sandboxName = os.path.basename(sandboxPath)
infoVariablesFile = sandboxName + ".info_variables.txt"
infoVariablesFile = "/".join([sandboxPath,"xLib",infoVariablesFile])

if os.path.isfile(infoVariablesFile):
    """
    all_info["sandboxName"] = sandboxName
    all_info["sandboxPath"] = sandboxPath
    all_info["infoVariablesFile"] = infoVariablesFile
    """
    print ".. registered infoVariablesFile as {}".format(infoVariablesFile)
else:
    # Error
    print "\n".join([
        "",
        "ERROR from {}".format(thisFile),
        ".. missing the infoVariablesFile as {}".format(infoVariablesFile)
        ])

thisDir = "../../../xBed/xLib" 
os.chdir(thisDir)
thisDir = os.getcwd()
sys.path.extend([thisDir])
os.chdir(workDir)
print ".. returning to directory \n{}".format(workDir)
print time.strftime("%a %b %d %H:%M:%S %Z %Y")

import P
import core
import util
#import P.lop
#import P.coord
