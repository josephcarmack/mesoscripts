'''Script for running peso postprocessors on all the simulations generated by
a parstubuilder parametric study.'''

# ------------------------------------------------------
# main program
# ------------------------------------------------------

import os
import errno
import subprocess as sp

studyDir = '/absolute/path/to/parstubuilder/parametric/study/directory'

cwd = os.getcwd()

# get list of simulation directories from parStudyInfo.txt file (generated by parstubuilder)
simList = []
with open(studyDir+'/parStudyInfo.txt') as fin:
    for line in fin:
        if 'Unique parameter set' in line:
            break
    for line in fin:
        simList.append(line.rstrip())
fin.close()

# run peso on each simulation directory
cmd = ['/home/joseph/projects/peso/bin/peso'] # modify this path to where you are storing the peso executable
for pth in sorted(simList):
    # copy peso input file to simulation directory
    simPath = studyDir+'/'+pth
    os.system('cp inputPeso.dat '+simPath)
    # change to simulation directory
    os.chdir(simPath)
    # run peso
    print("processing " + simPath)
    flag = sp.check_output(cmd)