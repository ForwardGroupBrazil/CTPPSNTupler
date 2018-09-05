#!/usr/bin/env python
import os, re
import commands
import math, time
import sys

queue = "1nw" # give bsub queue -- 8nm (8 minutes), 1nh (1 hour), 8nh, 1nd (1day), 2nd, 1nw (1 week), 2nw 

command = [	##'hadd /eos/cms/store/group/phys_pps/MissingMassSearch/NTuples/ElectronB.root /eos/cms/store/group/phys_pps/MissingMassSearch/ElectronB/*/*/*/*/*.root',
		#'hadd /eos/cms/store/group/phys_pps/MissingMassSearch/NTuples/ElectronC.root /eos/cms/store/group/phys_pps/MissingMassSearch/ElectronC/*/*/*/*/*.root',
		##'hadd /eos/cms/store/group/phys_pps/MissingMassSearch/NTuples/ElectronD.root /eos/cms/store/group/phys_pps/MissingMassSearch/ElectronD/*/*/*/*/*.root',
		#'hadd /eos/cms/store/group/phys_pps/MissingMassSearch/NTuples/ElectronE.root /eos/cms/store/group/phys_pps/MissingMassSearch/ElectronE/*/*/*/*/*.root',
		##'hadd /eos/cms/store/group/phys_pps/MissingMassSearch/NTuples/ElectronF.root /eos/cms/store/group/phys_pps/MissingMassSearch/ElectronF/*/*/*/*/*.root',
		##'hadd /eos/cms/store/group/phys_pps/MissingMassSearch/NTuples/MuonB.root /eos/cms/store/group/phys_pps/MissingMassSearch/MuonB/*/*/*/*/*.root',
		'hadd /eos/cms/store/group/phys_pps/MissingMassSearch/NTuples/MuonC.root /eos/cms/store/group/phys_pps/MissingMassSearch/MuonC/*/*/*/*/*.root',
		#'hadd /eos/cms/store/group/phys_pps/MissingMassSearch/NTuples/MuonD.root /eos/cms/store/group/phys_pps/MissingMassSearch/MuonD/*/*/*/*/*.root',
		#'hadd /eos/cms/store/group/phys_pps/MissingMassSearch/NTuples/MuonE.root /eos/cms/store/group/phys_pps/MissingMassSearch/MuonE/*/*/*/*/*.root',
		##'hadd /eos/cms/store/group/phys_pps/MissingMassSearch/NTuples/MuonF.root /eos/cms/store/group/phys_pps/MissingMassSearch/MuonF/*/*/*/*/*.root'
	  ]	

print '\nSending Lxbatch jobs to merge NTuples for Missing Mass CTPPS Analysis\n\n'
path = os.getcwd()

i = 0
while i < len(command):
	with open('job.csh', 'w') as fout:
		fout.write("#!/bin/csh\n")
		fout.write("echo\n")
		fout.write("echo 'START---------------'\n")
		fout.write("cd "+str(path)+"\n")
		fout.write("eval `scramv1 runtime -csh`\n")
		fout.write(command[i]+"\n")
		fout.write("echo 'STOP---------------'\n")
		fout.write("echo\n")
	os.system("chmod 755 job.csh")
	os.system("bsub -q "+queue+" -o logsmerge -J jobmerge_"+str(i)+" < job.csh")
	print "job nr " + str(i) + " submitted"
	i += 1

print
print "your jobs:"
os.system("bjobs")
print
print 'END'
print
