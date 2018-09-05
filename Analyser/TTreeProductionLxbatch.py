#!/usr/bin/env python
import os, re
import commands
import math, time
import sys

queue = "1nw" # give bsub queue -- 8nm (8 minutes), 1nh (1 hour), 8nh, 1nd (1day), 2nd, 1nw (1 week), 2nw 

command = [	#'./DarkMatterSearchLeadingMuonsJets -1 dataMuMu_2016 B 2016',
		#'./DarkMatterSearchLeadingMuonsJets -1 dataMuMu_2016 C 2016',
		#'./DarkMatterSearchLeadingMuonsJets -1 dataMuMu_2016 G 2016',
		'./DarkMatterSearchLeadingMuonsJets -1 dataMuMu_2017 B 2017',
		'./DarkMatterSearchLeadingMuonsJets -1 dataMuMu_2017 C 2017',
		'./DarkMatterSearchLeadingMuonsJets -1 dataMuMu_2017 D 2017',
		'./DarkMatterSearchLeadingMuonsJets -1 dataMuMu_2017 E 2017',
		'./DarkMatterSearchLeadingMuonsJets -1 dataMuMu_2017 F 2017',
		#'./DarkMatterSearchLeadingElectronsJets -1 dataEE_2016 B 2016',
		#'./DarkMatterSearchLeadingElectronsJets -1 dataEE_2016 C 2016',
		#'./DarkMatterSearchLeadingElectronsJets -1 dataEE_2016 G 2016',
		'./DarkMatterSearchLeadingElectronsJets -1 dataEE_2017 B 2017',
		'./DarkMatterSearchLeadingElectronsJets -1 dataEE_2017 C 2017',
		'./DarkMatterSearchLeadingElectronsJets -1 dataEE_2017 D 2017',
		'./DarkMatterSearchLeadingElectronsJets -1 dataEE_2017 E 2017',
		'./DarkMatterSearchLeadingElectronsJets -1 dataEE_2017 F 2017'
		#'./DarkMatterSearchLeadingMuonsFinn -1 dataMuMu_2017_finn C 2017',
		#'./DarkMatterSearchLeadingMuonsFinn -1 dataMuMu_2017_finn D 2017',
		#'./DarkMatterSearchLeadingMuonsFinn -1 dataMuMu_2017_finn E 2017',
		#'./DarkMatterSearchLeadingMuonsFinn -1 dataMuMu_2017_finn F 2017'	
	  ]	

print '\nSending Lxbatch jobs to produce NTuples for Missing Mass CTPPS Analysis\n\n'
path = os.getcwd()

i = 0
while i < len(command):
	with open('job.sh', 'w') as fout:
		fout.write("#!/bin/sh\n")
		fout.write("echo\n")
		fout.write("echo 'START---------------'\n")
		fout.write("\n")
		fout.write("source /afs/cern.ch/sw/lcg/external/gcc/4.9/x86_64-slc6-gcc49-opt/setup.sh\n")
		fout.write("source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.06.06/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh\n")
		fout.write("cd "+str(path)+"\n")
		fout.write(command[i]+"\n")
		fout.write("echo 'STOP---------------'\n")
		fout.write("echo\n")
	os.system("chmod 755 job.sh")
	os.system("bsub -q "+queue+" -o logs2nd -J job2nd_"+str(i)+" < job.sh")
	print "job nr " + str(i) + " submitted"
	i += 1

print
print "your jobs:"
os.system("bjobs")
print
print 'END'
print
