#!/usr/bin/env python

################
# Merger Files #
################

import time, os

os.system("crab getoutput -d crab_ElectronB/crab_ElectronB/ > getoutputJob.txt")
os.system("crab getoutput -d crab_ElectronC/crab_ElectronC/ >> getoutputJob.txt")
os.system("crab getoutput -d crab_ElectronD/crab_ElectronD/ >> getoutputJob.txt")
os.system("crab getoutput -d crab_ElectronE/crab_ElectronE/ >> getoutputJob.txt")
os.system("crab getoutput -d crab_ElectronF/crab_ElectronF/ >> getoutputJob.txt")
os.system("crab getoutput -d crab_MuonB/crab_MuonB/ >> getoutputJob.txt")
os.system("crab getoutput -d crab_MuonC/crab_MuonC/ >> getoutputJob.txt")
os.system("crab getoutput -d crab_MuonD/crab_MuonD/ >> getoutputJob.txt")
os.system("crab getoutput -d crab_MuonE/crab_MuonE/ >> getoutputJob.txt")
os.system("crab getoutput -d crab_MuonF/crab_MuonF/ >> getoutputJob.txt")
