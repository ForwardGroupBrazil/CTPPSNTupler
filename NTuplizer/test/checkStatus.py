#!/usr/bin/env python

################
# Merger Files #
################

import time, os

os.system("crab status -d crab_ElectronB/crab_ElectronB/ > statusJob.txt")
os.system("crab status -d crab_ElectronC/crab_ElectronC/ >> statusJob.txt")
os.system("crab status -d crab_ElectronD/crab_ElectronD/ >> statusJob.txt")
os.system("crab status -d crab_ElectronE/crab_ElectronE/ >> statusJob.txt")
os.system("crab status -d crab_ElectronF/crab_ElectronF/ >> statusJob.txt")
os.system("crab status -d crab_MuonB/crab_MuonB/ >> statusJob.txt")
os.system("crab status -d crab_MuonC/crab_MuonC/ >> statusJob.txt")
os.system("crab status -d crab_MuonD/crab_MuonD/ >> statusJob.txt")
os.system("crab status -d crab_MuonE/crab_MuonE/ >> statusJob.txt")
os.system("crab status -d crab_MuonF/crab_MuonF/ >> statusJob.txt")
