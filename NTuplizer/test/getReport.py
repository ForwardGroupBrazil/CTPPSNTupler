#!/usr/bin/env python

################
# Merger Files #
################

import time, os

os.system("crab report -d crab_ElectronB/crab_ElectronB/ > reportJob.txt")
os.system("crab report -d crab_ElectronC/crab_ElectronC/ >> reportJob.txt")
os.system("crab report -d crab_ElectronD/crab_ElectronD/ >> reportJob.txt")
os.system("crab report -d crab_ElectronE/crab_ElectronE/ >> reportJob.txt")
os.system("crab report -d crab_ElectronF/crab_ElectronF/ >> reportJob.txt")
os.system("crab report -d crab_MuonB/crab_MuonB/ >> reportJob.txt")
os.system("crab report -d crab_MuonC/crab_MuonC/ >> reportJob.txt")
os.system("crab report -d crab_MuonD/crab_MuonD/ >> reportJob.txt")
os.system("crab report -d crab_MuonE/crab_MuonE/ >> reportJob.txt")
os.system("crab report -d crab_MuonF/crab_MuonF/ >> reportJob.txt")
