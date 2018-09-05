#!/usr/bin/env python

################
# Merger Files #
################

import time, os

os.system("cp crab_ElectronB/crab_ElectronB/results/processedLumis.json /eos/cms/store/group/phys_pps/MissingMassSearch/NTuples/Luminosity/ElectronB_processed_lumis.json")
os.system("cp crab_ElectronC/crab_ElectronC/results/processedLumis.json /eos/cms/store/group/phys_pps/MissingMassSearch/NTuples/Luminosity/ElectronC_processed_lumis.json")
os.system("cp crab_ElectronD/crab_ElectronD/results/processedLumis.json /eos/cms/store/group/phys_pps/MissingMassSearch/NTuples/Luminosity/ElectronD_processed_lumis.json")
os.system("cp crab_ElectronE/crab_ElectronE/results/processedLumis.json /eos/cms/store/group/phys_pps/MissingMassSearch/NTuples/Luminosity/ElectronE_processed_lumis.json")
os.system("cp crab_ElectronF/crab_ElectronF/results/processedLumis.json /eos/cms/store/group/phys_pps/MissingMassSearch/NTuples/Luminosity/ElectronF_processed_lumis.json")

os.system("cp crab_MuonB/crab_MuonB/results/processedLumis.json /eos/cms/store/group/phys_pps/MissingMassSearch/NTuples/Luminosity/MuonB_processed_lumis.json")
os.system("cp crab_MuonC/crab_MuonC/results/processedLumis.json /eos/cms/store/group/phys_pps/MissingMassSearch/NTuples/Luminosity/MuonC_processed_lumis.json")
os.system("cp crab_MuonD/crab_MuonD/results/processedLumis.json /eos/cms/store/group/phys_pps/MissingMassSearch/NTuples/Luminosity/MuonD_processed_lumis.json")
os.system("cp crab_MuonE/crab_MuonE/results/processedLumis.json /eos/cms/store/group/phys_pps/MissingMassSearch/NTuples/Luminosity/MuonE_processed_lumis.json")
os.system("cp crab_MuonF/crab_MuonF/results/processedLumis.json /eos/cms/store/group/phys_pps/MissingMassSearch/NTuples/Luminosity/MuonF_processed_lumis.json")

os.system("cp reportJob.txt /eos/cms/store/group/phys_pps/MissingMassSearch/NTuples/Luminosity/.")


