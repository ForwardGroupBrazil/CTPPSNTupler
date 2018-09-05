#!/bin/csh
echo
echo 'START---------------'
cd /afs/cern.ch/work/d/dmf/private/CMSPhysicsAnalysis/DarkMatterSearch/NTuplizer/CMSSW_9_4_0/src/DiffractiveForwardAnalysis/GammaGammaLeptonLepton/test
eval `scramv1 runtime -csh`
hadd /eos/cms/store/group/phys_pps/MissingMassSearch/NTuples/MuonC.root /eos/cms/store/group/phys_pps/MissingMassSearch/MuonC/*/*/*/*/*.root
echo 'STOP---------------'
echo
