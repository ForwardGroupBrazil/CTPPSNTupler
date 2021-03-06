#!/usr/bin/env python

##############
# Multi CRAB #
##############

import time 

if __name__ == '__main__':
  from CRABAPI.RawCommand import crabCommand

def submit(config):
  res = crabCommand('submit', config = config)

from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

dataset = {
  'MuonB' : '/DoubleMuon/Run2017B-17Nov2017-v1/AOD',
  'MuonC' : '/DoubleMuon/Run2017C-17Nov2017-v1/AOD',
  'MuonD' : '/DoubleMuon/Run2017D-17Nov2017-v1/AOD',
  'MuonE' : '/DoubleMuon/Run2017E-17Nov2017-v1/AOD',
  'MuonF' : '/DoubleMuon/Run2017F-17Nov2017-v1/AOD',
  'ElectronB' : '/DoubleEG/Run2017B-17Nov2017-v1/AOD',
  'ElectronC' : '/DoubleEG/Run2017C-17Nov2017-v1/AOD',
  'ElectronD' : '/DoubleEG/Run2017D-17Nov2017-v1/AOD',
  'ElectronE' : '/DoubleEG/Run2017E-17Nov2017-v1/AOD',
  'ElectronF' : '/DoubleEG/Run2017F-17Nov2017-v1/AOD',
}

mode = "Muon" 
filesPerJob = 1

config.General.transferLogs = True
config.General.transferOutputs = True
config.JobType.pluginName = 'Analysis'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.publication = False
config.Site.storageSite = 'T2_CH_CERN'
config.JobType.psetName = 'RunGammaGammaLeptonLepton_cfg.py'
#config.Data.ignoreLocality = True 
#config.JobType.disableAutomaticOutputCollection = False

def doSubmit(listOfSamples):
  for sample in listOfSamples:
    config.JobType.outputFiles = ['output.root']
    config.General.workArea = 'crab_'+ sample
    config.General.requestName = sample
    config.Data.inputDataset = dataset[sample]
    config.Data.unitsPerJob = filesPerJob
    config.Data.lumiMask = 'combined_RPIN_CMS_2017.json'
    config.Data.outputDatasetTag = sample
    config.Data.outLFNDirBase = '/store/group/phys_pps/MissingMassSearch/' + sample
    config.Site.storageSite = 'T2_CH_CERN'
    config.JobType.pyCfgParams = ["Lepton="+mode]
    submit(config)

# MuonB
mode = "Muon"
filesPerJob = 200
listOfSamples = ['MuonB']
doSubmit(listOfSamples)

# MuonC
mode = "Muon"
filesPerJob = 200
listOfSamples = ['MuonC']
doSubmit(listOfSamples)

# MuonD
mode = "Muon"
filesPerJob = 200
listOfSamples = ['MuonD']
doSubmit(listOfSamples)

# MuonE
mode = "Muon"
filesPerJob = 200
listOfSamples = ['MuonE']
doSubmit(listOfSamples)

# MuonF
mode = "Muon"
filesPerJob = 200
listOfSamples = ['MuonF']
doSubmit(listOfSamples)

# ElectronB
mode = "Electron"
filesPerJob = 200
listOfSamples = ['ElectronB']
doSubmit(listOfSamples)

# ElectronC
mode = "Electron"
filesPerJob = 200
listOfSamples = ['ElectronC']
doSubmit(listOfSamples)

# ElectronD
mode = "Electron"
filesPerJob = 200
listOfSamples = ['ElectronD']
doSubmit(listOfSamples)

# ElectronE
mode = "Electron"
filesPerJob = 200
listOfSamples = ['ElectronE']
doSubmit(listOfSamples)

# ElectronF
mode = "Electron"
filesPerJob = 200
listOfSamples = ['ElectronF']
doSubmit(listOfSamples)

