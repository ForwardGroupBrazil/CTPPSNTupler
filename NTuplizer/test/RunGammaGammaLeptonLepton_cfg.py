#flake8: noqa

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

# Setting Input Parameters from Line Command
options = VarParsing ('analysis')
options.register('Lepton','Muon',VarParsing.multiplicity.singleton, VarParsing.varType.string,"Option to Run: Muon, ElectronMuon, Electron")
options.parseArguments()

process = cms.Process("ggll")

runOnMC = False
useAOD = True # AOD or MiniAOD?

if options.Lepton == "Muon":
  print("")
  print("#########")
  print("Data Muon")
  print("#########")
  print("")
  triggerlist = 'HLT_IsoMu27_v*', 'HLT_DoubleMu43NoFiltersNoVtx_v*', 'HLT_DoubleMu48NoFiltersNoVtx_v*'
 
if options.Lepton == "Electron":
  print("")
  print("#############")
  print("Data Electron")
  print("#############")
  print("")
  triggerlist = 'HLT_DoubleEle33_CaloIdL_MW_v*','HLT_Ele27_WPTight_Gsf_v*'

#########################
#    General options    #
#########################

process.load("FWCore.MessageService.MessageLogger_cfi")
process.options   = cms.untracked.PSet(
    #wantSummary = cms.untracked.bool(True),
    #SkipEvent = cms.untracked.vstring('ProductNotFound'),
    allowUnscheduled = cms.untracked.bool(True),
)

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#########################
#      Input files      #
#########################

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#'/store/data/Run2017B/DoubleMuon/AOD/17Nov2017-v1/30000/1A3141E3-20D6-E711-9BE7-02163E01A2E8.root',
'/store/data/Run2017C/DoubleEG/AOD/17Nov2017-v1/70000/FED9F7EB-93E4-E711-9116-FA163EE2AD3F.root'
    ),
    #firstEvent = cms.untracked.uint32(0)
)

#########################
#        Triggers       #
#########################
process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.HLTFilter_cfi")
process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltFilter.HLTPaths = (triggerlist)

#########################
#      Preskimming      #
#########################
process.load("Configuration.StandardSequences.GeometryDB_cff") ## FIXME need to ensure that this is the good one
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data')

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")


#########################
#     PAT-ification     #
#########################
## Look at https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATTools#Core_Tools for more information

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('PATuple.root'),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_offline*PrimaryVertices*_*_*',
        'keep *_selectedPatMuons*_*_*',
        'keep *_*lectron*_*_*',
        'keep *_selectedPatElectrons*_*_*',
        'keep *_selectedPat*Photons*_*_*',
        'keep *_selectedPatJets*_*_*',
        'keep *_*MET*_*_*',
        'keep *_*particleFlow*_*_*',
    ),
)
from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask
patAlgosToolsTask = getPatAlgosToolsTask(process)

process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
patAlgosToolsTask.add(process.patCandidatesTask)

process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")
patAlgosToolsTask.add(process.selectedPatCandidatesTask)

from PhysicsTools.PatAlgos.tools.coreTools import runOnData
if not runOnMC:
    runOnData( process )

## add trigger information to the configuration
#from PhysicsTools.PatAlgos.tools.trigTools import *
#switchOnTrigger( process )

#########################
#      Electron ID      #
#########################

#from PhysicsTools.SelectorUtils.tools.vid_id_tools import switchOnVIDElectronIdProducer, setupVIDElectronSelection, setupAllVIDIdsInModule, DataFormat
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

switchOnVIDElectronIdProducer(process, DataFormat.AOD)
#setupAllVIDIdsInModule(process, 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff', setupVIDElectronSelection)
setupAllVIDIdsInModule(process, 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff', setupVIDElectronSelection)

#########################
#       Photon ID       #
#########################

switchOnVIDPhotonIdProducer(process, DataFormat.AOD)
setupAllVIDIdsInModule(process, 'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring16_nonTrig_V1_cff', setupVIDPhotonSelection)

#########################
#       Analysis        #
#########################

process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.GammaGammaLL_cfi")

process.ggll_aod.triggersList = process.hltFilter.HLTPaths
process.ggll_aod.leptonsType = cms.string(options.Lepton)
#process.ggll_aod.leptonsType = cms.string('Muon')
#process.ggll_aod.leptonsType = cms.string('ElectronMuon')
#process.ggll_aod.leptonsType = cms.string('Electron')
process.ggll_aod.runOnMC = cms.bool(runOnMC)
process.ggll_aod.fetchProtons = cms.bool(True)

# E/gamma identification
process.ggll_aod.eleIdLabels = cms.PSet(
    mediumLabel = cms.InputTag('mvaEleID-Spring16-GeneralPurpose-V1-wp90'),
    tightLabel = cms.InputTag('mvaEleID-Spring16-GeneralPurpose-V1-wp80'),
)
process.ggll_aod.phoIdLabels = cms.PSet(
    mediumLabel = cms.InputTag('mvaPhoID-Spring16-nonTrig-V1-wp90'),
    tightLabel = cms.InputTag('mvaPhoID-Spring16-nonTrig-V1-wp80'),
)
#process.ggll_aod.eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp90")
#process.ggll_aod.eleTightIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp80")
#process.ggll_aod.phoMediumIdMap = cms.InputTag("egmPhotonIDs:mvaPhoID-Spring16-nonTrig-V1-wp90")
#process.ggll_aod.phoTightIdMap = cms.InputTag("egmPhotonIDs:mvaPhoID-Spring16-nonTrig-V1-wp80")

# prepare the output file
process.TFileService = cms.Service('TFileService',
    fileName = cms.string('output.root'),
    closeFileFast = cms.untracked.bool(True)
)

process.p = cms.Path(
    process.hltFilter*
    process.egmPhotonIDSequence*
    process.egmGsfElectronIDSequence*
    process.ggll_aod
)

#process.outpath = cms.EndPath(process.out, patAlgosToolsTask)
process.outpath = cms.EndPath(patAlgosToolsTask)

