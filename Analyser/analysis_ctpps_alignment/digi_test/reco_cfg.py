import FWCore.ParameterSet.Config as cms
import fnmatch

process = cms.Process("Reco")

# minimum of logs
process.MessageLogger = cms.Service("MessageLogger",
    statistics = cms.untracked.vstring(),
    destinations = cms.untracked.vstring('cout'),
    cout = cms.untracked.PSet(
        threshold = cms.untracked.string('WARNING')
    )
)

# raw data source
process.load("TotemRawData.Readers.TotemStandaloneRawDataSource_cfi")
process.source.fileNames = cms.untracked.vstring(
		"root://eostotem.cern.ch//eos/totem/data/rawdata/2016/run_10077_EVB11_1.000.srs",
#		"root://eostotem.cern.ch//eos/totem/data/rawdata/2016/run_10077_EVB11_1.001.srs",
#		"root://eostotem.cern.ch//eos/totem/data/rawdata/2016/run_10077_EVB11_1.002.srs",
#		"root://eostotem.cern.ch//eos/totem/data/rawdata/2016/run_10077_EVB11_1.003.srs",
#		"root://eostotem.cern.ch//eos/totem/data/rawdata/2016/run_10077_EVB11_1.004.srs",
#		"root://eostotem.cern.ch//eos/totem/data/rawdata/2016/run_10077_EVB11_1.005.srs",
#		"root://eostotem.cern.ch//eos/totem/data/rawdata/2016/run_10077_EVB11_1.006.srs",
#		"root://eostotem.cern.ch//eos/totem/data/rawdata/2016/run_10077_EVB11_1.007.srs",
#		"root://eostotem.cern.ch//eos/totem/data/rawdata/2016/run_10077_EVB11_1.008.srs",
#		"root://eostotem.cern.ch//eos/totem/data/rawdata/2016/run_10077_EVB11_1.009.srs"
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# raw-to-digi conversion
process.load("CondFormats.TotemReadoutObjects.TotemDAQMappingESSourceXML_cfi")
process.TotemDAQMappingESSourceXML.mappingFileNames.append("CondFormats/TotemReadoutObjects/xml/totem_rp_2016_210_mapping.xml")

#process.load("EventFilter.TotemRawToDigi.totemTriggerRawToDigi_cfi")
#process.totemTriggerRawToDigi.rawDataTag = cms.InputTag("source")
#process.totemTriggerRawToDigi.fedId = 12345

process.load("EventFilter.TotemRawToDigi.totemRPRawToDigi_cfi")
process.totemRPRawToDigi.rawDataTag = cms.InputTag("source")
process.totemRPRawToDigi.fedIds = cms.vuint32(437, 438, 445, 446)

# local RP reconstruction chain with standard settings
process.load("RecoCTPPS.Configuration.recoCTPPS_cff")

for item in fnmatch.filter(process.XMLIdealGeometryESSource_CTPPS.geomXMLFiles, "*/RP_Dist_Beam_Cent.xml"):
    process.XMLIdealGeometryESSource_CTPPS.geomXMLFiles.remove(item)
process.XMLIdealGeometryESSource_CTPPS.geomXMLFiles.append("Geometry/VeryForwardData/data/2016_04_20_fill4828/RP_Dist_Beam_Cent.xml")

process.load("Geometry.VeryForwardGeometryBuilder.TotemRPIncludeAlignments_cfi")
process.TotemRPIncludeAlignments.RealFiles = cms.vstring(
    "/afs/cern.ch/work/j/jkaspar/software/offline/704/src/TotemAlignment/RPData/LHC/2016_04_20_fill4828/version1/sr+el/45.xml",
    "/afs/cern.ch/work/j/jkaspar/software/offline/704/src/TotemAlignment/RPData/LHC/2016_04_20_fill4828/version1/sr+el/56.xml"
)

process.p = cms.Path(
    #process.totemTriggerRawToDigi *
    process.totemRPRawToDigi *
    process.recoCTPPS
)

# output configuration
process.output = cms.OutputModule("PoolOutputModule",
  fileName = cms.untracked.string("reco_test.root"),
  outputCommands = cms.untracked.vstring(
    'drop *_*_*_*',
    'keep TotemFEDInfos_totemRPRawToDigi_*_*',
    'keep TotemTriggerCounters_totemTriggerRawToDigi_*_*',
    'keep TotemRPDigiedmDetSetVector_totemRPRawToDigi_*_*',
    'keep TotemVFATStatusedmDetSetVector_totemRPRawToDigi_*_*',
    'keep TotemRPClusteredmDetSetVector_totemRPClusterProducer_*_*',
    'keep TotemRPRecHitedmDetSetVector_totemRPRecHitProducer_*_*',
    'keep TotemRPUVPatternedmDetSetVector_totemRPUVPatternFinder_*_*',
    'keep TotemRPLocalTrackedmDetSetVector_totemRPLocalTrackFitter_*_*'
  )
)

process.outpath = cms.EndPath(process.output)
