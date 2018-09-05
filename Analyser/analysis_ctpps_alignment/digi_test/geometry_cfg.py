import FWCore.ParameterSet.Config as cms
import fnmatch

process = cms.Process("GeomInfo")

# minimum of logs
process.MessageLogger = cms.Service("MessageLogger",
    statistics = cms.untracked.vstring(),
    destinations = cms.untracked.vstring('cout'),
    cout = cms.untracked.PSet(
        threshold = cms.untracked.string('WARNING')
    )
)

# raw data source
process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

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

process.geometryInfo = cms.EDAnalyzer("GeometryInfoModule")

process.p = cms.Path(
    process.geometryInfo
)
