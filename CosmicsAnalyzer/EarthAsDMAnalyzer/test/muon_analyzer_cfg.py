# import FWCore.ParameterSet.Config as cms

# #process = cms.Process("MuonAnalyzer")
# process = cms.Process("EarthAsDMAnalyzer")

# process.source = cms.Source("PoolSource",
#     fileNames = cms.untracked.vstring("file:TRK-Run3Winter23Reco-00009.root")
# )

# process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))


# ## This is for MuonAnalyzer
# # process.muonPhiAnalyzer = cms.EDAnalyzer("MyAnalyzer",
# #     muonCollection = cms.InputTag("muons"),
# #     muonCollection2 = cms.InputTag("muons1Leg"),
# #     muonCollection3 = cms.InputTag("muonsBeamHaloEndCapsOnly"),
# #     muonCollection4 = cms.InputTag("muonsNoRPC"),
# #     muonCollection5 = cms.InputTag("muonsWitht0Correction"),
# #     muonCollection6 = cms.InputTag("splitMuons"),
# #     muonCollection7 = cms.InputTag("lhcSTAMuons"),
# # )

# ## This is for EarthAsDMAnalyzer
# process.muonPhiAnalyzer = cms.EDAnalyzer("MyAnalyzer",
#     genParticles = cms.InputTag("muons"),
#     muonCollection = cms.InputTag("muons1Leg"),
#     muonTimeCollection = cms.InputTag("muonsBeamHaloEndCapsOnly"),
#     cscSegments = cms.InputTag("muonsNoRPC"),
#     dt4DSegments = cms.InputTag("muonsWitht0Correction"),
# )




# process.TFileService = cms.Service("TFileService",
#     fileName = cms.string("MCVariables.root")
# )

# process.p = cms.Path(process.muonPhiAnalyzer)

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import sys, os

options = VarParsing('analysis')
options.register('GTAG', '106X_upgrade2018_realistic_v11BasedCandidateTmp_2022_08_09_01_32_34',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Global Tag"
)
options.parseArguments()

process = cms.Process("MyAnalyzer")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Configuration.StandardSequences.Services_cff')

process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring("file:/ceph/cms/store/user/lbrennan/EarthAsDM/Cosmics/crab_RAWtoReco-0to75Theta-4to3000GeV-126X_mcRun3_2022cosmics_realistic_deco_v1_v4/230811_215634/0000/3RR-0to75Theta-4to3000GeV_35.root")
)
#    fileNames = cms.untracked.vstring("file:00f26807-549d-45cf-844f-351e3a270f0e.root")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))


process.MessageLogger.cerr.FwkReport.reportEvery = 100

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, options.GTAG, '')

process.muonPhiAnalyzer = cms.EDAnalyzer("EarthAsDMAnalyzer",
    muonCollection = cms.InputTag("splitMuons"), #muons, lhcSTAMuons, splitMuons, or muons1Leg
#    splitMuonsCollection = cms.InputTag("splitMuons"),
#    muons1LegCollection = cms.InputTag("muons1Leg"),
#    lhcSTAMuonsCollection = cms.InputTag("lhcSTAMuons"),

#    muonCollection="muon1Leg"
# I can add other collections here later
    isData = cms.untracked.int32(0)
    # isData = cms.untracked.int32(1)
)

# process.options.SkipEvent = cms.untracked.vstring('runNumber 1 lsNumber 1 eventNumber 20858')
# process.options.SkipEvent = cms.untracked.vstring('NotFound')
# process.options.SkipEvent = cms.untracked.vstring('Analyzing track 0')

process.TFileService = cms.Service("TFileService",
    #fileName = cms.string("MCntuple-91to180-splitMuons-3Triggers-WithAdditionalVariables-test13.root")
    # fileName = cms.string("ntuple_data.root")
    # fileName = cms.string("ntuple_MC_theta0to75_100kEvts.root")
    # fileName = cms.string("ntuple_MC_theta0to75_100Evts.root")
    fileName = cms.string("ntuple_MC_RR-0to75Theta-4to3000GeV_35_v4.root")
)
#ImportantSampleProductions
#good100,000Samples1
process.p = cms.Path(process.muonPhiAnalyzer)






