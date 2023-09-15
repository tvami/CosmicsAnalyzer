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
    # fileNames = cms.untracked.vstring("file:2ad63d9f-234b-4c68-8c18-49d485d42bc7.root")
#    fileNames = cms.untracked.vstring("file:00f26807-549d-45cf-844f-351e3a270f0e.root")
    # fileNames = cms.untracked.vstring("file:/home/users/dazhang/works/CosmicMuonSim/RAW-RECO/TRK-Run3Winter23Reco-00009-0to75-100000events.root")
    # fileNames = cms.untracked.vstring("file:/home/users/dazhang/works/CosmicMuonSim/CMSSW_12_6_5/src/TRK-Run3Winter23Reco-00009.root")
    fileNames = cms.untracked.vstring("file:/ceph/cms/store/user/lbrennan/EarthAsDM/Cosmics/crab_RAWtoReco-0to75Theta-4to3000GeV-126X_mcRun3_2022cosmics_realistic_deco_v1_v4/230811_215634/0000/3RR-0to75Theta-4to3000GeV_35.root")
)

#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(15000))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.MessageLogger.cerr.FwkReport.reportEvery = 100

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, options.GTAG, '')

process.muonPhiAnalyzer = cms.EDAnalyzer("EarthAsDMAnalyzer",
    isData = cms.untracked.int32(0)
    # isData = cms.untracked.int32(1)
)

process.TFileService = cms.Service("TFileService",
    # fileName = cms.string("ntuple_data.root")
    # fileName = cms.string("ntuple_MC_theta0to75_100kEvts.root")
    # fileName = cms.string("ntuple_MC_theta0to75_100Evts.root")
    fileName = cms.string("ntuple_MC_RR-0to75Theta-4to3000GeV_35_v4.root")
)

process.p = cms.Path(process.muonPhiAnalyzer)

