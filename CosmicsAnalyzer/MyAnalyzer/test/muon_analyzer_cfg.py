import FWCore.ParameterSet.Config as cms

process = cms.Process("MyAnalyzer")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring("file:2ad63d9f-234b-4c68-8c18-49d485d42bc7.root")
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.muonPhiAnalyzer = cms.EDAnalyzer("MyAnalyzer",
    muonCollection = cms.InputTag("muons")
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("muon_phi_ntuple.root")
)

process.p = cms.Path(process.muonPhiAnalyzer)

