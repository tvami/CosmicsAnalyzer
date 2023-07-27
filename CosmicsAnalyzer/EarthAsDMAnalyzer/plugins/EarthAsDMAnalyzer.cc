// -*- C++ -*-
//
// Package:    CosmicsAnalyzer/EarthAsDMAnalyzer
// Class:      EarthAsDMAnalyzer
//
/**\class EarthAsDMAnalyzer EarthAsDMAnalyzer.cc CosmicsAnalyzer/EarthAsDMAnalyzer/plugins/EarthAsDMAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Tamas Almos Vami
//         Created:  Fri, 21 Jul 2023 03:41:44 GMT
//
//

// system include files
#include <memory>
#include <vector>

// user include files
#include "AnalysisDataFormats/SUSYBSMObjects/interface/MuonSegment.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtra.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtraMap.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"

#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"

#include "TFile.h"
#include "TTree.h"
//
// class declaration
//

class EarthAsDMAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit EarthAsDMAnalyzer(const edm::ParameterSet&);
  ~EarthAsDMAnalyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<std::vector<reco::Muon>> muonToken_;
  edm::EDGetTokenT<reco::MuonTimeExtraMap> muonTimeToken_;
  edm::EDGetTokenT< CSCSegmentCollection > m_cscSegmentToken;
  edm::EDGetTokenT< DTRecSegment4DCollection > m_dtSegmentToken;
  edm::ESGetToken<CSCGeometry, MuonGeometryRecord> muonCSCGeomToken_;
  
  edm::ESGetToken<DTGeometry, MuonGeometryRecord> muonDTGeomToken_;
  const DTGeometry *muonDTGeom;
  
  TFile* outputFile_;
  TTree* outputTree_;
  unsigned int runNumber_;
  unsigned int lsNumber_;
  uint32_t eventNumber_;
  std::vector<float> muonPhis_;
  std::vector<float> muonPt_;
  std::vector<float> muonCombnDof_;
  std::vector<float> muonCombTimeAtIpInOut_;
  std::vector<float> muonCombTimeAtIpInOutErr_;
  std::vector<float> muonCombTimeAtIpOutIn_;
  std::vector<float> muonCombTimeAtIpOutInErr_;
  std::vector<float> muonCombinedInvBeta_;
  std::vector<float> muonCombinedFreeInvBeta_;

};

//
// constructors and destructor
//
EarthAsDMAnalyzer::EarthAsDMAnalyzer(const edm::ParameterSet& iConfig) :
 muonToken_(consumes<std::vector<reco::Muon>>(iConfig.getParameter<edm::InputTag>("muonCollection"))),
 muonTimeToken_(consumes<reco::MuonTimeExtraMap>(iConfig.getParameter<edm::InputTag>("muonTimeCollection"))),
 m_cscSegmentToken(consumes< CSCSegmentCollection >( edm::InputTag("cscSegments") )),
 m_dtSegmentToken(consumes< DTRecSegment4DCollection >( edm::InputTag("dt4DSegments") )),
 muonDTGeomToken_(esConsumes())
{}

EarthAsDMAnalyzer::~EarthAsDMAnalyzer()  = default;

//
// member functions
//

// ------------ method called for each event  ------------
void EarthAsDMAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  int verbose = 0;
  
  runNumber_ = iEvent.id().run();
  lsNumber_ = iEvent.id().luminosityBlock();
  eventNumber_ = iEvent.id().event();
  
  using namespace std;
  edm::Handle<reco::MuonCollection> muonCollectionHandle;
  iEvent.getByToken(muonToken_,muonCollectionHandle);
  
  edm::Handle<reco::MuonTimeExtraMap> tofMap;
  iEvent.getByToken(muonTimeToken_, tofMap);
  
  edm::Handle<CSCSegmentCollection> cscSegments;
  iEvent.getByToken(m_cscSegmentToken, cscSegments);
  
  edm::Handle<DTRecSegment4DCollection> dtSegments;
  iEvent.getByToken(m_dtSegmentToken, dtSegments);
  
  muonDTGeom = &iSetup.getData(muonDTGeomToken_);
  
//  cout << "This event has " << muonCollectionHandle->size() << " muons and " << dtSegments->size() << " segments" << endl;

//  for (const auto& muon : *muons) {
  for (unsigned int i = 0; i < muonCollectionHandle->size(); i++) {
    reco::MuonRef muon  = reco::MuonRef( muonCollectionHandle, i );
    float phi = muon->phi();
    float pt = muon->pt();
    muonPhis_.push_back(phi);
    muonPt_.push_back(pt);
    
//    const reco::MuonTime time = muon->time();
//    const reco::MuonTime rpcTime = muon->rpcTime();
    float nDof = 0.;
    float timeAtIpInOut = 9999.;
    float timeAtIpInOutErr = 9999.;
    float timeAtIpOutIn = 9999.;
    float timeAtIpOutInErr = 9999.;
    float combinedInvBeta = 9999.;
    float combinedFreeInvBeta = 9999.;
    if (tofMap.isValid()) {
      const reco::MuonTimeExtra* combinedTimeExtra = NULL;
      combinedTimeExtra = &tofMap->get(muon.key());

      
    nDof = combinedTimeExtra->nDof();
    timeAtIpInOut = combinedTimeExtra->timeAtIpInOut();
    timeAtIpInOutErr = combinedTimeExtra->timeAtIpInOutErr();
    timeAtIpOutIn = combinedTimeExtra->timeAtIpOutIn();
    timeAtIpOutInErr = combinedTimeExtra->timeAtIpOutInErr();
    combinedInvBeta = combinedTimeExtra->inverseBeta();
    combinedFreeInvBeta = combinedTimeExtra->freeInverseBeta();
      // Sign convention for muonCombinedFreeInvBeta_:
      //   positive - outward moving particle
      //   negative - inward moving particle
    }
    
    muonCombnDof_.push_back(nDof);
    muonCombTimeAtIpInOut_.push_back(timeAtIpInOut);
    muonCombTimeAtIpInOutErr_.push_back(timeAtIpInOutErr);
    muonCombTimeAtIpOutIn_.push_back(timeAtIpOutIn);
    muonCombTimeAtIpOutInErr_.push_back(timeAtIpOutInErr);
    muonCombinedInvBeta_.push_back(combinedInvBeta);
    muonCombinedFreeInvBeta_.push_back(combinedFreeInvBeta);
    
//    if (nDof > 1000) {
//      std::cout << "for i= " << i << " nDof is high! " << nDof << " phi " << phi << " and pt " << pt << std::endl;
//    std::cout << "combinedInvBeta: " << combinedInvBeta << " combinedFreeInvBeta: " << combinedFreeInvBeta << " timeAtIpInOut: " << timeAtIpInOut << " timeAtIpInOutErr:  "
//    << timeAtIpInOutErr << " timeAtIpOutIn: "  << timeAtIpOutIn << " timeAtIpOutInErr " << timeAtIpOutInErr << std::endl;
//    }
//    if (timeAtIpInOutErr > timeAtIpOutInErr)
//      return OutsideIn;
  }

  for (unsigned int c=0; c<cscSegments->size(); c++) {
    CSCSegmentRef segRef  = CSCSegmentRef( cscSegments, c );
    susybsm::MuonSegment muonSegment;
    muonSegment.setCSCSegmentRef(segRef);
//    const GeomDet* cscDet = cscGeom->idToDet(SegRef->geographicalId());
//    GlobalPoint point = cscDet->toGlobal(SegRef->localPosition());
//    cout << " X: " << point.x() << endl;
  }
  
  for (unsigned int d=0; d<dtSegments->size(); d++) {
    DTRecSegment4DRef SegRef  = DTRecSegment4DRef( dtSegments, d );
    susybsm::MuonSegment muonSegment;
    muonSegment.setDTSegmentRef(SegRef);
    
    const GeomDet* dtDet = muonDTGeom->idToDet(SegRef->geographicalId());
    GlobalPoint point = dtDet->toGlobal(SegRef->localPosition());
    if (verbose > 1) cout << " X: " << point.x() << endl;
  }

}

// ------------ method called once each job just before starting event loop  ------------
void EarthAsDMAnalyzer::beginJob() {
  outputFile_ = new TFile("ntuple.root", "RECREATE");
  outputTree_ = new TTree("tree", "Tree for this analysis");
  outputTree_->Branch("runNumber", &runNumber_);
  outputTree_->Branch("lsNumber", &lsNumber_);
  outputTree_->Branch("eventNumber", &eventNumber_);
  outputTree_->Branch("muonPhi", &muonPhis_);
  outputTree_->Branch("muonPt", &muonPt_);
  outputTree_->Branch("muonCombnDof", &muonCombnDof_);
  outputTree_->Branch("muonCombTimeAtIpInOut", &muonCombTimeAtIpInOut_);
  outputTree_->Branch("muonCombTimeAtIpInOutErr", &muonCombTimeAtIpInOutErr_);
  outputTree_->Branch("muonCombTimeAtIpOutIn", &muonCombTimeAtIpOutIn_);
  outputTree_->Branch("muonCombTimeAtIpOutInErr", &muonCombTimeAtIpOutInErr_);
  outputTree_->Branch("muonCombinedInvBeta", &muonCombinedInvBeta_);
  outputTree_->Branch("muonCombinedFreeInvBeta", &muonCombinedFreeInvBeta_);
  
}

// ------------ method called once each job just after ending the event loop  ------------
void EarthAsDMAnalyzer::endJob() {
  outputTree_->Fill();
  outputFile_->Write();
  outputFile_->Close();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void EarthAsDMAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setComment("Analyzer for cosmics searches");
  desc.add("muonCollection", edm::InputTag("muons"))
  ->setComment("Muon collection");
  desc.add("muonTimeCollection", edm::InputTag("muons1Leg", "combined"))
  ->setComment("Input collection for combined muon timing information");
  descriptions.add("EarthAsDMAnalyzer",desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(EarthAsDMAnalyzer);
