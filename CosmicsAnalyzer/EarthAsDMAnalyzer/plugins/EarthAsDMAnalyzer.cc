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
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"

#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"

#include "TFile.h"
#include "TTree.h"
//
// class declaration
//

const int kTrackNMax = 10000;
const int kGenNMax = 10000;
const int kMuonNMax = 10000;

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
  edm::EDGetTokenT< std::vector<reco::GenParticle> > genParticlesToken_;
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

  
  int      gen_n_;
  int      gen_pdg_[kGenNMax];
  float    gen_pt_[kGenNMax];
  float    gen_eta_[kGenNMax];
  float    gen_phi_[kGenNMax];
  float    gen_mass_[kGenNMax];
  bool     gen_isHardProcess_[kGenNMax];
  int      gen_status_[kGenNMax];
  int      gen_moth_pdg_[kGenNMax];
  int      gen_daughter_n_[kGenNMax];
  int      gen_daughter_pdg_[kGenNMax];
  
  int      muon_n_;
  float    muon_pt_[kMuonNMax];
  float    muon_p_[kMuonNMax];
  float    muon_eta_[kMuonNMax];
  float    muon_phi_[kMuonNMax];
  float    muon_comb_ndof_[kMuonNMax];
  float    muon_comb_timeAtIpInOut_[kMuonNMax];
  float    muon_comb_timeAtIpInOutErr_[kMuonNMax];
  float    muon_comb_timeAtIpOutIn_[kMuonNMax];
  float    muon_comb_timeAtIpOutInErr_[kMuonNMax];
  float    muon_comb_invBeta_[kMuonNMax];
  float    muon_comb_freeInvBeta_[kMuonNMax];

};

//
// constructors and destructor
//
EarthAsDMAnalyzer::EarthAsDMAnalyzer(const edm::ParameterSet& iConfig) :
 genParticlesToken_(consumes< std::vector<reco::GenParticle> >( edm::InputTag("genParticles") )),
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
  bool is_data_ = true;
  
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
  
  edm::Handle< std::vector<reco::GenParticle> > genColl;
  gen_n_ = 0;
  if (!is_data_) {
    iEvent.getByToken(genParticlesToken_, genColl);
    for (unsigned int i = 0; i < genColl->size(); i++) {
//      const reco::GenParticle* genCand = &(*genColl)[i];
      gen_n_++;
    }
  }
  
  muon_n_ = 0;
  for (unsigned int i = 0; i < muonCollectionHandle->size(); i++) {
    reco::MuonRef muon  = reco::MuonRef( muonCollectionHandle, i );
    muon_pt_[muon_n_] = muon->pt();
    muon_p_[muon_n_] = muon->p();
    muon_eta_[muon_n_] = muon->eta();
    muon_phi_[muon_n_] = muon->phi();
    
//    const reco::MuonTime time = muon->time();
//    const reco::MuonTime rpcTime = muon->rpcTime();

    if (tofMap.isValid()) {
      const reco::MuonTimeExtra* combinedTimeExtra = NULL;
      combinedTimeExtra = &tofMap->get(muon.key());

      muon_comb_ndof_[muon_n_] = combinedTimeExtra->nDof();
      muon_comb_timeAtIpInOut_[muon_n_] = combinedTimeExtra->timeAtIpInOut();
      muon_comb_timeAtIpInOutErr_[muon_n_] = combinedTimeExtra->timeAtIpInOutErr();
      muon_comb_timeAtIpOutIn_[muon_n_] = combinedTimeExtra->timeAtIpOutIn();
      muon_comb_timeAtIpOutInErr_[muon_n_] = combinedTimeExtra->timeAtIpOutInErr();
      muon_comb_invBeta_[muon_n_] = combinedTimeExtra->inverseBeta();
      muon_comb_freeInvBeta_[muon_n_] = combinedTimeExtra->freeInverseBeta();
      // Sign convention for muonCombinedFreeInvBeta_:
      //   positive - outward moving particle
      //   negative - inward moving particle
    }
    
//    if (nDof > 1000) {
//      std::cout << "for i= " << i << " nDof is high! " << nDof << " phi " << phi << " and pt " << pt << std::endl;
//    std::cout << "combinedInvBeta: " << combinedInvBeta << " combinedFreeInvBeta: " << combinedFreeInvBeta << " timeAtIpInOut: " << timeAtIpInOut << " timeAtIpInOutErr:  "
//    << timeAtIpInOutErr << " timeAtIpOutIn: "  << timeAtIpOutIn << " timeAtIpOutInErr " << timeAtIpOutInErr << std::endl;
//    }
//    if (timeAtIpInOutErr > timeAtIpOutInErr)
//      return OutsideIn;
    muon_n_++ ;
  }

  for (unsigned int c=0; c<cscSegments->size(); c++) {
    CSCSegmentRef segRef  = CSCSegmentRef( cscSegments, c );
    susybsm::MuonSegment muonSegment;
    muonSegment.setCSCSegmentRef(segRef);
//    const GeomDet* cscDet = cscGeom->idToDet(SegRef->geographicalId());
//    GlobalPoint point = cscDet->toGlobal(SegRef->localPosition());
//    cout << " X: " << point.x() << endl;
  }
  
  //  cout << "This event has " << muonCollectionHandle->size() << " muons and " << dtSegments->size() << " segments" << endl;
  for (unsigned int d=0; d<dtSegments->size(); d++) {
    DTRecSegment4DRef SegRef  = DTRecSegment4DRef( dtSegments, d );
    susybsm::MuonSegment muonSegment;
    muonSegment.setDTSegmentRef(SegRef);
    
    const GeomDet* dtDet = muonDTGeom->idToDet(SegRef->geographicalId());
    GlobalPoint point = dtDet->toGlobal(SegRef->localPosition());
    if (verbose > 1) cout << " X: " << point.x() << endl;
  }
  outputTree_->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void EarthAsDMAnalyzer::beginJob() {
//  outputFile_ = new TFile("ntuple.root", "RECREATE");
//  outputTree_ = new TTree("tree", "Tree for this analysis");
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  outputTree_ = fs->make<TTree>("tree", "tree");
  outputTree_->Branch("run",    &runNumber_);
  outputTree_->Branch("ls",     &lsNumber_);
  outputTree_->Branch("event",  &eventNumber_);
  
  outputTree_ -> Branch ( "gen_n",            &gen_n_) ;
  outputTree_ -> Branch ( "gen_pdg",          gen_pdg_,          "gen_pdg[gen_n]/I");
  outputTree_ -> Branch ( "gen_pt",           gen_pt_,           "gen_pt[gen_n]/F");
  outputTree_ -> Branch ( "gen_eta",          gen_eta_,          "gen_eta[gen_n]/F");
  outputTree_ -> Branch ( "gen_phi",          gen_phi_,          "gen_phi[gen_n]/F");
  outputTree_ -> Branch ( "gen_mass",         gen_mass_,         "gen_mass[gen_n]/F");
  outputTree_ -> Branch ( "gen_isHardProcess",gen_isHardProcess_,"gen_isHardProcess[gen_n]/O");
  outputTree_ -> Branch ( "gen_status",       gen_status_,       "gen_status[gen_n]/I");
  outputTree_ -> Branch ( "gen_moth_pdg",     gen_moth_pdg_,     "gen_moth_pdg[gen_n]/I");
  outputTree_ -> Branch ( "gen_daughter_n",   gen_daughter_n_,   "gen_daughter_n[gen_n]/I");
  outputTree_ -> Branch ( "gen_daughter_pdg", gen_daughter_pdg_, "gen_daughter_pdg[gen_n]/I");
  
  outputTree_ -> Branch ( "muon_n",                     &muon_n_) ;
  outputTree_ -> Branch ( "muon_pt",                    muon_pt_,                   "muon_pt[muon_n]/F");
  outputTree_ -> Branch ( "muon_p",                     muon_p_,                    "muon_p[muon_n]/F");
  outputTree_ -> Branch ( "muon_eta",                   muon_eta_,                  "muon_eta[muon_n]/F");
  outputTree_ -> Branch ( "muon_phi",                   muon_phi_,                  "muon_phi[muon_n]/F");
  outputTree_ -> Branch ( "muon_comb_ndof",             muon_comb_ndof_,            "muon_comb_ndof[muon_n]/F");
  outputTree_ -> Branch ( "muon_comb_timeAtIpInOut",    muon_comb_timeAtIpInOut_,   "muon_comb_timeAtIpInOut[muon_n]/F");
  outputTree_ -> Branch ( "muon_comb_timeAtIpInOutErr", muon_comb_timeAtIpInOutErr_,"muon_comb_timeAtIpInOutErr[muon_n]/F");
  outputTree_ -> Branch ( "muon_comb_timeAtIpOutIn",    muon_comb_timeAtIpOutIn_,   "muon_comb_timeAtIpOutIn[muon_n]/F");
  outputTree_ -> Branch ( "muon_comb_timeAtIpOutInErr", muon_comb_timeAtIpOutInErr_,"muon_comb_timeAtIpOutInErr[muon_n]/F");
  outputTree_ -> Branch ( "muon_comb_invBeta",          muon_comb_invBeta_,         "muon_comb_invBeta[muon_n]/F");
  outputTree_ -> Branch ( "muon_comb_freeInvBeta",      muon_comb_freeInvBeta_,     "muon_comb_freeInvBeta[muon_n]/F");
  
}

// ------------ method called once each job just after ending the event loop  ------------
void EarthAsDMAnalyzer::endJob() {
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
