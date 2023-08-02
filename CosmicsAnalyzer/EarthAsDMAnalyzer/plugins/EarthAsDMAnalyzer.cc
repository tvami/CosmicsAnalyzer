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
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/DTRecHit/interface/DTRecHitCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtra.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtraMap.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"


#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
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
const int kMuonNMax = 100;
const int kSegmentNMax = 20;

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
  int verbose_;
  edm::EDGetTokenT< std::vector<reco::GenParticle> > genParticlesToken_;
  edm::EDGetTokenT<std::vector<reco::Muon>> muonToken_;
  edm::EDGetTokenT<reco::MuonTimeExtraMap> muonTimeToken_;
  edm::EDGetTokenT< CSCSegmentCollection > cscSegmentToken_;
  edm::EDGetTokenT< DTRecSegment4DCollection > dtSegmentToken_;
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
  int      muon_signOfGlobY_[kMuonNMax];
  
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
 verbose_(iConfig.getUntrackedParameter<int>("verbosityLevel")),
 genParticlesToken_(consumes< std::vector<reco::GenParticle> >( edm::InputTag("genParticles") )),
 muonToken_(consumes<std::vector<reco::Muon>>(iConfig.getParameter<edm::InputTag>("muonCollection"))),
 muonTimeToken_(consumes<reco::MuonTimeExtraMap>(iConfig.getParameter<edm::InputTag>("muonTimeCollection"))),
 cscSegmentToken_(consumes< CSCSegmentCollection >( edm::InputTag("cscSegments") )),
 dtSegmentToken_(consumes< DTRecSegment4DCollection >( edm::InputTag("dt4DSegments") )),
 muonDTGeomToken_(esConsumes())
{}

EarthAsDMAnalyzer::~EarthAsDMAnalyzer()  = default;

//
// member functions
//

// ------------ method called for each event  ------------
void EarthAsDMAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace std;
  using namespace edm;
  static constexpr const char* const MOD = "Analyzer";
  
  bool is_data_ = true;
  
  runNumber_ = iEvent.id().run();
  lsNumber_ = iEvent.id().luminosityBlock();
  eventNumber_ = iEvent.id().event();
  
  if (verbose_ > 1) LogPrint(MOD) << "Analyzing runNumber " << runNumber_ << " lsNumber " << lsNumber_ << " eventNumber " << eventNumber_;
  

  edm::Handle<reco::MuonCollection> muonCollectionHandle;
  iEvent.getByToken(muonToken_,muonCollectionHandle);
  
  edm::Handle<reco::MuonTimeExtraMap> tofMap;
  iEvent.getByToken(muonTimeToken_, tofMap);
  
  edm::Handle<CSCSegmentCollection> cscSegments;
  iEvent.getByToken(cscSegmentToken_, cscSegments);
  
  edm::Handle<DTRecSegment4DCollection> dtSegments;
  iEvent.getByToken(dtSegmentToken_, dtSegments);
  
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
    if (verbose_ > 2) LogPrint(MOD) << "  Analyzing track " << i ;
    reco::MuonRef muon  = reco::MuonRef( muonCollectionHandle, i );
    muon_pt_[muon_n_] = muon->pt();
    muon_p_[muon_n_] = muon->p();
    muon_eta_[muon_n_] = muon->eta();
    muon_phi_[muon_n_] = muon->phi();
    
    if (verbose_ > 2) LogPrint(MOD) << "  >> muon_pt_ " << muon->pt() << " muon_eta_ " << muon->eta() << " muon_phi_ " << muon->phi();
    
//    const reco::MuonTime time = muon->time();
//    const reco::MuonTime rpcTime = muon->rpcTime();

    // expectedNnumberOfMatchedStations
    int signOfGlobY = 0;
    if (muon->isMatchesValid()) {
      std::vector<reco::MuonChamberMatch>::const_iterator chamberMatch;
      for ( chamberMatch = muon->matches().begin(); chamberMatch != muon->matches().end(); ++chamberMatch) {
        const vector<reco::MuonSegmentMatch> matchedSegments = chamberMatch->segmentMatches;
        vector<reco::MuonSegmentMatch>::const_iterator segment;
        for (segment = matchedSegments.begin(); segment != matchedSegments.end(); ++segment) {
          edm::Ref<DTRecSegment4DCollection> dtSegment = segment->dtSegmentRef;
          if (dtSegment.isNull()) continue;
          LocalPoint segmentLocalPosition = dtSegment->localPosition();
          const GeomDet* dtDet = muonDTGeom->idToDet(dtSegment->geographicalId());
          GlobalPoint globalPoint = dtDet->toGlobal(segmentLocalPosition);
          if ( globalPoint.y() > 0) {
            signOfGlobY = +1;
          } else if ( globalPoint.y() < 0) {
            if (signOfGlobY > 0) cout << "!!!!!We have a problem, the track flipped CMS top to bottom." << endl;
            signOfGlobY = -1;
          }
        } // end loop on segments
      } // end loop on chamber matches
    } // end condition on muon having valid match
    muon_signOfGlobY_[muon_n_] = signOfGlobY;
    if (verbose_ > 3) LogPrint(MOD) << "  >> signOfGlobY " << signOfGlobY ;
    
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
      if (verbose_ > 4) LogPrint(MOD) << "  muon_comb_timeAtIpInOut_ " << combinedTimeExtra->timeAtIpInOut()  << " +/- " <<  combinedTimeExtra->timeAtIpOutInErr() << " muon_comb_timeAtIpOutIn_ " << combinedTimeExtra->timeAtIpOutIn() << " +/- " <<  combinedTimeExtra->timeAtIpInOutErr() << " muon_comb_invBeta_ " << combinedTimeExtra->inverseBeta() << " muon_comb_freeInvBeta_ " << combinedTimeExtra->freeInverseBeta();
      
//      if ((signOfGlobY*combinedTimeExtra->freeInverseBeta()) > 0) {
//        cout << "    >> This is signal " << endl;
//      } else if ((signOfGlobY*combinedTimeExtra->freeInverseBeta()) < 0) {
//        cout << "    >> This is background " << endl;
//      }
      if (fabs(combinedTimeExtra->timeAtIpInOutErr()) < fabs(combinedTimeExtra->timeAtIpOutInErr())) {
        if (signOfGlobY > 0) {cout << "    >> This is background " << endl;}
        else if (signOfGlobY < 0) {cout << "    >> This is signal " << endl;}
      } else {
        if (signOfGlobY < 0) {cout << "    >> This is background " << endl;}
        else if (signOfGlobY > 0) {cout << "    >> This is signal " << endl;}
      }
        
    }
    
    muon_n_++ ;
  }
    
//  for (unsigned int c=0; c<cscSegments->size(); c++) {
//    CSCSegmentRef segRef  = CSCSegmentRef( cscSegments, c );
//    const GeomDet* cscDet = cscGeom->idToDet(SegRef->geographicalId());
//    GlobalPoint point = cscDet->toGlobal(SegRef->localPosition());
//    cout << " X: " << point.x() << endl;
//  }
  
  // Fill the tree at the end of every event
  outputTree_->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void EarthAsDMAnalyzer::beginJob() {
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
  
  outputTree_ -> Branch ( "muon_n",                     &muon_n_);
  outputTree_ -> Branch ( "muon_pt",                    muon_pt_,                   "muon_pt[muon_n]/F");
  outputTree_ -> Branch ( "muon_p",                     muon_p_,                    "muon_p[muon_n]/F");
  outputTree_ -> Branch ( "muon_eta",                   muon_eta_,                  "muon_eta[muon_n]/F");
  outputTree_ -> Branch ( "muon_phi",                   muon_phi_,                  "muon_phi[muon_n]/F");
  outputTree_ -> Branch ( "muon_signOfGlobY",           muon_signOfGlobY_,          "muon_signOfGlobY[muon_n]/I");
  outputTree_ -> Branch ( "muon_comb_ndof",             muon_comb_ndof_,            "muon_comb_ndof[muon_n]/F");
  outputTree_ -> Branch ( "muon_comb_timeAtIpInOut",    muon_comb_timeAtIpInOut_,   "muon_comb_timeAtIpInOut[muon_n]/F");
  outputTree_ -> Branch ( "muon_comb_timeAtIpInOutErr", muon_comb_timeAtIpInOutErr_,"muon_comb_timeAtIpInOutErr[muon_n]/F");
  outputTree_ -> Branch ( "muon_comb_timeAtIpOutIn",    muon_comb_timeAtIpOutIn_,   "muon_comb_timeAtIpOutIn[muon_n]/F");
  outputTree_ -> Branch ( "muon_comb_timeAtIpOutInErr", muon_comb_timeAtIpOutInErr_,"muon_comb_timeAtIpOutInErr[muon_n]/F");
  outputTree_ -> Branch ( "muon_comb_invBeta",          muon_comb_invBeta_,         "muon_comb_invBeta[muon_n]/F");
  outputTree_ -> Branch ( "muon_comb_freeInvBeta",      muon_comb_freeInvBeta_,     "muon_comb_freeInvBeta[muon_n]/F");
  
}

// ------------ method called once each job just after ending the event loop  ------------
void EarthAsDMAnalyzer::endJob() { }

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void EarthAsDMAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setComment("Analyzer for cosmics searches");
  desc.addUntracked("verbosityLevel", 6)
  ->setComment("Higher the integer more verbose");
  desc.add("muonCollection", edm::InputTag("splitMuons"))
//  desc.add("muonCollection", edm::InputTag("lhcSTAMuons"))
//  desc.add("muonCollection", edm::InputTag("splitMuons"))
  ->setComment("Muon collection");
//  desc.add("muonTimeCollection", edm::InputTag("muons", "combined"))
//  desc.add("muonTimeCollection", edm::InputTag("lhcSTAMuons", "combined"))
//  desc.add("muonTimeCollection", edm::InputTag("muons1Leg", "combined"))
//  desc.add("muonTimeCollection", edm::InputTag("splitMuons", "combined"))
  desc.add("muonTimeCollection", edm::InputTag("splitMuons", "combined"))
  ->setComment("Input collection for combined muon timing information");
  descriptions.add("EarthAsDMAnalyzer",desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(EarthAsDMAnalyzer);

/*
// loop on the dt segments to get coordinates out
dtSeg_n_ = 0;
for (unsigned int d=0; d<dtSegments->size(); d++) {
  DTRecSegment4DRef segRef  = DTRecSegment4DRef( dtSegments, d );
  LocalPoint segmentLocalPosition = segRef->localPosition();
  LocalVector segmentLocalDirection = segRef->localDirection();
  LocalError segmentLocalPositionError = segRef->localPositionError();
  LocalError segmentLocalDirectionError = segRef->localDirectionError();
  const GeomDet* dtDet = muonDTGeom->idToDet(segRef->geographicalId());
  GlobalPoint globalPoint = dtDet->toGlobal(segmentLocalPosition);
  bool segmentFound = false;
  int segmentFoundMuonIndex = 0;
  
  for (unsigned int i = 0; i < muonCollectionHandle->size(); i++) {
    reco::MuonRef muon  = reco::MuonRef( muonCollectionHandle, i );
      // no point in looking if there are no matches
    if (!muon->isMatchesValid()) continue;
    
      // expectedNnumberOfMatchedStations
    for (std::vector<reco::MuonChamberMatch>::const_iterator chamberMatch = muon->matches().begin();
         chamberMatch != muon->matches().end();
         ++chamberMatch) {
      for (std::vector<reco::MuonSegmentMatch>::const_iterator segmentMatch = chamberMatch->segmentMatches.begin();
           segmentMatch != chamberMatch->segmentMatches.end();
           ++segmentMatch) {
        if (fabs(segmentMatch->x - segmentLocalPosition.x()) < 1E-6 &&
            fabs(segmentMatch->y - segmentLocalPosition.y()) < 1E-6 &&
            fabs(segmentMatch->dXdZ - segmentLocalDirection.x() / segmentLocalDirection.z()) < 1E-6 &&
            fabs(segmentMatch->dYdZ - segmentLocalDirection.y() / segmentLocalDirection.z()) < 1E-6 &&
            fabs(segmentMatch->xErr - sqrt(segmentLocalPositionError.xx())) < 1E-6 &&
            fabs(segmentMatch->yErr - sqrt(segmentLocalPositionError.yy())) < 1E-6 &&
            fabs(segmentMatch->dXdZErr - sqrt(segmentLocalDirectionError.xx())) < 1E-6 &&
            fabs(segmentMatch->dYdZErr - sqrt(segmentLocalDirectionError.yy())) < 1E-6) {
          segmentFound = true;
          segmentFoundMuonIndex = i;
          break;
        }
      }  // end loop on segments
      if (segmentFound) break;
    }  // end loop on chambers
    if (segmentFound)   break;
  }  // end loop on muon
  if (segmentFound) {
      //      cout << " Segment found for dtSeg_n_ index " << dtSeg_n_ << " and muon index " << segmentFoundMuonIndex << endl;
    dtSeg_muon_found_[dtSeg_n_][segmentFoundMuonIndex] = 1;
    dtSeg_muon_globX_[dtSeg_n_][segmentFoundMuonIndex] = globalPoint.x();
    dtSeg_muon_globY_[dtSeg_n_][segmentFoundMuonIndex] = globalPoint.y();
    dtSeg_muon_globZ_[dtSeg_n_][segmentFoundMuonIndex] = globalPoint.z();
  }
  
    //    if (segRef->hasPhi()) {
    //      const auto& dtPhiSegment = segRef->phiSegment();
    //      const auto& recHits = dtPhiSegment->specificRecHits();
    //      for (const auto& recHit : recHits) {
    //        auto digiTime = recHit.digiTime();
    //        cout << "1D-Phi hit digiTime: " << digiTime << endl;
    //      }
    //    }
  
    //    if (segRef->hasZed()) {
    //      const auto& dtZSegment = segRef->zSegment();
    //      const auto& recHits = dtZSegment->specificRecHits();
    //      for (const auto& recHit : recHits) {
    //        auto digiTime = recHit.digiTime();
    //        cout << "1D-Z hit digiTime: " << digiTime << endl;
    //      }
    //      }
  
    //    }
  dtSeg_n_++;
}  // dt segment
 
 int      dtSeg_n_;
 int      dtSeg_muon_found_[kSegmentNMax][kMuonNMax];
 float    dtSeg_muon_globX_[kSegmentNMax][kMuonNMax];
 float    dtSeg_muon_globY_[kSegmentNMax][kMuonNMax];
 float    dtSeg_muon_globZ_[kSegmentNMax][kMuonNMax];
 
 outputTree_ -> Branch ( "dtSeg_n",               &dtSeg_n_);
 outputTree_ -> Branch ( "dtSeg_muon_found",      dtSeg_muon_found_,     "dtSeg_muon_found[dtSeg_n][100]/I");
 outputTree_ -> Branch ( "dtSeg_muon_globX",      dtSeg_muon_globX_,     "dtSeg_muon_globX[dtSeg_n][100]/F");
 outputTree_ -> Branch ( "dtSeg_muon_globY",      dtSeg_muon_globY_,     "dtSeg_muon_globY[dtSeg_n][100]/F");
 outputTree_ -> Branch ( "dtSeg_muon_globZ",      dtSeg_muon_globZ_,     "dtSeg_muon_globZ[dtSeg_n][100]/F");
*/


