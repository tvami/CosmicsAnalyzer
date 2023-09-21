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
#include <string>
#include <map>
#include <exception>
#include <unordered_map>

// ~~~~~~~~~ ROOT include files ~~~~~~~~~
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TTree.h"
#include "TObject.h"
#include "TVector3.h"
#include "TChain.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TProfile.h"
#include "TLorentzVector.h"
//#include "TCanvas.h"

// user include files
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/DTRecHit/interface/DTRecHitCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/DTRecHit/interface/DTRecHitCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
//muons
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtra.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtraMap.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"


#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"

//new
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/PFIsolation.h"
// #include "FWCore/Framework/interface/EDGetTokenT.h"




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


// #include "DataFormats/Common/interface/TriggerResults.h"
// #include "FWCore/Common/interface/TriggerNames.h"

#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"



#include "TFile.h"
#include "TTree.h"
//
// class declaration
//

const int kTrackNMax = 100;
const int kGenNMax = 100;
const int kMuonNMax = 100;
const int kSegmentNMax = 100;

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
  int isData_;
  edm::EDGetTokenT< std::vector<reco::GenParticle> > genParticlesToken_;
  edm::EDGetTokenT<std::vector<reco::Muon>> muonToken_;
  edm::EDGetTokenT<reco::MuonTimeExtraMap> muonTimeToken_;
  edm::EDGetTokenT< CSCSegmentCollection > cscSegmentToken_;
  edm::EDGetTokenT< DTRecSegment4DCollection > dtSegmentToken_;
  edm::ESGetToken<CSCGeometry, MuonGeometryRecord> muonCSCGeomToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> offlinePrimaryVerticesToken_;

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
  float    muon_energy_[kMuonNMax];

//New Test Variables
  bool    muon_IsLoose_[kMuonNMax];
  bool    muon_IsMedium_[kMuonNMax];
  bool    muon_IsTight_[kMuonNMax];
  bool    muon_isHighPtMuon_[kMuonNMax];
  bool    muon_isTrackerHighPtMuon_[kMuonNMax];

  unsigned int    muon_Type_;
  unsigned int    muon_Quality_;

  float    muon_d0_[kMuonNMax];
  float    muon_d0Err_[kMuonNMax];
  float    muon_charge_[kMuonNMax];
  float    muon_dZ_[kMuonNMax];

  float    muon_pileupIso_[kMuonNMax];
  float    muon_chargedIso_[kMuonNMax];
  float    muon_photonIso_[kMuonNMax];
  float    muon_neutralHadIso_[kMuonNMax];

  float    muon_validFractionTrackerHits_[kMuonNMax];
  float    muon_normChi2_[kMuonNMax];
  float    muon_chi2LocalPosition_[kMuonNMax];
  float    muon_kinkFinder_[kMuonNMax];
  float    muon_segmentCompatability_[kMuonNMax];
  float    muon_trkIso_[kMuonNMax];

  float    muon_tuneP_Pt_[kMuonNMax];
  float    muon_tuneP_PtErr_[kMuonNMax];
  float    muon_tuneP_Eta_[kMuonNMax];
  float    muon_tuneP_Phi_[kMuonNMax];
  float    muon_tuneP_MuonBestTrackType_[kMuonNMax];

//End Test Variables

  float    muon_comb_ndof_[kMuonNMax];
  float    muon_comb_timeAtIpInOut_[kMuonNMax];
  float    muon_comb_timeAtIpInOutErr_[kMuonNMax];
  float    muon_comb_timeAtIpOutIn_[kMuonNMax];
  float    muon_comb_timeAtIpOutInErr_[kMuonNMax];
  float    muon_comb_invBeta_[kMuonNMax];
  float    muon_comb_freeInvBeta_[kMuonNMax];
  int      muon_tofMap_found_[kMuonNMax];
  

  float    muon_dtSeg_x_[kMuonNMax][kSegmentNMax];
  float    muon_dtSeg_y_[kMuonNMax][kSegmentNMax];
  float    muon_dtSeg_z_[kMuonNMax][kSegmentNMax];
  
  int      muon_dtSeg_n_[kMuonNMax];
  int      dtSeg_n_;
  float    muon_dtSeg_t0timing_[kMuonNMax][kSegmentNMax];
  int      muon_dtSeg_found_[kMuonNMax][kSegmentNMax];
  float    muon_dtSeg_globX_[kMuonNMax][kSegmentNMax];
  float    muon_dtSeg_globY_[kMuonNMax][kSegmentNMax];
  float    muon_dtSeg_globZ_[kMuonNMax][kSegmentNMax];

  bool     trig_HLT_L1SingleMu18_v4_;
  bool     trig_HLT_L1SingleMu25_v3_;
  bool     trig_HLT_L1SingleMuCosmics_v2_;


};

//
// constructors and destructor
//
EarthAsDMAnalyzer::EarthAsDMAnalyzer(const edm::ParameterSet& iConfig) :
 verbose_(iConfig.getUntrackedParameter<int>("verbosityLevel")),
 isData_(iConfig.getUntrackedParameter<int>("isData")),
 genParticlesToken_(consumes< std::vector<reco::GenParticle> >( edm::InputTag("genParticles") )),
 muonToken_(consumes<std::vector<reco::Muon>>(iConfig.getParameter<edm::InputTag>("muonCollection"))),
 muonTimeToken_(consumes<reco::MuonTimeExtraMap>(iConfig.getParameter<edm::InputTag>("muonTimeCollection"))),
 cscSegmentToken_(consumes< CSCSegmentCollection >( edm::InputTag("cscSegments") )),
 dtSegmentToken_(consumes< DTRecSegment4DCollection >( edm::InputTag("dt4DSegments") )),
 triggerResultsToken_(consumes<edm::TriggerResults>(edm::InputTag(std::string("TriggerResults"),std::string(""),std::string("HLT")))),
 //offlinePrimaryVerticesToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("offlinePrimaryVertices"))),
 offlinePrimaryVerticesToken_(consumes<std::vector<reco::Vertex>>(edm::InputTag("offlinePrimaryVertices"))),
 //offlinePrimaryVerticesToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("offlinePrimaryVertices"))), 
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
    
  runNumber_ = iEvent.id().run();
  lsNumber_ = iEvent.id().luminosityBlock();
  eventNumber_ = iEvent.id().event();
  
  if (verbose_ > 1) LogPrint(MOD) << "Analyzing runNumber " << runNumber_ << " lsNumber " << lsNumber_ << " eventNumber " << eventNumber_;
  
  //------------------------------------------------------------------
  // Get trigger results for this event
  //------------------------------------------------------------------
  const edm::Handle<edm::TriggerResults> triggerH = iEvent.getHandle(triggerResultsToken_);
  // iEvent.getByToken(triggerResultsToken_,triggerH);
  const auto triggerNames = iEvent.triggerNames(*triggerH);

  edm::Handle<reco::MuonCollection> muonCollectionHandle;
  iEvent.getByToken(muonToken_,muonCollectionHandle);


  //vector<reco::Muon> muonColl = iEvent.get(muonToken_);
  
  edm::Handle<reco::MuonTimeExtraMap> tofMap;
  iEvent.getByToken(muonTimeToken_, tofMap);
  
  edm::Handle<CSCSegmentCollection> cscSegments;
  iEvent.getByToken(cscSegmentToken_, cscSegments);
  
  edm::Handle<DTRecSegment4DCollection> dtSegments;
  iEvent.getByToken(dtSegmentToken_, dtSegments);
  
  muonDTGeom = &iSetup.getData(muonDTGeomToken_);
  
  edm::Handle< std::vector<reco::GenParticle> > genColl;
  // vector<reco::Muon> muonColl = iEvent.get(muonToken_);

  //new

  edm::Handle<reco::VertexCollection> vertexColl;
  iEvent.getByToken(offlinePrimaryVerticesToken_, vertexColl);





  gen_n_ = 0;
  
  if (!isData_) {
    iEvent.getByToken(genParticlesToken_, genColl);
    for (unsigned int i = 0; i < genColl->size(); i++) {
      const reco::GenParticle* genCand = &(*genColl)[i];
      gen_pdg_[i] = genCand->pdgId();
      gen_pt_[i] = genCand->pt();
      gen_eta_[i] = genCand->eta();
      gen_phi_[i] = genCand->phi();
      gen_mass_[i] = genCand->mass();
      gen_isHardProcess_[i] = genCand->isHardProcess();
      gen_status_[i] = genCand->status();
      if (genCand->numberOfMothers() > 0) {
        gen_moth_pdg_[i] = genCand->mother()->pdgId();
      } else {
        gen_moth_pdg_[i] = 0;
      }
      gen_daughter_n_[i] = genCand->numberOfDaughters();
      if (genCand->numberOfDaughters() > 0) {
        gen_daughter_pdg_[i] = genCand->daughter(0)->pdgId();
      } else {
        gen_daughter_pdg_[i] = 0;
      }

      gen_n_++;
    }
  }
  if (gen_n_ > kGenNMax) cout << "!!!! please increase kGenNMax to " << gen_n_ << endl;
  

  //int highestSumPt2VertexIndex = -1;
  const reco::Vertex& highestSumPt2Vertex = vertexColl->front();


  muon_n_ = 0;
  for (unsigned int i = 0; i < muonCollectionHandle->size(); i++) {
    if (verbose_ > 2) LogPrint(MOD) << "\n  Analyzing track " << i ;
    reco::MuonRef muon  = reco::MuonRef( muonCollectionHandle, i );
    muon_pt_[muon_n_] = muon->pt();
    muon_p_[muon_n_] = muon->p();
    muon_eta_[muon_n_] = muon->eta();
    muon_phi_[muon_n_] = muon->phi();
    muon_energy_[muon_n_] = muon->energy();
    muon_charge_[muon_n_] = muon->charge();

    //New Variables
    //bool Variables
    bool isMedium = muon::isMediumMuon(*muon);
    muon_IsMedium_[muon_n_] = isMedium;

    bool isTight = muon::isTightMuon(*muon, highestSumPt2Vertex);
    muon_IsTight_[muon_n_] = isTight;

    bool isHigh = muon::isHighPtMuon(*muon, highestSumPt2Vertex);
    muon_isHighPtMuon_[muon_n_] = isHigh;

    bool isTrackerHighPtMuon = muon::isTrackerHighPtMuon(*muon, highestSumPt2Vertex);
    muon_isTrackerHighPtMuon_[muon_n_] = isTrackerHighPtMuon;


    //Positional Variables

    // Calculate the transverse impact parameter (d0) of the muon with respect to the highest sum pt vertex
    const reco::TrackRef bestTrack = muon->muonBestTrack();
    double dxy = -bestTrack->dxy(highestSumPt2Vertex.position());
    double dxyError = -bestTrack->dxyError();
    double dz = -bestTrack->dz(highestSumPt2Vertex.position());
    muon_d0_[muon_n_] = dxy;
    muon_d0Err_[muon_n_] = dxyError;
    muon_dZ_[muon_n_] = dz;

    muon_validFractionTrackerHits_[muon_n_] = (muon->innerTrack().isNonnull() ? muon->track()->validFraction() : -99.0);



    //only give 0s
    muon_pileupIso_[muon_n_] = muon->pfIsolationR04().sumPUPt;
    muon_chargedIso_[muon_n_] = muon->pfIsolationR04().sumChargedHadronPt;
    muon_photonIso_[muon_n_] = muon->pfIsolationR04().sumPhotonEt;
    muon_neutralHadIso_[muon_n_] = muon->pfIsolationR04().sumNeutralHadronEt;
    muon_trkIso_[muon_n_] = muon->isolationR03().sumPt;
    muon_chi2LocalPosition_[muon_n_] = muon->combinedQuality().chi2LocalPosition;
    muon_kinkFinder_[muon_n_] = muon->combinedQuality().trkKink;
    //only give 0s



    //not sure if working
    muon_Type_ = (muon->isMuon() + 2*muon->isGlobalMuon() + 4*muon->isTrackerMuon() + 8*muon->isStandAloneMuon()
        + 16*muon->isCaloMuon() + 32*muon->isPFMuon() + 64*muon->isRPCMuon());


    muon_tuneP_Pt_[muon_n_] = muon->tunePMuonBestTrack()->pt();
    muon_tuneP_PtErr_[muon_n_] = muon->tunePMuonBestTrack()->ptError();
    muon_tuneP_Eta_[muon_n_] = muon->tunePMuonBestTrack()->eta();
    muon_tuneP_Phi_[muon_n_] = muon->tunePMuonBestTrack()->phi();
    muon_tuneP_MuonBestTrackType_[muon_n_] = muon->tunePMuonBestTrackType();
    muon_segmentCompatability_[muon_n_] = (muon::segmentCompatibility(*muon));

    muon_Quality_ = (
        muon::isGoodMuon(*muon,muon::All)
        + pow(2,1)*muon::isGoodMuon(*muon,muon::AllGlobalMuons)
      + pow(2,2)*muon::isGoodMuon(*muon,muon::AllStandAloneMuons)
      + pow(2,3)*muon::isGoodMuon(*muon,muon::AllTrackerMuons)
      + pow(2,4)*muon::isGoodMuon(*muon,muon::TrackerMuonArbitrated)
      + pow(2,5)*muon::isGoodMuon(*muon,muon::AllArbitrated)
      + pow(2,6)*muon::isGoodMuon(*muon,muon::GlobalMuonPromptTight)
      + pow(2,7)*muon::isGoodMuon(*muon,muon::TMLastStationLoose)
      + pow(2,8)*muon::isGoodMuon(*muon,muon::TMLastStationTight)
      + pow(2,9)*muon::isGoodMuon(*muon,muon::TM2DCompatibilityLoose)
      + pow(2,10)*muon::isGoodMuon(*muon,muon::TM2DCompatibilityTight)
      + pow(2,11)*muon::isGoodMuon(*muon,muon::TMOneStationLoose)
      + pow(2,12)*muon::isGoodMuon(*muon,muon::TMOneStationTight)
      + pow(2,13)*muon::isGoodMuon(*muon,muon::TMLastStationOptimizedLowPtLoose)
      + pow(2,14)*muon::isGoodMuon(*muon,muon::TMLastStationOptimizedLowPtTight)
      + pow(2,15)*muon::isGoodMuon(*muon,muon::GMTkChiCompatibility)
      + pow(2,16)*muon::isGoodMuon(*muon,muon::GMStaChiCompatibility)
      + pow(2,17)*muon::isGoodMuon(*muon,muon::GMTkKinkTight)
      + pow(2,18)*muon::isGoodMuon(*muon,muon::TMLastStationAngLoose)
      + pow(2,19)*muon::isGoodMuon(*muon,muon::TMLastStationAngTight)
      + pow(2,20)*muon::isGoodMuon(*muon,muon::TMOneStationAngLoose)
      + pow(2,21)*muon::isGoodMuon(*muon,muon::TMOneStationAngTight)
      + pow(2,22)*muon::isGoodMuon(*muon,muon::TMLastStationOptimizedBarrelLowPtLoose)
      + pow(2,23)*muon::isGoodMuon(*muon,muon::TMLastStationOptimizedBarrelLowPtTight)
      + pow(2,24)*muon::isGoodMuon(*muon,muon::RPCMuLoose)      // //This is the soft muon ID
      + pow(2,25)*( muon::isGoodMuon(*muon,muon::TMOneStationTight)
        && muon->innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5
        && muon->innerTrack()->hitPattern().pixelLayersWithMeasurement() > 0
        && muon->innerTrack()->quality(reco::TrackBase::highPurity)
        && fabs(muon->innerTrack()->dxy(highestSumPt2Vertex.position())) < 0.3
        && fabs(muon->innerTrack()->dz(highestSumPt2Vertex.position())) < 20.
      ));



    muon_normChi2_[muon_n_] = (muon::isGoodMuon(*muon,muon::AllGlobalMuons) ? muon->globalTrack()->normalizedChi2() : -99.0);
    // muonQuality_ =
    if (verbose_) {
    std::cout << "Verbose output: muon_normChi2_[" << muon_n_ << "] = " << muon_normChi2_[muon_n_] << std::endl;
}







    
    if (verbose_ > 2) LogPrint(MOD) << "  >> muon_pt_ " << muon->pt() << " muon_eta_ " << muon->eta() << " muon_phi_ " << muon->phi();
    
//    const reco::MuonTime time = muon->time();
//    const reco::MuonTime rpcTime = muon->rpcTime();
    
    if (muon->isMatchesValid()) {
      // Loop on the chambers belonging to this muon
      int dtChamb_n_ = 0;
      std::vector<reco::MuonChamberMatch>::const_iterator chamberMatch;
      dtSeg_n_ = 0;
      
      for ( chamberMatch = muon->matches().begin(); chamberMatch != muon->matches().end(); ++chamberMatch) {
        if (verbose_ > 3) LogPrint(MOD)  << "    >> Chamber index " << dtChamb_n_;
        const vector<reco::MuonSegmentMatch> matchedSegments = chamberMatch->segmentMatches;
        vector<reco::MuonSegmentMatch>::const_iterator segment;
        // Now loop on the segments in the chamber
        for (segment = matchedSegments.begin(); segment != matchedSegments.end(); ++segment) {
          edm::Ref<DTRecSegment4DCollection> dtSegment = segment->dtSegmentRef;
          int found = 0;
          float t0timing = 9999;
          float dtGlobalPointX = 9999;
          float dtGlobalPointY = 9999;
          float dtGlobalPointZ = 9999;
          
          if (!dtSegment.isNull()) {
            found = 1;
            if (verbose_ > 3) LogPrint(MOD)  << "      >> DT segment index " << dtSeg_n_;
            LocalPoint segmentLocalPosition = dtSegment->localPosition();
            if (dtSegment->hasPhi()) {
              const auto& dtPhiSegment = dtSegment->phiSegment();
              t0timing = dtPhiSegment->t0();
              if (verbose_ > 4) LogPrint(MOD) << "        >> t0timing: " << t0timing;
            } else {
              if (verbose_ > 5) LogPrint(MOD) << "        >> This 4D segment does not have a phi segment: ";
              if (dtSegment->hasZed()) {
                if (verbose_ > 5) LogPrint(MOD) << "          >> But it has a zed segment: ";
              } else {
                if (verbose_ > 5) LogPrint(MOD) << "          >> Neither does it has a zed segment: ";
              }
            }
            const GeomDet* dtDet = muonDTGeom->idToDet(dtSegment->geographicalId());
            GlobalPoint globalPoint = dtDet->toGlobal(segmentLocalPosition);
            dtGlobalPointX = globalPoint.x();
            dtGlobalPointY = globalPoint.y();
            dtGlobalPointZ = globalPoint.z();
          }
          muon_dtSeg_found_[muon_n_][dtSeg_n_] = found;
          muon_dtSeg_t0timing_[muon_n_][dtSeg_n_] = t0timing;
          muon_dtSeg_globX_[muon_n_][dtSeg_n_] = dtGlobalPointX;
          muon_dtSeg_globY_[muon_n_][dtSeg_n_] = dtGlobalPointY;
          muon_dtSeg_globZ_[muon_n_][dtSeg_n_] = dtGlobalPointZ;

          dtSeg_n_++;
        } // end loop on segments
        dtChamb_n_++;
      } // end loop on chamber matches
      muon_dtSeg_n_[muon_n_] = dtSeg_n_;
      if (verbose_ > 3) LogPrint(MOD)  << "  >> This track had " << dtSeg_n_ << " segments";
    } // end condition on muon having valid match
  
    // if (tofMap.isValid()) {
    //   const reco::MuonTimeExtra* combinedTimeExtra = NULL;
    //   combinedTimeExtra = &tofMap->get(muon.key());
    // if (verbose_ > 2) LogPrint(MOD) << "  >> muon_pt_ " << muon->pt() << " muon_eta_ " << muon->eta() << " muon_phi_ " << muon->phi();
    
//    const reco::MuonTime time = muon->time();
//    const reco::MuonTime rpcTime = muon->rpcTime();
    
    // if (muon->isMatchesValid()) {
    //   // Loop on the chambers belonging to this muon
    //   int dtChamb_n_ = 0;
    //   std::vector<reco::MuonChamberMatch>::const_iterator chamberMatch;
    //   dtSeg_n_ = 0;
      
    //   for ( chamberMatch = muon->matches().begin(); chamberMatch != muon->matches().end(); ++chamberMatch) {
    //     if (verbose_ > 3) LogPrint(MOD)  << "    >> Chamber index " << dtChamb_n_;
    //     const vector<reco::MuonSegmentMatch> matchedSegments = chamberMatch->segmentMatches;
    //     vector<reco::MuonSegmentMatch>::const_iterator segment;
    //     // Now loop on the segments in the chamber
    //     for (segment = matchedSegments.begin(); segment != matchedSegments.end(); ++segment) {
    //       edm::Ref<DTRecSegment4DCollection> dtSegment = segment->dtSegmentRef;
    //       int found = 0;
    //       float t0timing = 9999;
    //       float dtGlobalPointX = 9999;
    //       float dtGlobalPointY = 9999;
    //       float dtGlobalPointZ = 9999;
          
    //       if (!dtSegment.isNull()) {
    //         found = 1;
    //         if (verbose_ > 3) LogPrint(MOD)  << "      >> DT segment index " << dtSeg_n_;
    //         LocalPoint segmentLocalPosition = dtSegment->localPosition();
    //         if (dtSegment->hasPhi()) {
    //           const auto& dtPhiSegment = dtSegment->phiSegment();
    //           t0timing = dtPhiSegment->t0();
    //           if (verbose_ > 4) LogPrint(MOD) << "        >> t0timing: " << t0timing;
    //         } else {
    //           if (verbose_ > 5) LogPrint(MOD) << "        >> This 4D segment does not have a phi segment: ";
    //           if (dtSegment->hasZed()) {
    //             if (verbose_ > 5) LogPrint(MOD) << "          >> But it has a zed segment: ";
    //           } else {
    //             if (verbose_ > 5) LogPrint(MOD) << "          >> Neither does it has a zed segment: ";
    //           }
    //         }
    //         const GeomDet* dtDet = muonDTGeom->idToDet(dtSegment->geographicalId());
    //         GlobalPoint globalPoint = dtDet->toGlobal(segmentLocalPosition);
    //         dtGlobalPointX = globalPoint.x();
    //         dtGlobalPointY = globalPoint.y();
    //         dtGlobalPointZ = globalPoint.z();
    //       }
    //       muon_dtSeg_found_[muon_n_][dtSeg_n_] = found;
    //       muon_dtSeg_t0timing_[muon_n_][dtSeg_n_] = t0timing;
    //       muon_dtSeg_globX_[muon_n_][dtSeg_n_] = dtGlobalPointX;
    //       muon_dtSeg_globY_[muon_n_][dtSeg_n_] = dtGlobalPointY;
    //       muon_dtSeg_globZ_[muon_n_][dtSeg_n_] = dtGlobalPointZ;

    //       dtSeg_n_++;
    //     } // end loop on segments
    //     dtChamb_n_++;
    //   // } // end loop on chamber matches
    //   // muon_dtSeg_n_[muon_n_] = dtSeg_n_;
    //   // if (verbose_ > 3) LogPrint(MOD)  << "  >> This track had " << dtSeg_n_ << " segments";
    // } // end condition on muon having valid match
  
    if (tofMap.isValid()) {
      const reco::MuonTimeExtra* combinedTimeExtra = NULL;
      combinedTimeExtra = &tofMap->get(muon.key());
      muon_tofMap_found_[muon_n_] = 1;

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
    }
    
    muon_n_++ ;
  
  
    
  }
    // trigInfo_ = 0;


  trig_HLT_L1SingleMu18_v4_ = false;
  trig_HLT_L1SingleMu25_v3_ = false;
  trig_HLT_L1SingleMuCosmics_v2_ = false;


  for (unsigned int i = 0; i < triggerH->size(); i++) {

    if (TString(triggerNames.triggerName(i)).Contains("HLT_L1SingleMu18_v") && triggerH->accept(i))
       trig_HLT_L1SingleMu18_v4_ = true;
        //cout << " HLT_Mu50 True? " << HLT_Mu50 << endl; 
    if (TString(triggerNames.triggerName(i)).Contains("HLT_L1SingleMu25_v") && triggerH->accept(i))
       trig_HLT_L1SingleMu25_v3_ = true;
    if (TString(triggerNames.triggerName(i)).Contains("HLT_L1SingleMuCosmics_v") && triggerH->accept(i))
       trig_HLT_L1SingleMuCosmics_v2_ = true;
  }

  









      
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
  outputTree_ -> Branch ( "muon_energy",                muon_energy_,               "muon_energy[muon_n]/F");

//new vairables
//floats
  outputTree_ -> Branch ( "muon_charge",                muon_charge_,               "muon_charge[muon_n]/F");


//bool
  outputTree_ -> Branch ( "muon_IsLoose",                muon_IsLoose_,               "muon_IsLoose[muon_n]/F");
  outputTree_ -> Branch ( "muon_IsMedium",                muon_IsMedium_,               "muon_IsMedium[muon_n]/F");
  outputTree_ -> Branch ( "muon_IsTight",                muon_IsTight_,               "muon_IsTight[muon_n]/F");


  outputTree_ -> Branch ( "muon_d0",                muon_d0_,               "muon_d0[muon_n]/F");
  outputTree_ -> Branch ( "muon_d0Err",                muon_d0Err_,               "muon_d0Err[muon_n]/F");
  outputTree_ -> Branch ( "muon_dZ",                muon_dZ_,               "muon_dZ[muon_n]/F");
  outputTree_ -> Branch ( "muon_Type",                &muon_Type_,               "muon_Type[muon_n]/I");

  outputTree_ -> Branch ( "muon_pileupIso",                muon_pileupIso_,               "muon_pileupIso[muon_n]/F");
  outputTree_ -> Branch ( "muon_chargedIso",                muon_chargedIso_,               "muon_chargedIso[muon_n]/F");
  outputTree_ -> Branch ( "muon_photonIso",                muon_photonIso_,               "muon_photonIso[muon_n]/F");
  outputTree_ -> Branch ( "muon_neutralHadIso",                muon_neutralHadIso_,               "muon_neutralHadIso[muon_n]/F");
  outputTree_ -> Branch ( "muon_validFractionTrackerHits",                muon_validFractionTrackerHits_,               "muon_validFractionTrackerHits[muon_n]/F");

  outputTree_ -> Branch ( "muon_tuneP_Pt",                muon_tuneP_Pt_,               "muon_tuneP_Pt[muon_n]/F");
  outputTree_ -> Branch ( "muon_tuneP_PtErr",                muon_tuneP_PtErr_,               "muon_tuneP_PtErr[muon_n]/F");
  outputTree_ -> Branch ( "muon_tuneP_Eta",                muon_tuneP_Eta_,               "muon_tuneP_Eta[muon_n]/F");
  outputTree_ -> Branch ( "muon_tuneP_Phi",                muon_tuneP_Phi_,               "muon_tuneP_Phi[muon_n]/F");
  outputTree_ -> Branch ( "muon_tuneP_MuonBestTrackType",                muon_tuneP_MuonBestTrackType_,               "muon_tuneP_MuonBestTrackType[muon_n]/F");

  outputTree_ -> Branch ( "muon_trkIso",                muon_trkIso_,               "muon_trkIso[muon_n]/F");
  outputTree_ -> Branch ( "muon_normChi2",                muon_normChi2_,               "muon_normChi2[muon_n]/F");
  outputTree_ -> Branch ( "muon_chi2LocalPosition",                muon_chi2LocalPosition_,               "muon_chi2LocalPosition[muon_n]/F");
  outputTree_ -> Branch ( "muon_kinkFinder",                muon_kinkFinder_,               "muon_kinkFinder[muon_n]/F");
  outputTree_ -> Branch ( "muon_segmentCompatability",                muon_segmentCompatability_,               "muon_segmentCompatability[muon_n]/F");

  outputTree_ -> Branch ( "muon_isTrackerHighPtMuon",                muon_isTrackerHighPtMuon_,               "muon_isTrackerHighPtMuon[muon_n]/O");
  outputTree_ -> Branch ( "muon_isHighPtMuon",                muon_isHighPtMuon_,               "muon_isHighPtMuon[muon_n]/O");
  outputTree_ -> Branch ( "muon_Quality",                &muon_Quality_,               "muon_Quality[muon_n]/I");

//End new Variables

  outputTree_ -> Branch ( "muon_dtSeg_n",          muon_dtSeg_n_,         "muon_dtSeg_n[muon_n]/I");
  outputTree_ -> Branch ( "muon_dtSeg_t0timing",   muon_dtSeg_t0timing_,  "muon_dtSeg_t0timing[muon_n][100]/F");
  outputTree_ -> Branch ( "muon_dtSeg_found",      muon_dtSeg_found_,     "muon_dtSeg_found[muon_n][100]/I");
  outputTree_ -> Branch ( "muon_dtSeg_globX",      muon_dtSeg_globX_,     "muon_dtSeg_globX[muon_n][100]/F");
  outputTree_ -> Branch ( "muon_dtSeg_globY",      muon_dtSeg_globY_,     "muon_dtSeg_globY[muon_n][100]/F");
  outputTree_ -> Branch ( "muon_dtSeg_globZ",      muon_dtSeg_globZ_,     "muon_dtSeg_globZ[muon_n][100]/F");
  
  outputTree_ -> Branch ( "muon_comb_ndof",             muon_comb_ndof_,            "muon_comb_ndof[muon_n]/F");
  outputTree_ -> Branch ( "muon_comb_timeAtIpInOut",    muon_comb_timeAtIpInOut_,   "muon_comb_timeAtIpInOut[muon_n]/F");
  outputTree_ -> Branch ( "muon_comb_timeAtIpInOutErr", muon_comb_timeAtIpInOutErr_,"muon_comb_timeAtIpInOutErr[muon_n]/F");
  outputTree_ -> Branch ( "muon_comb_timeAtIpOutIn",    muon_comb_timeAtIpOutIn_,   "muon_comb_timeAtIpOutIn[muon_n]/F");
  outputTree_ -> Branch ( "muon_comb_timeAtIpOutInErr", muon_comb_timeAtIpOutInErr_,"muon_comb_timeAtIpOutInErr[muon_n]/F");
  outputTree_ -> Branch ( "muon_comb_invBeta",          muon_comb_invBeta_,         "muon_comb_invBeta[muon_n]/F");
  outputTree_ -> Branch ( "muon_comb_freeInvBeta",      muon_comb_freeInvBeta_,     "muon_comb_freeInvBeta[muon_n]/F");
  
  outputTree_ -> Branch ( "muon_dtSeg_x",      muon_dtSeg_x_,     "muon_dtSeg_x[muon_n][20]/F");
  outputTree_ -> Branch ( "muon_dtSeg_y",      muon_dtSeg_y_,     "muon_dtSeg_y[muon_n][20]/F");
  outputTree_ -> Branch ( "muon_dtSeg_z",      muon_dtSeg_z_,     "muon_dtSeg_z[muon_n][20]/F");
  outputTree_ -> Branch ( "muon_tofMap_found",          muon_tofMap_found_,         "muon_tofMap_found[muon_n]/F");

  outputTree_ -> Branch ( "trig_HLT_L1SingleMu18_v4",     &trig_HLT_L1SingleMu18_v4_, "trig_HLT_L1SingleMu18_v4/O") ;
  outputTree_ -> Branch ( "trig_HLT_L1SingleMu25_v3",     &trig_HLT_L1SingleMu25_v3_, "trig_HLT_L1SingleMu25_v3/O") ;
  outputTree_ -> Branch ( "trig_HLT_L1SingleMuCosmics_v2",     &trig_HLT_L1SingleMuCosmics_v2_, "trig_HLT_L1SingleMuCosmics_v2/O") ;




}

// ------------ method called once each job just after ending the event loop  ------------
void EarthAsDMAnalyzer::endJob() { }

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void EarthAsDMAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setComment("Analyzer for cosmics searches");
  desc.addUntracked("verbosityLevel", 6)
  ->setComment("Higher the integer more verbose");
  desc.addUntracked("isData", 0)
  ->setComment("0 means MC, 1 means data");
//  desc.add("muonCollection", edm::InputTag("splitMuons")) //muons1Leg
//  desc.add("muonCollection", edm::InputTag("lhcSTAMuons"))
  desc.add("muonCollection", edm::InputTag("splitMuons"))
  ->setComment("Muon collection");
  desc.add("muonTimeCollection", edm::InputTag("splitMuons", "dt"))
  ->setComment("Input collection for combined muon timing information");
  // desc.add("TriggerResults", edm::InputTag("TriggerResults","","HLT"))
  //   ->setComment("A");
  desc.addUntracked("Trigger_Mu", std::vector<std::string>{"HLT_L1SingleMu18_v", "HLT_L1SingleMu25_v"})
//  desc.addUntracked("Trigger_Mu", std::vector<std::string>{"HLT_Mu50_v","HLT_OldMu100_v","HLT_TkMu100_v"})
  ->setComment("Add the list of muon triggers");



//   desc.addUntracked("Trigger_Mu", std::vector<std::string>{"HLT_Mu50_v", "HLT_IsoMu24_v"})
// //  desc.addUntracked("Trigger_Mu", std::vector<std::string>{"HLT_Mu50_v","HLT_OldMu100_v","HLT_TkMu100_v"})
//   ->setComment("Add the list of muon triggers");

  descriptions.add("EarthAsDMAnalyzer",desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(EarthAsDMAnalyzer);

  //  for (unsigned int c=0; c<cscSegments->size(); c++) {
  //    CSCSegmentRef segRef  = CSCSegmentRef( cscSegments, c );
  //    const GeomDet* cscDet = cscGeom->idToDet(SegRef->geographicalId());
  //    GlobalPoint point = cscDet->toGlobal(SegRef->localPosition());
  //    cout << " X: " << point.x() << endl;
  //  }

//  // loop on the dt segments to get coordinates out
//dtSeg_n_ = 0;
//for (unsigned int d=0; d<dtSegments->size(); d++) {
//  DTRecSegment4DRef segRef  = DTRecSegment4DRef( dtSegments, d );
//  LocalPoint segmentLocalPosition = segRef->localPosition();
//  LocalVector segmentLocalDirection = segRef->localDirection();
//  LocalError segmentLocalPositionError = segRef->localPositionError();
//  LocalError segmentLocalDirectionError = segRef->localDirectionError();
//  const GeomDet* dtDet = muonDTGeom->idToDet(segRef->geographicalId());
//  GlobalPoint globalPoint = dtDet->toGlobal(segmentLocalPosition);
//  bool segmentFound = false;
//  int segmentFoundMuonIndex = 0;
//
//  for (unsigned int i = 0; i < muonCollectionHandle->size(); i++) {
//    reco::MuonRef muon  = reco::MuonRef( muonCollectionHandle, i );
//      //      no point in looking if there are no matches
//    if (!muon->isMatchesValid()) continue;
//
//      //      expectedNnumberOfMatchedStations
//    for (std::vector<reco::MuonChamberMatch>::const_iterator chamberMatch = muon->matches().begin();
//         chamberMatch != muon->matches().end();
//         ++chamberMatch) {
//      for (std::vector<reco::MuonSegmentMatch>::const_iterator segmentMatch = chamberMatch->segmentMatches.begin();
//           segmentMatch != chamberMatch->segmentMatches.end();
//           ++segmentMatch) {
//        if (fabs(segmentMatch->x - segmentLocalPosition.x()) < 1E-6 &&
//            fabs(segmentMatch->y - segmentLocalPosition.y()) < 1E-6 &&
//            fabs(segmentMatch->dXdZ - segmentLocalDirection.x() / segmentLocalDirection.z()) < 1E-6 &&
//            fabs(segmentMatch->dYdZ - segmentLocalDirection.y() / segmentLocalDirection.z()) < 1E-6 &&
//            fabs(segmentMatch->xErr - sqrt(segmentLocalPositionError.xx())) < 1E-6 &&
//            fabs(segmentMatch->yErr - sqrt(segmentLocalPositionError.yy())) < 1E-6 &&
//            fabs(segmentMatch->dXdZErr - sqrt(segmentLocalDirectionError.xx())) < 1E-6 &&
//            fabs(segmentMatch->dYdZErr - sqrt(segmentLocalDirectionError.yy())) < 1E-6) {
//          segmentFound = true;
//          segmentFoundMuonIndex = i;
//          break;
//        }
//      }  // end loop on segments
//      if (segmentFound) break;
//    } //  end loop on chambers
//    if (segmentFound)   break;
//  }  // end loop on muon
//  if (segmentFound) {
//    cout << " Segment found for dtSeg_n_ index " << dtSeg_n_ << " and muon index " << segmentFoundMuonIndex << " with Y =" << globalPoint.y() << endl;
//
//    if (segRef->hasPhi()) {
//      const auto& dtPhiSegment = segRef->phiSegment();
//      cout << "Timing: " << dtPhiSegment->t0() << endl;
//      const auto& recHits = dtPhiSegment->specificRecHits();
//      for (const auto& recHit : recHits) {
//        auto digiTime = recHit.digiTime();
//        cout << "1D-Phi hit digiTime: " << digiTime << endl;
//      }
//    } // end of segm has Phi
//
//    if (segRef->hasZed()) {
//      const auto& dtZSegment = segRef->zSegment();
//      cout << "Timing: " << dtZSegment->t0() << endl;
//      const auto& recHits = dtZSegment->specificRecHits();
//      for (const auto& recHit : recHits) {
//        auto digiTime = recHit.digiTime();
//        cout << "1D-Z hit digiTime: " << digiTime << endl;
//      }
//    } // end of segm has Zed
//  } // end of segm found
//}  // end loop on dt segment
