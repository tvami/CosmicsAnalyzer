// -*- C++ -*-
//
// Package:    CosmicsAnalyzer/EarthAsDMAnalyzer
// Class:      EarthAsDMAnalyzer
//
/**\class EarthAsDMAnalyzer EarthAsDMAnalyzer.cc CosmicsAnalyzer/EarthAsDMAnalyzer/plugins/EarthAsDMAnalyzer.cc

 Description: Ntuplizer code ran on RAW-RECO files to extract muon and track information

 Implementation:
     Oct 8: moved everything to vectors instead of arrays
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
#include <iostream>
#include <cmath>

// ~~~~~~~~~ ROOT include files ~~~~~~~~~
#include "TDirectory.h"
#include "TFile.h"
#include "TTree.h"
#include "TObject.h"
#include "TVector3.h"
#include "TChain.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TLorentzVector.h"


// ~~~~~~~~~ ROOT include files ~~~~~~~~~

// ~~~~~~~~~ CMSSW include files ~~~~~~~~~
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/DTRecHit/interface/DTRecHitCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/DTLayerId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtra.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtraMap.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "DataFormats/PatCandidates/interface/PFIsolation.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "CondFormats/DTObjects/interface/DTT0.h"
#include "CondFormats/DataRecord/interface/DTT0Rcd.h"
#include "CondFormats/DTObjects/interface/DTTtrig.h"
#include "CondFormats/DataRecord/interface/DTTtrigRcd.h"

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
  int verbose_;
  int isData_;
  int hasSim_;
  edm::EDGetTokenT< std::vector<reco::GenParticle> > genParticlesToken_;
  edm::EDGetTokenT<std::vector<reco::Muon>> muonToken_;
  edm::EDGetTokenT<reco::MuonTimeExtraMap> muonTimeToken_;
  edm::EDGetTokenT<edm::PSimHitContainer> PSimHitContainerToken_;
  edm::EDGetTokenT< CSCSegmentCollection > cscSegmentToken_;
  edm::EDGetTokenT< DTRecSegment4DCollection > dtSegmentToken_;
  edm::ESGetToken<CSCGeometry, MuonGeometryRecord> muonCSCGeomToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> offlinePrimaryVerticesToken_;
  const DTGeometry *muonDTGeom;
  const DTT0* t0Map;
  const DTTtrig* tTrigMap;
  edm::ESGetToken<DTGeometry, MuonGeometryRecord> muonDTGeomToken_;
  edm::ESGetToken<DTT0, DTT0Rcd> t0Token_;
  edm::ESGetToken<DTTtrig, DTTtrigRcd> ttrigToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
  edm::EDGetTokenT<std::vector<reco::Track>> tracksToken_;
  
  TFile* outputFile_;
  TTree* outputTree_;

  unsigned int runNumber_;
  unsigned int lsNumber_;
  uint32_t eventNumber_;

  // trigger variable (event level)
  bool     HLT_L1SingleMu3_;
  bool     HLT_L1SingleMu5_;
  bool     HLT_L1SingleMu7_;
  bool     HLT_L1SingleMu18_;
  bool     HLT_L2Mu23NoVtx_2Cha_CosmicSeed_;
  bool     HLT_L1SingleMuOpen_;
  bool     HLT_L1SingleMuOpen_DT_;
  bool     HLT_L1SingleMuCosmics_;
  
  // Variables that are at event level are numbers (e.g. uint)
  unsigned int gen_n_;
  // The size of these vectors will be gen_n_
  std::vector<int>   gen_pdg_;
  std::vector<float> gen_energy_;
  std::vector<float> gen_pt_;
  std::vector<float> gen_eta_;
  std::vector<float> gen_phi_;
  std::vector<float> gen_mass_;
  std::vector<float> gen_vx_;
  std::vector<float> gen_vy_;
  std::vector<float> gen_vz_;
  std::vector<bool> gen_isHardProcess_;
  std::vector<int> gen_status_;
  std::vector<int> gen_moth_pdg_;
  std::vector<int> gen_daughter_n_;
  std::vector<int> gen_daughter_pdg_;
  
  unsigned int      muon_n_;
  std::vector<float> muon_pt_;
  std::vector<float> muon_p_;
  std::vector<float> muon_eta_;
  std::vector<float> muon_phi_;
  std::vector<float> muon_energy_;
  std::vector<float> muon_ptErr_;
  
  std::vector<bool>  muon_isLoose_;
  std::vector<bool>  muon_isMedium_;
  std::vector<bool>  muon_isTight_;
  std::vector<bool>  muon_isHighPtMuon_;
  std::vector<bool>  muon_isTrackerHighPtMuon_;
  std::vector<int>   muon_type_;
  std::vector<int>   muon_quality_;

  std::vector<bool> muon_hasMatchedGenTrack_;
  std::vector<float> muon_fromGenTrack_Pt_;
  std::vector<float> muon_fromGenTrack_PtErr_;
  std::vector<float> muon_fromGenTrack_Eta_;
  std::vector<float> muon_fromGenTrack_Phi_;

  std::vector<float> muon_d0_;
  std::vector<float> muon_d0Err_;
  std::vector<float> muon_charge_;
  std::vector<float> muon_dZ_;
  std::vector<float> muon_chi2_;

  std::vector<float> muon_pileupIso_;
  std::vector<float> muon_chargedIso_;
  std::vector<float> muon_photonIso_;
  std::vector<float> muon_neutralHadIso_;

  std::vector<float> muon_validFractionTrackerHits_;
  std::vector<float> muon_normChi2_;
  std::vector<float> muon_chi2LocalPosition_;
  std::vector<float> muon_kinkFinder_;
  std::vector<float> muon_segmentCompatability_;
  std::vector<float> muon_trkIso_;

  std::vector<float> muon_tuneP_Pt_;
  std::vector<float> muon_tuneP_PtErr_;
  std::vector<float> muon_tuneP_Eta_;
  std::vector<float> muon_tuneP_Phi_;
  std::vector<float> muon_tuneP_MuonBestTrackType_;
  
  std::vector<float> muon_avgEtaFromDTseg_;
  std::vector<float> muon_avgPhiFromDTseg_;
  
  std::vector<float> muon_comb_ndof_;
  std::vector<float> muon_comb_timeAtIpInOut_;
  std::vector<float> muon_comb_timeAtIpInOutErr_;
  std::vector<float> muon_comb_timeAtIpOutIn_;
  std::vector<float> muon_comb_timeAtIpOutInErr_;
  std::vector<float> muon_comb_invBeta_;
  std::vector<float> muon_comb_freeInvBeta_;
  std::vector<int>   muon_tofMap_found_;
  
  std::vector<int>   muon_dtSeg_n_;
  std::vector<int>   muon_dtSeg_found_;
  std::vector<float> muon_dtSeg_t0timing_;
  std::vector<float> muon_dtSeg_globX_;
  std::vector<float> muon_dtSeg_globY_;
  std::vector<float> muon_dtSeg_globZ_;
  std::vector<float> muon_dtSeg_eta_;
  std::vector<float> muon_dtSeg_phi_;

  unsigned int       simHit_n_;
  std::vector<float> simHit_globX_;
  std::vector<float> simHit_globY_;
  std::vector<float> simHit_globZ_;
  std::vector<float> simHit_tof_;

  std::vector<float> muon_r2_;
  std::vector<float> muon_dtSeg_rPhi_globY_; //dtGlobalPointYValues
  std::vector<float> muon_dtSeg_rPhi_t0timing_; //t0timingValues
  std::vector<float> muon_dtSeg_rPhi_t0timingCorrected_; //t0timingCorrectedValues
  std::vector<float> muon_dtSeg_rZ_globY_; //dtGlobalPointYValues_Ztiming
  std::vector<float> muon_dtSeg_rZ_t0timing_;

  std::vector<float> muon_rPhiSeg_correlationFactor_;
  std::vector<float> muon_rZSeg_correlationFactor_;
  
  std::vector<int> muon_dtSeg_Station_;
  std::vector<int> muon_dtSeg_Sector_;
  //std::vector<std::pair<int, int>> muon_dtSeg_rPhi_stationSector_;
  //std::vector<std::pair<int, int>> muon_dtSeg_rZ_stationSector_;

  // general tracks
  int                track_n_;
  std::vector<float> track_vx_;
  std::vector<float> track_vy_;
  std::vector<float> track_vz_;
  std::vector<float> track_phi_;
  std::vector<float> track_eta_;
  std::vector<float> track_pt_;
  std::vector<float> track_ptErr_;

};

//
// constructors and destructor
//
EarthAsDMAnalyzer::EarthAsDMAnalyzer(const edm::ParameterSet& iConfig) :
 verbose_(iConfig.getUntrackedParameter<int>("verbosityLevel")),
 isData_(iConfig.getUntrackedParameter<int>("isData")),
 hasSim_(iConfig.getUntrackedParameter<int>("hasSim")),
 genParticlesToken_(consumes< std::vector<reco::GenParticle> >( edm::InputTag("genParticles") )),
 muonToken_(consumes<std::vector<reco::Muon>>(iConfig.getParameter<edm::InputTag>("muonCollection"))),
 muonTimeToken_(consumes<reco::MuonTimeExtraMap>(iConfig.getParameter<edm::InputTag>("muonTimeCollection"))),
 PSimHitContainerToken_(consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("PSimHitContainer"))),
 cscSegmentToken_(consumes< CSCSegmentCollection >( edm::InputTag("cscSegments") )),
 dtSegmentToken_(consumes< DTRecSegment4DCollection >( edm::InputTag("dt4DSegments") )),
 offlinePrimaryVerticesToken_(consumes<std::vector<reco::Vertex>>(edm::InputTag("offlinePrimaryVertices"))),
 muonDTGeomToken_(esConsumes()),
 t0Token_(esConsumes()),
 ttrigToken_(esConsumes()),
 triggerResultsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
 tracksToken_(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("trackCollection")))
{}


EarthAsDMAnalyzer::~EarthAsDMAnalyzer()  = default;

//
// member functions
//

  float mean(const std::vector<float>& data);
  float calculateRSquared(const std::vector<float>& x, const std::vector<float>& y);
  float pearsonCorrelation(const std::vector<float>& x, const std::vector<float>& y);


// ------------ method called for each event  ------------
void EarthAsDMAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace std;
  using namespace edm;
  static constexpr const char* const MOD = "Analyzer";
  
  runNumber_ = iEvent.id().run();
  lsNumber_ = iEvent.id().luminosityBlock();
  eventNumber_ = iEvent.id().event();
  
  if (verbose_ > 1) LogPrint(MOD) << "\n\n---------------------------------------------------------------\nAnalyzing runNumber "
    << runNumber_ << " lsNumber " << lsNumber_ << " eventNumber " << eventNumber_;


  //------------------------------------------------------------------
  // Necessary handles to be used
  //------------------------------------------------------------------
  
  edm::Handle<reco::MuonCollection> muonCollectionHandle;
  iEvent.getByToken(muonToken_,muonCollectionHandle);

  //vector<reco::Muon> muonColl = iEvent.get(muonToken_);
  
  edm::Handle<reco::MuonTimeExtraMap> tofMap;
  iEvent.getByToken(muonTimeToken_, tofMap);
  
  edm::Handle<CSCSegmentCollection> cscSegments;
  iEvent.getByToken(cscSegmentToken_, cscSegments);

  edm::Handle<std::vector<reco::Track>> trackCollectionHandle;
  iEvent.getByToken(tracksToken_, trackCollectionHandle);
  
//  edm::Handle<DTRecSegment4DCollection> dtSegments;
//  iEvent.getByToken(dtSegmentToken_, dtSegments);
  
  muonDTGeom = &iSetup.getData(muonDTGeomToken_);

  ESHandle<DTT0> t0H;
  t0H = iSetup.getHandle(t0Token_);
  t0Map = &iSetup.getData(t0Token_);

  ESHandle<DTTtrig> tTrigH;
  tTrigH = iSetup.getHandle(ttrigToken_);
  tTrigMap = &iSetup.getData(ttrigToken_);
  
  //------------------------------------------------------------------
  // SimHits saved in the ntuple
  //------------------------------------------------------------------

  if (hasSim_) {
    edm::Handle<edm::PSimHitContainer> SimHitCollection;
    PSimHitContainer simhitTC;

    iEvent.getByToken(PSimHitContainerToken_, SimHitCollection);
    simhitTC = *(SimHitCollection.product());
    
    simHit_n_ = 0;
    for (PSimHitContainer::const_iterator simHit=simhitTC.begin(); simHit!=simhitTC.end(); simHit++){
      float simHitGlobalPointX = 9999;
      float simHitGlobalPointY = 9999;
      float simHitGlobalPointZ = 9999;

      LocalPoint simHitLocalPosition = simHit->localPosition();
      const GeomDet* simHitDet = muonDTGeom->idToDet(simHit->detUnitId());
      GlobalPoint simHitglobalPoint = simHitDet->toGlobal(simHitLocalPosition);
      simHitGlobalPointX = simHitglobalPoint.x();
      simHitGlobalPointY = simHitglobalPoint.y();
      simHitGlobalPointZ = simHitglobalPoint.z();

      simHit_globX_.push_back( simHitGlobalPointX);
      simHit_globY_.push_back( simHitGlobalPointY);
      simHit_globZ_.push_back( simHitGlobalPointZ);
      simHit_tof_.push_back(simHit->tof());
      simHit_n_++;
    }
  }

  //------------------------------------------------------------------
  // GenParticles to be saved in the ntuple
  //------------------------------------------------------------------
  edm::Handle<reco::VertexCollection> vertexColl;
  iEvent.getByToken(offlinePrimaryVerticesToken_, vertexColl);
  
  edm::Handle< std::vector<reco::GenParticle> > genColl;

  if (!isData_) {
    iEvent.getByToken(genParticlesToken_, genColl);
    gen_n_ = genColl->size();
    if (verbose_ > 2) LogPrint(MOD) << "The stable GenCadidate has PDG ID = " <<  (*genColl)[1].pdgId() << " and pT = " << (*genColl)[1].pt()
      << " eta = " << (*genColl)[1].eta() << " phi = "  << (*genColl)[1].phi() << " status = " << (*genColl)[1].status();
    for (const auto& genCand : *genColl) {
      gen_pdg_.push_back(genCand.pdgId());
      gen_energy_.push_back(genCand.energy());
      gen_pt_.push_back(genCand.pt());
      gen_eta_.push_back(genCand.eta());
      gen_phi_.push_back(genCand.phi());
      gen_mass_.push_back(genCand.mass());
      gen_vx_.push_back(genCand.vx());
      gen_vy_.push_back(genCand.vy());
      gen_vz_.push_back(genCand.vz());
      gen_isHardProcess_.push_back(genCand.isHardProcess());
      gen_status_.push_back(genCand.status());
      if (genCand.numberOfMothers() > 0) {
        gen_moth_pdg_.push_back(genCand.mother()->pdgId());
      } else {
        gen_moth_pdg_.push_back(0);
      }
      gen_daughter_n_.push_back(genCand.numberOfDaughters());
      if (genCand.numberOfDaughters() > 0) {
        gen_daughter_pdg_.push_back(genCand.daughter(0)->pdgId());
      } else {
        gen_daughter_pdg_.push_back(0);
      }
    }
  }
  
  const reco::Vertex& highestSumPt2Vertex = vertexColl->front();

  //------------------------------------------------------------------
  // Get trigger results for this event
  //------------------------------------------------------------------
  const edm::Handle<edm::TriggerResults> triggerH = iEvent.getHandle(triggerResultsToken_);
  const auto triggerNames = iEvent.triggerNames(*triggerH);

  //------------------------------------------------------------------
  // Save trigger decisions in array of booleans
  //------------------------------------------------------------------
  
  HLT_L1SingleMu3_                 = false;
  HLT_L1SingleMu5_                 = false;
  HLT_L1SingleMu7_                 = false;
  HLT_L1SingleMu18_                = false;
  HLT_L2Mu23NoVtx_2Cha_CosmicSeed_ = false;
  HLT_L1SingleMuOpen_              = false;
  HLT_L1SingleMuOpen_DT_           = false;
  HLT_L1SingleMuCosmics_           = false;

  if (verbose_ > 5) LogPrint(MOD) << "The following triggers pass in this event: ";
  for (unsigned int i = 0; i < triggerH->size(); i++) {
    if (TString(triggerNames.triggerName(i)).Contains("HLT_") && triggerH->accept(i)) {
      if (verbose_ > 5) cout << "  * " <<  triggerNames.triggerName(i) << endl;
    }
    if (TString(triggerNames.triggerName(i)).Contains("HLT_L1SingleMu3_v") && triggerH->accept(i))
      HLT_L1SingleMu3_ = true;
    if (TString(triggerNames.triggerName(i)).Contains("HLT_L1SingleMu5_v") && triggerH->accept(i))
      HLT_L1SingleMu5_ = true;
    if (TString(triggerNames.triggerName(i)).Contains("HLT_L1SingleMu7_v") && triggerH->accept(i))
      HLT_L1SingleMu7_ = true;
    if (TString(triggerNames.triggerName(i)).Contains("HLT_L1SingleMu18_v") && triggerH->accept(i))
      HLT_L1SingleMu18_ = true;
    if (TString(triggerNames.triggerName(i)).Contains("HLT_L2Mu23NoVtx_2Cha_CosmicSeed") && triggerH->accept(i))
      HLT_L2Mu23NoVtx_2Cha_CosmicSeed_ = true;
    if (TString(triggerNames.triggerName(i)).Contains("HLT_L1SingleMuOpen_v") && triggerH->accept(i))
      HLT_L1SingleMuOpen_ = true;
    if (TString(triggerNames.triggerName(i)).Contains("HLT_L1SingleMuOpen_DT_v") && triggerH->accept(i))
      HLT_L1SingleMuOpen_DT_ = true;
    if (TString(triggerNames.triggerName(i)).Contains("HLT_L1SingleMuCosmics_v") && triggerH->accept(i))
      HLT_L1SingleMuCosmics_ = true;
          
  }

  //------------------------------------------------------------------
  // Loop on the muon collection
  //------------------------------------------------------------------

  muon_n_ = 0;
  for (unsigned int i = 0; i < muonCollectionHandle->size(); i++) {
    if (verbose_ > 2) LogPrint(MOD) << "\n  Analyzing track " << i ;
    reco::MuonRef muon  = reco::MuonRef( muonCollectionHandle, i );
    muon_pt_.push_back( muon->pt());
    // muon_ptErr_.push_back( muon->ptError());  // reco::Muon doesn't have func ptErr()
    muon_p_.push_back( muon->p());
    muon_eta_.push_back( muon->eta());
    muon_phi_.push_back( muon->phi());
    muon_energy_.push_back( muon->energy());
    muon_charge_.push_back( muon->charge());
    if (muon->track().isNonnull()) {
        // The muon has a valid track
      muon_chi2_.push_back(muon->track()->chi2());
    // Add other track-related data here
    } else {
        // Handle the case where the track is null or invalid
        // You might want to store a default or special value
        muon_chi2_.push_back(0.0); // Replace with an appropriate default value
    }

    bool isMedium = muon::isMediumMuon(*muon);
    muon_isMedium_.push_back( isMedium);

    bool isTight = muon::isTightMuon(*muon, highestSumPt2Vertex);
    muon_isTight_.push_back( isTight);

    bool isHighPt = muon::isHighPtMuon(*muon, highestSumPt2Vertex);
    muon_isHighPtMuon_.push_back( isHighPt);

    bool isTrackerHighPtMuon = muon::isTrackerHighPtMuon(*muon, highestSumPt2Vertex);
    muon_isTrackerHighPtMuon_.push_back( isTrackerHighPtMuon);

    // Match muon inner tracks with general track
    bool hasMatchedGenTrack = 0;
    float genTrackPt = -9999.; 
    float genTrackPtErr = -9999.; 
    float genTrackEta = -9999.;
    float genTrackPhi = -9999.;
    reco::TrackRef innertrack = muon->innerTrack();
    if(innertrack.isNonnull()) {
      if (verbose_ > 2) LogPrint(MOD) << "  >> muon inner track exists: key " << innertrack.key()
                                      << " pT " << innertrack->pt()
                                      << " eta " << innertrack->eta()
                                      << " phi " << innertrack->phi();
      // loop on the general track
      for(unsigned int c=0;c<trackCollectionHandle->size();c++) {
        reco::TrackRef genTrackRef = reco::TrackRef( trackCollectionHandle.product(), c );
        // match keys
        if (genTrackRef.key() == innertrack.key()) {
          if (verbose_ > 2) LogPrint(MOD) << "    >> general track matches: key " << genTrackRef.key()
                                          << " pT " << genTrackRef->pt()
                                          << " eta " << genTrackRef->eta()
                                          << " phi " << genTrackRef->phi();
          hasMatchedGenTrack = 1;
          genTrackPt = genTrackRef->pt();
          genTrackPtErr = genTrackRef->ptError();
          genTrackEta = genTrackRef->eta();
          genTrackPhi = genTrackRef->phi();
        }
      } // end loop on general track
    }
    muon_hasMatchedGenTrack_.push_back( hasMatchedGenTrack );
    muon_fromGenTrack_Pt_.push_back( genTrackPt );
    muon_fromGenTrack_PtErr_.push_back( genTrackPtErr );
    muon_fromGenTrack_Eta_.push_back( genTrackEta );
    muon_fromGenTrack_Phi_.push_back( genTrackPhi );

    // Calculate the transverse impact parameter (d0) of the muon with respect to the highest sum pt vertex
    // TAV: FYI this has no meaning in cosmics
    const reco::TrackRef bestTrack = muon->muonBestTrack();
    muon_d0_.push_back( bestTrack->dxy(highestSumPt2Vertex.position()));
    muon_d0Err_.push_back( bestTrack->dxyError());
    muon_dZ_.push_back( bestTrack->dz(highestSumPt2Vertex.position()));

    muon_validFractionTrackerHits_.push_back( (muon->innerTrack().isNonnull() ? muon->track()->validFraction() : -99.0));

    muon_pileupIso_.push_back( muon->pfIsolationR04().sumPUPt);
    muon_chargedIso_.push_back( muon->pfIsolationR04().sumChargedHadronPt);
    muon_photonIso_.push_back( muon->pfIsolationR04().sumPhotonEt);
    muon_neutralHadIso_.push_back( muon->pfIsolationR04().sumNeutralHadronEt);
    muon_trkIso_.push_back( muon->isolationR03().sumPt);
    muon_chi2LocalPosition_.push_back( muon->combinedQuality().chi2LocalPosition);
    muon_kinkFinder_.push_back( muon->combinedQuality().trkKink);

    muon_type_.push_back((muon->isMuon() + 2*muon->isGlobalMuon() + 4*muon->isTrackerMuon() + 8*muon->isStandAloneMuon()
        + 16*muon->isCaloMuon() + 32*muon->isPFMuon() + 64*muon->isRPCMuon()));


    muon_tuneP_Pt_.push_back( muon->tunePMuonBestTrack()->pt());
    muon_tuneP_PtErr_.push_back( muon->tunePMuonBestTrack()->ptError());
    muon_ptErr_.push_back( muon->tunePMuonBestTrack()->ptError() );  // muons uses tuneP reco algo
    muon_tuneP_Eta_.push_back( muon->tunePMuonBestTrack()->eta());
    muon_tuneP_Phi_.push_back( muon->tunePMuonBestTrack()->phi());
    muon_tuneP_MuonBestTrackType_.push_back( muon->tunePMuonBestTrackType());
    muon_segmentCompatability_.push_back( (muon::segmentCompatibility(*muon)));

    muon_quality_.push_back((
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
      )));

    muon_normChi2_.push_back( (muon::isGoodMuon(*muon,muon::AllGlobalMuons) ? muon->globalTrack()->normalizedChi2() : -99.0));
    
    
    if (verbose_ > 2) LogPrint(MOD) << "  >> muon_pt_ " << muon->pt() << " muon_eta_ " << muon->eta() << " muon_phi_ " << muon->phi();
    if (verbose_ > 2) LogPrint(MOD) << "  >> muon_vx_ " << muon->vx() << " muon_vy_ " << muon->vy() << " muon_vz_ " << muon->vz();

        
//    const reco::MuonTime time = muon->time();
//    const reco::MuonTime rpcTime = muon->rpcTime();
    float muonAvgEtaFromDTseg = 9999.;
    float muonAvgPhiFromDTseg = 9999.;
    float muonSumEtaFromDTseg = 0.;
    float muonSumPhiFromDTseg = 0.;
    int dtChamb_n = 0;
    int dtSeg_n = 0;
    if (muon->isMatchesValid()) {
      // Loop on the chambers belonging to this muon
      std::vector<reco::MuonChamberMatch>::const_iterator chamberMatch;
      for ( chamberMatch = muon->matches().begin(); chamberMatch != muon->matches().end(); ++chamberMatch) {
        if (verbose_ > 4) LogPrint(MOD)  << "    >> Chamber index " << dtChamb_n;
        const vector<reco::MuonSegmentMatch> matchedSegments = chamberMatch->segmentMatches;
        vector<reco::MuonSegmentMatch>::const_iterator segment;
        // Now loop on the segments in the chamber
        // Each DT chamber is made of three (or two in MB4) Superlayers (SL)
        // each SL consisting of four layers of rectangular drift cells staggered by half a tube width
        // SL1 and SL3 measure the RPhi coordinate
        // SL2 measures the Z coordinate
        for (segment = matchedSegments.begin(); segment != matchedSegments.end(); ++segment) {
          edm::Ref<DTRecSegment4DCollection> dtSegment = segment->dtSegmentRef;
          int found = 0;
          float t0timing = 9999;
          float t0timingZed = 9999;
          float dtGlobalPointX = 9999;
          float dtGlobalPointY = 9999;
          float dtGlobalPointZ = 9999;
          float dtGlobalPointEta = 9999;
          float dtGlobalPointPhi = 9999;
          float station = 9999;
          float sector = 9999;
          
          if (!dtSegment.isNull()) {
            found = 1;
            if (verbose_ > 4) LogPrint(MOD)  << "      >> DT segment index " << dtSeg_n;
            LocalPoint segmentLocalPosition = dtSegment->localPosition();
            // Lets check if it has R-Phi segment
            if (dtSegment->hasPhi()) {
              const auto& dtPhiSegment = dtSegment->phiSegment();
              t0timing = dtPhiSegment->t0();
            } else {
              if (verbose_ > 5) LogPrint(MOD) << "        >> This 4D segment does not have a phi segment: ";
            }
            // Now let's look at the Zed
            if (dtSegment->hasZed()) {
              const auto& dtZedSegment = dtSegment->zSegment();
              t0timingZed = dtZedSegment->t0();
            } else {
              if (verbose_ > 5) LogPrint(MOD) << "        >> It does not it has a zed segment: ";
            }
            const GeomDet* dtDet = muonDTGeom->idToDet(dtSegment->geographicalId());
            DTChamberId dtChamberId(dtSegment->geographicalId());
            DTLayerId dtLayerId((dtSegment->geographicalId()));
            if (verbose_ > 4) LogPrint(MOD) << "        >> Wheel /  Station / Sector / Layer " << dtChamberId.wheel()
              << " / " << dtChamberId.station() << " / "  << dtChamberId.sector() << " / " << dtLayerId.layer();
            
            station = dtChamberId.station();
            sector = dtChamberId.sector();

            // Global Point coordinates
            GlobalPoint globalPoint = dtDet->toGlobal(segmentLocalPosition);
            dtGlobalPointX = globalPoint.x();
            dtGlobalPointY = globalPoint.y();
            dtGlobalPointZ = globalPoint.z();
            dtGlobalPointEta = globalPoint.eta();
            dtGlobalPointPhi = globalPoint.phi();
            if (verbose_ > 4) LogPrint(MOD) << "        >> t0timing (t0timingZed): " << t0timing << " (" << t0timingZed << ")";
            if (verbose_ > 4) LogPrint(MOD) << "        >> globalY = " << dtGlobalPointY;
            if (verbose_ > 5) LogPrint(MOD) << "        >> eta = " << dtGlobalPointEta << " phi = " << dtGlobalPointPhi;
            muonSumEtaFromDTseg += dtGlobalPointEta;
            muonSumPhiFromDTseg += dtGlobalPointPhi;
            muon_dtSeg_rPhi_t0timing_.push_back(t0timing);
            muon_dtSeg_rPhi_globY_.push_back(dtGlobalPointY);
            muon_dtSeg_rPhi_t0timingCorrected_.push_back(t0timing);
            muon_dtSeg_rZ_t0timing_.push_back(t0timingZed); // HERE
            muon_dtSeg_rZ_globY_.push_back(dtGlobalPointY);
            for (size_t i = 0; i < muon_dtSeg_rPhi_globY_.size(); ++i) {
              // cout << "dtGlobalPointYValues[" << i << "]: " << dtGlobalPointYValues[i] << endl;
              if (verbose_ > 2) LogPrint(MOD) << "muon_dtSeg_rPhi_globY[" << i << "]: " << muon_dtSeg_rPhi_globY_[i] << endl;
            }

            for (size_t i = 0; i < muon_dtSeg_rPhi_t0timing_.size(); ++i) {
              // cout << "t0timingValues[" << i << "]: " << t0timingValues[i] << endl;
              if (verbose_ > 2) LogPrint(MOD) << "muon_dtSeg_rPhi_t0timing[" << i << "]: " << muon_dtSeg_rPhi_t0timing_[i] << endl;
            }

            dtSeg_n++;
          }
          muon_dtSeg_found_.push_back( found);
          muon_dtSeg_t0timing_.push_back( t0timing);
          muon_dtSeg_globX_.push_back( dtGlobalPointX);
          muon_dtSeg_globY_.push_back( dtGlobalPointY);
          muon_dtSeg_globZ_.push_back( dtGlobalPointZ);
          muon_dtSeg_eta_.push_back( dtGlobalPointEta);
          muon_dtSeg_phi_.push_back( dtGlobalPointPhi);
          muon_dtSeg_Station_.push_back( station);
          muon_dtSeg_Sector_.push_back( sector);
        } // end loop on segments
        dtChamb_n++;
      } // end loop on chamber matches
      muon_dtSeg_n_.push_back( dtSeg_n);
      if (verbose_ > 3) LogPrint(MOD)  << "  >> This track had " << dtSeg_n << " segments";
    } // end condition on muon having valid match
    // muon_dtSeg_globX_.push_back( muon_dtSeg_t0timing_ ) ;
    if (dtSeg_n > 0) {
      muonAvgEtaFromDTseg = muonSumEtaFromDTseg / dtSeg_n;
      muonAvgPhiFromDTseg = muonSumPhiFromDTseg / dtSeg_n;
      
      float rSquared = calculateRSquared( muon_dtSeg_rPhi_globY_, muon_dtSeg_rPhi_t0timingCorrected_);
      // Store the R-squared value in the rSquaredValues vector
      muon_r2_.push_back( rSquared);

      float PearsonCorrelation = pearsonCorrelation(muon_dtSeg_rPhi_globY_, muon_dtSeg_rPhi_t0timingCorrected_);
        muon_rPhiSeg_correlationFactor_.push_back( PearsonCorrelation);

      // for (size_t i = 0; i < muon_dtSeg_Station_.size(); ++i) {
      //   muon_dtSeg_rPhi_stationSector_.push_back(std::make_pair(muon_dtSeg_Station_[i], muon_dtSeg_Sector_[i]));
      // }
      // Print the combined vector
      // cout << "Station and Sector Vector: ";
      // int muonIndex = 1; // Initialize the Muon index
      // if (verbose_ > 16) LogPrint(MOD)  << "Muon rPhi Station and Sector: ";
      // for (const auto& pair : muon_dtSeg_rPhi_stationSector_) {
      //   if (verbose_ > 16) { LogPrint(MOD) << "Muon " << muonIndex << " Station " << pair.first << " in Sector " << pair.second << " ";
      //   }
      //   muonIndex++;
      // }
      // cout << endl;
      }



    float PearsonCorrelation_Z = pearsonCorrelation(muon_dtSeg_rZ_globY_, muon_dtSeg_rZ_t0timing_);
      muon_rZSeg_correlationFactor_.push_back( PearsonCorrelation_Z);


      // for (size_t i = 0; i < muon_dtSeg_Station_.size(); ++i) {
      //   muon_dtSeg_rZ_stationSector_.push_back(std::make_pair(muon_dtSeg_Station_[i], muon_dtSeg_Sector_[i]));
      // }
      // // Print the combined vector
      // if (verbose_ > 16) LogPrint(MOD)  << "Muon rZ Station and Sector: ";
      // int muonIndex = 1;
      // for (const auto& pair : muon_dtSeg_rZ_stationSector_) {
      //   if (verbose_ > 16) { LogPrint(MOD) << "Muon " << muonIndex   << " Station " << pair.first << " in Sector " << pair.second << " ";
      //   }
      //   muonIndex++;
      // }

      cout << endl;

      if (verbose_ > 2) LogPrint(MOD) << "muon_rPhiSeg_correlationFactor:";
      string PearsonValuesStr;
      for (float value : muon_rPhiSeg_correlationFactor_) {
          PearsonValuesStr += to_string(value) + " ";
      }
      if (verbose_ > 2) LogPrint(MOD) << PearsonValuesStr << endl;



      if (verbose_ > 2) LogPrint(MOD) << "muon_rZSeg_correlationFactor:";
      string PearsonValues_Z_Str;
      for (float value : muon_rZSeg_correlationFactor_) {
          PearsonValues_Z_Str += to_string(value) + " ";
      }
      if (verbose_ > 2) LogPrint(MOD) << PearsonValues_Z_Str << endl;


      if (verbose_ > 2) LogPrint(MOD) << "muon_r2:";
      string rSquaredStr;
      for (float value : muon_r2_) {
          rSquaredStr += to_string(value) + " ";
      }
      if (verbose_ > 2) LogPrint(MOD) << rSquaredStr << endl;



      muon_dtSeg_rPhi_globY_.clear();
      muon_dtSeg_rPhi_t0timing_.clear();
      muon_dtSeg_rPhi_t0timingCorrected_.clear();
      muon_dtSeg_rZ_globY_.clear();
      muon_dtSeg_rZ_t0timing_.clear();

    
    if (!isData_) {
      float dRGenMuonFromAvg = deltaR((*genColl)[1].eta(),(*genColl)[1].phi(),muonAvgEtaFromDTseg,muonAvgPhiFromDTseg);
      if (verbose_ > 3) LogPrint(MOD)  << "  >> muonAvgEtaFromDTseg = " << muonAvgEtaFromDTseg << " muonAvgPhiFromDTseg = "  << muonAvgPhiFromDTseg << " dRGenMuonFromAvg = " << dRGenMuonFromAvg;
    }
    muon_avgEtaFromDTseg_.push_back(muonAvgEtaFromDTseg);
    muon_avgPhiFromDTseg_.push_back(muonAvgPhiFromDTseg);
    
    if (tofMap.isValid()) {
      const reco::MuonTimeExtra* combinedTimeExtra = NULL;
      float nDof = 9999;
      float timeAtIpInOut = 9999;
      float timeAtIpInOutErr = 9999;
      float timeAtIpOutIn = 9999;
      float timeAtIpOutInErr = 9999;
      float invBeta = 9999;
      float freeInvBeta = 9999;

      combinedTimeExtra = &tofMap->get(muon.key());
      if (combinedTimeExtra != NULL) {
        nDof = combinedTimeExtra->nDof();
        muon_tofMap_found_.push_back(1);
        timeAtIpInOut = combinedTimeExtra->timeAtIpInOut();
        timeAtIpInOutErr = combinedTimeExtra->timeAtIpInOutErr();
        timeAtIpOutIn = combinedTimeExtra->timeAtIpOutIn();
        timeAtIpOutInErr = combinedTimeExtra->timeAtIpOutInErr();
        invBeta = combinedTimeExtra->inverseBeta();
        freeInvBeta = combinedTimeExtra->freeInverseBeta();

        // Sign convention for muonCombinedFreeInvBeta_:
        //   positive - outward moving particle
        //   negative - inward moving particle
        if (verbose_ > 4) LogPrint(MOD) << "combinedTimeExtra->nDof() " << combinedTimeExtra->nDof() << "  timeAtIpInOut_ " << combinedTimeExtra->timeAtIpInOut()  << " +/- " <<  combinedTimeExtra->timeAtIpOutInErr() << " timeAtIpOutIn_ " << combinedTimeExtra->timeAtIpOutIn() << " +/- " <<  combinedTimeExtra->timeAtIpInOutErr() << " invBeta_ " << combinedTimeExtra->inverseBeta() << " freeInvBeta_ " << combinedTimeExtra->freeInverseBeta();
      }
      else {
        muon_tofMap_found_.push_back(0);
      }
      muon_comb_ndof_.push_back(nDof);
      muon_comb_timeAtIpInOut_.push_back(timeAtIpInOut);
      muon_comb_timeAtIpInOutErr_.push_back(timeAtIpInOutErr);
      muon_comb_timeAtIpOutIn_.push_back(timeAtIpOutIn);
      muon_comb_timeAtIpOutInErr_.push_back(timeAtIpOutInErr);
      muon_comb_invBeta_.push_back(invBeta);
      muon_comb_freeInvBeta_.push_back(freeInvBeta);
    }    
    muon_n_++ ;
  } // end of muon collection loop

  //------------------------------------------------------------------
  // General track collection
  //------------------------------------------------------------------
  
  track_n_ = 0;
  for (const auto& track : *trackCollectionHandle) {
    track_vx_.push_back(track.vx());
    track_vy_.push_back(track.vy());
    track_vz_.push_back(track.vz());
    track_phi_.push_back(track.phi());
    track_eta_.push_back(track.eta());
    track_pt_.push_back(track.pt());
    track_ptErr_.push_back(track.ptError());

    track_n_++;
    
  }

  // Fill the tree
  outputTree_->Fill();
  
  // Clear the vectors
  gen_pdg_.clear();
  gen_energy_.clear();
  gen_pt_.clear();
  gen_eta_.clear();
  gen_phi_.clear();
  gen_mass_.clear();
  gen_vx_.clear();
  gen_vy_.clear();
  gen_vz_.clear();
  gen_isHardProcess_.clear();
  gen_status_.clear();
  gen_moth_pdg_.clear();
  gen_daughter_n_.clear();
  gen_daughter_pdg_.clear();
  
  muon_pt_.clear();
  muon_ptErr_.clear();
  muon_p_.clear();
  muon_eta_.clear();
  muon_phi_.clear();
  muon_energy_.clear();
  muon_charge_.clear();
  
  muon_isLoose_.clear();
  muon_isMedium_.clear();
  muon_isTight_.clear();
  muon_isHighPtMuon_.clear();
  muon_isTrackerHighPtMuon_.clear();
  muon_type_.clear();
  muon_quality_.clear();
  muon_chi2_.clear();

  muon_hasMatchedGenTrack_.clear();
  muon_fromGenTrack_Pt_.clear();
  muon_fromGenTrack_PtErr_.clear();
  muon_fromGenTrack_Eta_.clear();
  muon_fromGenTrack_Phi_.clear();
  
  muon_d0_.clear();
  muon_d0Err_.clear();
  muon_charge_.clear();
  muon_dZ_.clear();
  
  muon_pileupIso_.clear();
  muon_chargedIso_.clear();
  muon_photonIso_.clear();
  muon_neutralHadIso_.clear();
  
  muon_validFractionTrackerHits_.clear();
  muon_normChi2_.clear();
  muon_chi2LocalPosition_.clear();
  muon_kinkFinder_.clear();
  muon_segmentCompatability_.clear();
  muon_trkIso_.clear();
  
  muon_tuneP_Pt_.clear();
  muon_tuneP_PtErr_.clear();
  muon_tuneP_Eta_.clear();
  muon_tuneP_Phi_.clear();
  muon_tuneP_MuonBestTrackType_.clear();
  
  muon_avgEtaFromDTseg_.clear();
  muon_avgPhiFromDTseg_.clear();
  
  muon_comb_ndof_.clear();
  muon_comb_timeAtIpInOut_.clear();
  muon_comb_timeAtIpInOutErr_.clear();
  muon_comb_timeAtIpOutIn_.clear();
  muon_comb_timeAtIpOutInErr_.clear();
  muon_comb_invBeta_.clear();
  muon_comb_freeInvBeta_.clear();
  muon_tofMap_found_.clear();
  
  muon_dtSeg_n_.clear();
  muon_dtSeg_found_.clear();
  muon_dtSeg_t0timing_.clear();
  muon_dtSeg_globX_.clear();
  muon_dtSeg_globY_.clear();
  muon_dtSeg_globZ_.clear();
  muon_dtSeg_eta_.clear();
  muon_dtSeg_phi_.clear();

  simHit_globX_.clear();
  simHit_globY_.clear();
  simHit_globZ_.clear();
  simHit_tof_.clear();
  
  muon_r2_.clear();
  muon_rPhiSeg_correlationFactor_.clear();
  muon_rZSeg_correlationFactor_.clear();

  track_vx_.clear();
  track_vy_.clear();
  track_vz_.clear();
  track_phi_.clear();
  track_eta_.clear();
  track_pt_.clear();
  track_ptErr_.clear();

  muon_dtSeg_Station_.clear();
  muon_dtSeg_Sector_.clear();
  //muon_dtSeg_rPhi_stationSector_.clear();
  //muon_dtSeg_rZ_stationSector_.clear();
  
}

// ------------ method called once each job just before starting event loop  ------------
void EarthAsDMAnalyzer::beginJob() {
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  outputTree_ = fs->make<TTree>("tree", "tree");
  
  // Setup the branches of the the tree
  outputTree_->Branch("run",    &runNumber_);
  outputTree_->Branch("ls",     &lsNumber_);
  outputTree_->Branch("event",  &eventNumber_);
  
  outputTree_ -> Branch ( "HLT_L1SingleMu3",       &HLT_L1SingleMu3_) ;
  outputTree_ -> Branch ( "HLT_L1SingleMu5",       &HLT_L1SingleMu5_) ;
  outputTree_ -> Branch ( "HLT_L1SingleMu7",       &HLT_L1SingleMu7_) ;
  outputTree_ -> Branch ( "HLT_L1SingleMu18",      &HLT_L1SingleMu18_) ;
  outputTree_ -> Branch ( "HLT_L2Mu23NoVtx_2Cha_CosmicSeed",
                                                   &HLT_L2Mu23NoVtx_2Cha_CosmicSeed_) ;
  outputTree_ -> Branch ( "HLT_L1SingleMuOpen",    &HLT_L1SingleMuOpen_) ;
  outputTree_ -> Branch ( "HLT_L1SingleMuOpen_DT", &HLT_L1SingleMuOpen_DT_) ;
  outputTree_ -> Branch ( "HLT_L1SingleMuCosmics", &HLT_L1SingleMuCosmics_) ;
  
  outputTree_ -> Branch ( "gen_n",            &gen_n_) ;
  outputTree_ -> Branch ( "gen_pdg",          &gen_pdg_);
  outputTree_ -> Branch ( "gen_pt",           &gen_pt_);
  outputTree_ -> Branch ( "gen_eta",          &gen_eta_);
  outputTree_ -> Branch ( "gen_phi",          &gen_phi_);
  outputTree_ -> Branch ( "gen_mass",         &gen_mass_);
  outputTree_ -> Branch ( "gen_vx",           &gen_vx_);
  outputTree_ -> Branch ( "gen_vy",           &gen_vy_);
  outputTree_ -> Branch ( "gen_vz",           &gen_vz_);
  outputTree_ -> Branch ( "gen_isHardProcess",&gen_isHardProcess_);
  outputTree_ -> Branch ( "gen_status",       &gen_status_);
  outputTree_ -> Branch ( "gen_moth_pdg",     &gen_moth_pdg_);
  outputTree_ -> Branch ( "gen_daughter_n",   &gen_daughter_n_);
  outputTree_ -> Branch ( "gen_daughter_pdg", &gen_daughter_pdg_);
  
  outputTree_ -> Branch ( "muon_n",           &muon_n_);
  outputTree_ -> Branch ( "muon_pt",          &muon_pt_);
  outputTree_ -> Branch ( "muon_ptErr",       &muon_ptErr_);
  outputTree_ -> Branch ( "muon_p",           &muon_p_);
  outputTree_ -> Branch ( "muon_eta",         &muon_eta_);
  outputTree_ -> Branch ( "muon_phi",         &muon_phi_);
  outputTree_ -> Branch ( "muon_energy",      &muon_energy_);
  outputTree_ -> Branch ( "muon_charge",      &muon_charge_);
  outputTree_ -> Branch ( "muon_chi2",        &muon_chi2_);
  outputTree_ -> Branch ( "muon_r2",           &muon_r2_);



  outputTree_ -> Branch ( "muon_isLoose",     &muon_isLoose_);
  outputTree_ -> Branch ( "muon_isMedium",    &muon_isMedium_);
  outputTree_ -> Branch ( "muon_isTight",     &muon_isTight_);
  outputTree_ -> Branch ( "muon_isTrackerHighPtMuon",
                                              &muon_isTrackerHighPtMuon_);
  outputTree_ -> Branch ( "muon_isHighPtMuon",&muon_isHighPtMuon_);
  outputTree_ -> Branch ( "muon_type",        &muon_type_);
  outputTree_ -> Branch ( "muon_quality",     &muon_quality_);

  outputTree_ -> Branch ( "muon_hasMatchedGenTrack", &muon_hasMatchedGenTrack_);
  outputTree_ -> Branch ( "muon_fromGenTrack_Pt",    &muon_fromGenTrack_Pt_);
  outputTree_ -> Branch ( "muon_fromGenTrack_PtErr", &muon_fromGenTrack_PtErr_);
  outputTree_ -> Branch ( "muon_fromGenTrack_Eta",   &muon_fromGenTrack_Eta_);
  outputTree_ -> Branch ( "muon_fromGenTrack_Phi",   &muon_fromGenTrack_Phi_);

  outputTree_ -> Branch ( "muon_d0",          &muon_d0_);
  outputTree_ -> Branch ( "muon_d0Err",       &muon_d0Err_);
  outputTree_ -> Branch ( "muon_dZ",          &muon_dZ_);

  outputTree_ -> Branch ( "muon_pileupIso",     &muon_pileupIso_);
  outputTree_ -> Branch ( "muon_chargedIso",    &muon_chargedIso_);
  outputTree_ -> Branch ( "muon_photonIso",     &muon_photonIso_);
  outputTree_ -> Branch ( "muon_neutralHadIso", &muon_neutralHadIso_);
  outputTree_ -> Branch ( "muon_validFractionTrackerHits",
                                                &muon_validFractionTrackerHits_);

  outputTree_ -> Branch ( "muon_tuneP_Pt",      &muon_tuneP_Pt_);
  outputTree_ -> Branch ( "muon_tuneP_PtErr",   &muon_tuneP_PtErr_);
  outputTree_ -> Branch ( "muon_tuneP_Eta",     &muon_tuneP_Eta_);
  outputTree_ -> Branch ( "muon_tuneP_Phi",     &muon_tuneP_Phi_);
  outputTree_ -> Branch ( "muon_tuneP_MuonBestTrackType",
                                                &muon_tuneP_MuonBestTrackType_);

  outputTree_ -> Branch ( "muon_trkIso",           &muon_trkIso_);
  outputTree_ -> Branch ( "muon_normChi2",         &muon_normChi2_);
  outputTree_ -> Branch ( "muon_chi2LocalPosition",&muon_chi2LocalPosition_);
  outputTree_ -> Branch ( "muon_kinkFinder",       &muon_kinkFinder_);
  outputTree_ -> Branch ( "muon_segmentCompatability",
                                                   &muon_segmentCompatability_);

  outputTree_ -> Branch ( "muon_dtSeg_n",          &muon_dtSeg_n_);
  outputTree_ -> Branch ( "muon_dtSeg_t0timing",   &muon_dtSeg_t0timing_);
  outputTree_ -> Branch ( "muon_dtSeg_found",      &muon_dtSeg_found_);
  outputTree_ -> Branch ( "muon_dtSeg_globX",      &muon_dtSeg_globX_);
  outputTree_ -> Branch ( "muon_dtSeg_globY",      &muon_dtSeg_globY_);
  outputTree_ -> Branch ( "muon_dtSeg_globZ",      &muon_dtSeg_globZ_);
  outputTree_ -> Branch ( "muon_dtSeg_eta",        &muon_dtSeg_eta_);
  outputTree_ -> Branch ( "muon_dtSeg_phi",        &muon_dtSeg_phi_);
  outputTree_ -> Branch( "muon_dtSeg_Station_",    &muon_dtSeg_Station_);
  outputTree_ -> Branch( "muon_dtSeg_Sector_",      &muon_dtSeg_Sector_);

  outputTree_ -> Branch ( "simHit_n",          &simHit_n_);
  outputTree_ -> Branch ( "simHit_globX",      &simHit_globX_);
  outputTree_ -> Branch ( "simHit_globY",      &simHit_globY_);
  outputTree_ -> Branch ( "simHit_globZ",      &simHit_globZ_);
  outputTree_ -> Branch ( "simHit_tof",        &simHit_tof_);

  outputTree_ -> Branch ( "muon_avgEtaFromDTseg",       &muon_avgEtaFromDTseg_);
  outputTree_ -> Branch ( "muon_avgPhiFromDTseg",       &muon_avgPhiFromDTseg_);

  outputTree_ -> Branch( "muon_rPhiSeg_correlationFactor",        &muon_rPhiSeg_correlationFactor_);
  outputTree_ -> Branch( "muon_rZSeg_correlationFactor",          &muon_rZSeg_correlationFactor_);
  
  outputTree_ -> Branch ( "muon_comb_ndof",             &muon_comb_ndof_);
  outputTree_ -> Branch ( "muon_comb_timeAtIpInOut",    &muon_comb_timeAtIpInOut_);
  outputTree_ -> Branch ( "muon_comb_timeAtIpInOutErr", &muon_comb_timeAtIpInOutErr_);
  outputTree_ -> Branch ( "muon_comb_timeAtIpOutIn",    &muon_comb_timeAtIpOutIn_);
  outputTree_ -> Branch ( "muon_comb_timeAtIpOutInErr", &muon_comb_timeAtIpOutInErr_);
  outputTree_ -> Branch ( "muon_comb_invBeta",          &muon_comb_invBeta_);
  outputTree_ -> Branch ( "muon_comb_freeInvBeta",      &muon_comb_freeInvBeta_);

  outputTree_ -> Branch( "track_n",          &track_n_);
  outputTree_ -> Branch( "track_vx",         &track_vx_);
  outputTree_ -> Branch( "track_vy",         &track_vy_);
  outputTree_ -> Branch( "track_vz",         &track_vz_);
  outputTree_ -> Branch( "track_phi",        &track_phi_);
  outputTree_ -> Branch( "track_eta",        &track_eta_);
  outputTree_ -> Branch( "track_pt",         &track_pt_);
  outputTree_ -> Branch( "track_ptErr",      &track_ptErr_);
//
}

// ------------ method called once each job just after ending the event loop  ------------
void EarthAsDMAnalyzer::endJob() { }

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void EarthAsDMAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setComment("Analyzer for cosmics searches");
  desc.addUntracked("verbosityLevel", 15)
  ->setComment("Higher the integer more verbose");
  desc.addUntracked("isData", 0)
  ->setComment("0 means MC, 1 means data");
  desc.addUntracked("hasSim", 1)
  ->setComment("1 means SimHits exists, 0 means not");
//  desc.add("muonCollection", edm::InputTag("splitMuons")) //muons1Leg
//  desc.add("muonCollection", edm::InputTag("lhcSTAMuons"))
  desc.add("muonCollection", edm::InputTag("splitMuons"))
  ->setComment("Muon collection");
  desc.add("muonTimeCollection", edm::InputTag("splitMuons", "dt"))
  ->setComment("Input collection for combined muon timing information");
  desc.add("PSimHitContainer", edm::InputTag("g4SimHits","MuonDTHits"))
  ->setComment("Input collection for g4SimHits Hit information");
  desc.add("TriggerResults", edm::InputTag("TriggerResults","","HLT"))
  ->setComment("HLTrigger results");
  desc.add("trackCollection", edm::InputTag("tevMuons:default:RECO"))
  ->setComment("TeV muon collection");

  descriptions.add("EarthAsDMAnalyzer",desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(EarthAsDMAnalyzer);


// using namespace std;

  // Function to calculate the mean of a vector
float mean(const std::vector<float>& data) {
    float sum = 0.0;
    for (size_t i = 0; i < data.size(); ++i) {
        sum += data[i];
    }
    return sum / data.size();
}

// Function to calculate the R-squared value
float calculateRSquared(const std::vector<float>& x, const std::vector<float>& y) {
    if (x.size() != y.size() || x.size() == 0) {
        return 9999; // Handle invalid input
    }

    float xMean = mean(x);
    float yMean = mean(y);
    float ssr = 0.0; // Sum of squared residuals
    float sst_x = 0.0; // Total sum of squares for x
    float sst_y = 0.0; // Total sum of squares for y

    for (size_t i = 0; i < x.size(); ++i) {
        float xDeviation = x[i] - xMean;
        float yDeviation = y[i] - yMean;
        ssr += xDeviation * yDeviation;
        sst_x += xDeviation * xDeviation;
        sst_y += yDeviation * yDeviation;
    }

    if (sst_x == 0.0 || sst_y == 0.0) {
        return 0.0; // Avoid division by zero, and R-squared is 0 in this case
    }

    float rSquared = (ssr * ssr) / (sst_x * sst_y);
    return rSquared;
}

// Function to calculate Pearson's correlation coefficient
float pearsonCorrelation(const std::vector<float>& x, const std::vector<float>& y) {
    if (x.size() != y.size()) {
        std::cerr << "Input vectors must have the same size." << std::endl;
        return 9999; // Return an error value
    }

    float sumXY = 0.0;
    float sumX2 = 0.0;
    float sumY2 = 0.0;

    float meanX = mean(x);
    float meanY = mean(y);

    for (size_t i = 0; i < x.size(); i++) {
        float deviationX = x[i] - meanX;
        float deviationY = y[i] - meanY;
        sumXY += deviationX * deviationY;
        sumX2 += deviationX * deviationX;
        sumY2 += deviationY * deviationY;
    }

    return sumXY / (sqrt(sumX2) * sqrt(sumY2));
}

//  for (unsigned int c=0; c<cscSegments->size(); c++) {
//    CSCSegmentRef segRef  = CSCSegmentRef( cscSegments, c );
//    const GeomDet* cscDet = cscGeom->idToDet(SegRef->geographicalId());
//    GlobalPoint point = cscDet->toGlobal(SegRef->localPosition());
//    cout << " X: " << point.x() << endl;
//  }

//  // loop on the dt segments to get coordinates out
//dtSeg_n = 0;
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
//    cout << " Segment found for dtSeg_n index " << dtSeg_n << " and muon index " << segmentFoundMuonIndex << " with Y =" << globalPoint.y() << endl;
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


/*
 
 TH1F *GenPt, *RecoPt, *GenPhi, *RecoPhi;
 TH1F *RelDiffTrackPtAndTruthPt;
 TH1F *RelDiffTrackPhiAndTruthPhi;
 TH2F *GenPtVsRecoPt;
 TH2F *GenPhiVsRecoPhi;
 
 TH1F *GenPt_TeVMuonsMatch, *RecoPt_TeVMuonsMatch, *GenPhi_TeVMuonsMatch, *RecoPhi_TeVMuonsMatch;
 TH1F *RelDiffTrackPtAndTruthPt_TeVMuonsMatch;
 TH1F *RelDiffTrackPhiAndTruthPhi_TeVMuonsMatch;
 TH2F *GenPtVsRecoPt_TeVMuonsMatch;
 TH2F *GenPhiVsRecoPhi_TeVMuonsMatch;

 // Create subdirectory for histograms
 TFileDirectory subDir = fs->mkdir( "Histograms" );
 // Reco track collection: splitMuons
 float upperPt = 3000.; // 3000.
 GenPt                      = subDir.make<TH1F>( "GenPt" , "GenPt", 100,0.,upperPt);
 RecoPt                     = subDir.make<TH1F>( "RecoPt" , "RecoPt", 100,0.,upperPt);
 GenPhi                     = subDir.make<TH1F>( "GenPhi" , "GenPhi", 50, -3.14, 3.14);
 RecoPhi                    = subDir.make<TH1F>( "RecoPhi" , "RecoPhi", 50, -3.14, 3.14);
 RelDiffTrackPtAndTruthPt   = subDir.make<TH1F>( "RelDiffTrackPtAndTruthPt" , "RelDiffTrackPtAndTruthPt", 200, -1., 1.);
 RelDiffTrackPhiAndTruthPhi = subDir.make<TH1F>( "RelDiffTrackPhiAndTruthPhi" , "RelDiffTrackPhiAndTruthPhi", 200, -2., 2.);
 GenPtVsRecoPt              = subDir.make<TH2F>( "GenPtVsRecoPt" , "GenPtVsRecoPt", 100,0.,upperPt,100,0.,upperPt );
 GenPhiVsRecoPhi            = subDir.make<TH2F>( "GenPhiVsRecoPhi" , "GenPhiVsRecoPhi", 50, -3.14, 3.14,50, -3.14, 3.14);
 // Reco track collection: tevMuons
 upperPt = 6000.;
 GenPt_TeVMuonsMatch                      = subDir.make<TH1F>( "GenPt_TeVMuonsMatch" , "GenPt_TeVMuonsMatch", 100,0.,upperPt);
 RecoPt_TeVMuonsMatch                     = subDir.make<TH1F>( "RecoPt_TeVMuonsMatch" , "RecoPt_TeVMuonsMatch", 100,0.,upperPt);
 GenPhi_TeVMuonsMatch                     = subDir.make<TH1F>( "GenPhi_TeVMuonsMatch" , "GenPhi_TeVMuonsMatch", 50, -3.14, 3.14);
 RecoPhi_TeVMuonsMatch                    = subDir.make<TH1F>( "RecoPhi_TeVMuonsMatch" , "RecoPhi_TeVMuonsMatch", 50, -3.14, 3.14);
 RelDiffTrackPtAndTruthPt_TeVMuonsMatch   = subDir.make<TH1F>( "RelDiffTrackPtAndTruthPt_TeVMuonsMatch" , "RelDiffTrackPtAndTruthPt_TeVMuonsMatch", 200, -1., 1.);
 RelDiffTrackPhiAndTruthPhi_TeVMuonsMatch = subDir.make<TH1F>( "RelDiffTrackPhiAndTruthPhi_TeVMuonsMatch" , "RelDiffTrackPhiAndTruthPhi_TeVMuonsMatch", 200, -2., 2.);
 GenPtVsRecoPt_TeVMuonsMatch              = subDir.make<TH2F>( "GenPtVsRecoPt_TeVMuonsMatch" , "GenPtVsRecoPt_TeVMuonsMatch", 100,0.,upperPt,100,0.,upperPt );
 GenPhiVsRecoPhi_TeVMuonsMatch            = subDir.make<TH2F>( "GenPhiVsRecoPhi_TeVMuonsMatch" , "GenPhiVsRecoPhi_TeVMuonsMatch", 50, -3.14, 3.14,50, -3.14, 3.14);
 */
