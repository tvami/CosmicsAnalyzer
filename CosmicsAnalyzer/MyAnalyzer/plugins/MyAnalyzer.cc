// -*- C++ -*-
//
// Package:    CosmicsAnalyzer/MyAnalyzer
// Class:      MyAnalyzer
//
/**\class MyAnalyzer MyAnalyzer.cc CosmicsAnalyzer/MyAnalyzer/plugins/MyAnalyzer.cc

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
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TFile.h"
#include "TTree.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;

class MyAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit MyAnalyzer(const edm::ParameterSet&);
  ~MyAnalyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<std::vector<reco::Muon>> muonToken_;
  TFile* outputFile_;
  TTree* outputTree_;
  std::vector<double> muonPhis_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MyAnalyzer::MyAnalyzer(const edm::ParameterSet& iConfig) :
 muonToken_(consumes<std::vector<reco::Muon>>(iConfig.getParameter<edm::InputTag>("muonCollection"))) {}

MyAnalyzer::~MyAnalyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void MyAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& setup) {
  edm::Handle<std::vector<reco::Muon>> muons;
  event.getByToken(muonToken_, muons);

  for (const auto& muon : *muons) {
    double phi = muon.phi();
    muonPhis_.push_back(phi);
  }
}

// ------------ method called once each job just before starting event loop  ------------
void MyAnalyzer::beginJob() {
  outputFile_ = new TFile("muon_phi_ntuple.root", "RECREATE");
  outputTree_ = new TTree("MuonPhiTree", "Muon Phi Distribution");
  outputTree_->Branch("muonPhi", &muonPhis_);
}

// ------------ method called once each job just after ending the event loop  ------------
void MyAnalyzer::endJob() {
  outputTree_->Fill();
  outputFile_->Write();
  outputFile_->Close();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void MyAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MyAnalyzer);
