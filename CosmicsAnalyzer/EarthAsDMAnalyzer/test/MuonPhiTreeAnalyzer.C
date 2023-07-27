#define MuonPhiTreeAnalyzer_cxx
#include "MuonPhiTreeAnalyzer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void MuonPhiTreeAnalyzer::Loop()
{
   TH1::SetDefaultSumw2(kTRUE);
   TH2::SetDefaultSumw2(kTRUE);
   TH2D *histogram = new TH2D("histogram", "Muon #Phi vs Muon p_{T};Muon #Phi;Muon p_{T}", 100, -3.5, 3.5, 200, 0, 1000);

   if (fChain == 0) return;

   Long64_t nEntries = fChain->GetEntriesFast();
  
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nEntries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
       fChain->GetEntry(jentry);
       for (size_t iMuon = 0; iMuon < muonPhi->size(); ++iMuon) {
         histogram->Fill(muonPhi->at(iMuon), muonPt->at(iMuon));
      
     }
          
     TFile outputFile("outputFileName.root", "RECREATE");
     histogram->Write();
     outputFile.Close();
     delete histogram;
   }
}

void MuonPhiTreeAnalyzer() {
  class MuonPhiTreeAnalyzer t;
  t.Loop();
}
