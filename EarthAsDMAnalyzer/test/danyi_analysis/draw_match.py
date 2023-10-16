import argparse
import glob
from ROOT import TChain, TH1F, TCanvas, gStyle, TLegend
import tdrStyle
import commonFunctions

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='')
    # parser.add_argument('-d','--indirs', nargs='+', action='store', dest='indirs', help='Director(y/ies) of input files, absolute path(s)')
    parser.add_argument('-i','--indir', action='store', dest='indir', help='Directory of the input ntuples')
    parser.add_argument('-o','--outdir', action='store', dest='outdir', help='Directory of output files, absolute path')
    parser.add_argument('-t','--tree', action='store', dest='treename', help='Name of the tree in the ntuple')
    parser.add_argument('-v','--verbose', action='store', dest='verbose', help='Set for extra printout')
    args = parser.parse_args()
    
    if not args.indir:
        parser.error('provide input directory of ntuples, -i or --indir')
    
    if not args.outdir:
        parser.error('provide output directory, -o or --outdir')
        
    if not args.treename:
        parser.error('provide tree name, -t or --tree')
    
    if not args.verbose:
        parser.error('provide verbose level, -v or --verbose')
    
    treename = args.treename
    indir = args.indir
    outdir = args.outdir
    verbose = int(args.verbose)
    
    treeDict = {
        'muons_ctfWithMaterialTracksP5' : TChain(treename),
        # 'muons_standAloneMuons'         : TChain(treename),
        'muons_tevMuons'                : TChain(treename),
        'muons_cosmicMuons'             : TChain(treename),
        # 'splitMuons_tevMuons'           : TChain(treename),
        # 'standAloneMuons' : TChain(treename),
        # 'cosmicMuons'     : TChain(treename),
        # 'cosmicMuons1Leg' : TChain(treename)
    }
    
    etaHistDict = {
        'gen' : TH1F("genEta", "", 50, -2., 2.),
        # 'muons' : TH1F("muonsEta", "", 50, -2., 2.),
        # 'DT' : TH1F("DTEta", "", 50, -2., 2.),
        'tevMuons' : TH1F("tevMuonsEta", "", 50, -2., 2.),
        'ctfWithMaterialTracksP5' : TH1F("ctfWithMaterialTracksP5Eta", "", 50, -2., 2.),
        # 'standAloneMuons' : TH1F("standAloneMuonsEta", "", 50, -2., 2.),
        'cosmicMuons' : TH1F("cosmicMuonsEta", "", 50, -2., 2.),
    }
    
    phiHistDict = {
        'gen' : TH1F("genPhi", "", 50, -3.14, 3.14),
        # 'muons' : TH1F("muonsPhi", "", 50, -3.14, 3.14),
        # 'DT' : TH1F("DTPhi", "", 50, -3.14, 3.14),
        'tevMuons' : TH1F("tevMuonsPhi", "", 50, -3.14, 3.14),
        'ctfWithMaterialTracksP5' : TH1F("ctfWithMaterialTracksP5Phi", "", 50, -3.14, 3.14),
        # 'standAloneMuons' : TH1F("standAloneMuonsPhi", "", 50, -3.14, 3.14),
        'cosmicMuons' : TH1F("cosmicMuonsPhi", "", 50, -3.14, 3.14),
    }
    
    ptHistDict = {
        'gen' : TH1F("genPt", "", 100, 0, 5000),
        # 'muons' : TH1F("muonsPt", "", 100, 0, 5000),
        # 'DT' : TH1F("DTPt", "", 100, 0, 5000),
        'tevMuons' : TH1F("tevMuonsPt", "", 100, 0, 5000),
        'ctfWithMaterialTracksP5' : TH1F("ctfWithMaterialTracksP5Pt", "", 100, 0, 5000),
        # 'standAloneMuons' : TH1F("standAloneMuonsPt", "", 100, 0, 5000),
        'cosmicMuons' : TH1F("cosmicMuonsPt", "", 100, 0, 5000),
    }
    
    # Load the trees
    for key, tc in treeDict.items():
        # # test local ntuples
        for file_match in glob.glob(indir + '/*{}*.root'.format(key)):
            tc.Add(file_match)
        # crab ntuples
        # for file_match in glob.glob(indir + '/*1000GeV_{}*/*/*/*.root'.format(key)):
        #     tc.Add(file_match)
        commonFunctions.vprint(verbose, 3, "loaded tree nevents = {}".format(tc.GetEntries()))
        
    # tree loop
    for key, tc in treeDict.items():
        commonFunctions.vprint(verbose, 3, "!!!!! ===== {} ===== ".format(key))
        
        muon_key = key.split("_")[0]
        track_key = key.split("_")[1]
        commonFunctions.vprint(verbose, 3, "muon_key:~~{}~~".format(muon_key))
        commonFunctions.vprint(verbose, 3, "track_key:~~{}~~".format(track_key))
        
        # Reco-Gen match
        muon_match_index = commonFunctions.recoMuonGenMatch(tc, verbose, method='mindR', drawOpt=True, outdir=outdir, key=key+'_3-4TeV_sig')
        
        # Fill gen pT only for once
        if etaHistDict['gen'].GetEntries() == 0. or phiHistDict['gen'].GetEntries() == 0.:
            commonFunctions.vprint(verbose, 3, "  >> Filling in gen hists")
            for evt in tc:
                etaHistDict['gen'].Fill(evt.gen_eta[1])
                phiHistDict['gen'].Fill(evt.gen_phi[1])
                ptHistDict['gen'].Fill(evt.gen_pt[1])
        
        '''
        # Fill muons only for once
        if muon_key.find("muons")>=0:
            commonFunctions.vprint(verbose, 3, "  >> Filling in muons hists")
            if etaHistDict['muons'].GetEntries() == 0. or phiHistDict['muons'].GetEntries() == 0.:
                # tc.Draw("muon_eta >> {}Eta".format(muon_key), "", "hist")
                # tc.Draw("muon_phi >> {}Phi".format(muon_key), "", "hist")
                for ind, evt in enumerate(tc):
                    if muon_match_index[ind] != -1:
                        etaHistDict['muons'].Fill(evt.muon_fromGenTrack_Eta[muon_match_index[ind]])
                        phiHistDict['muons'].Fill(evt.muon_fromGenTrack_Phi[muon_match_index[ind]])
        
        if etaHistDict['DT'].GetEntries() == 0. or phiHistDict['DT'].GetEntries() == 0.:
            # tc.Draw("muon_avgEtaFromDTseg >> {}Eta".format("DT"), "", "hist")
            # tc.Draw("muon_avgPhiFromDTseg >> {}Phi".format("DT"), "", "hist")
            for ind, evt in enumerate(tc):
                if muon_match_index[ind] != -1:
                    etaHistDict['DT'].Fill(evt.muon_avgEtaFromDTseg[muon_match_index[ind]])
                    phiHistDict['DT'].Fill(evt.muon_avgPhiFromDTseg[muon_match_index[ind]])
        '''
        
        # Fill tracks
        for hkey in etaHistDict:
            if track_key.find(hkey)>=0 and etaHistDict[hkey].GetEntries() == 0.:
                commonFunctions.vprint(verbose, 3, "  >> Filling in {} hists".format(track_key))
                # tc.Draw("track_eta >> {}Eta".format(track_key), "", "hist")
                # tc.Draw("track_phi >> {}Phi".format(track_key), "", "hist")
                for ind, evt in enumerate(tc):
                    if muon_match_index[ind] != -1:
                        etaHistDict[track_key].Fill(evt.muon_avgEtaFromDTseg[muon_match_index[ind]])
                        phiHistDict[track_key].Fill(evt.muon_avgPhiFromDTseg[muon_match_index[ind]])
                        ptHistDict[track_key].Fill(evt.muon_fromGenTrack_Pt[muon_match_index[ind]])
        
    
    commonFunctions.drawNormHistDict(etaHistDict, "#eta", "eta_compare_match_sig", outdir)
    commonFunctions.drawNormHistDict(phiHistDict, "#phi", "phi_compare_match_sig", outdir)
    commonFunctions.drawNormHistDict(ptHistDict, "p_{T}", "pt_compare_match_sig", outdir, logy=True)
    

    
if __name__ == "__main__":
    main()