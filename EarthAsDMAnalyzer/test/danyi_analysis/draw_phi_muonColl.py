import argparse
import glob
from ROOT import TChain, TH1F, TCanvas, gStyle, TLegend
from tdrStyle import SetTDRStyle
from draw_pt_muonColl import recoMuonGenMatch


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
        'muons'      : TChain(treename),
        'splitMuons' : TChain(treename),
        # 'standAloneMuons' : TChain(treename),
        # 'cosmicMuons'     : TChain(treename),
        # 'cosmicMuons1Leg' : TChain(treename)
    }
    
    histDict = {
        'genPhi' : TH1F("genPhi", "", 50, -3.14, 3.14),
        'recoMuonPhi_muons' : TH1F("recoMuonPhi_muons", "", 50, -3.14, 3.14),
        'recoMuonPhi_splitMuons' : TH1F("recoMuonPhi_splitMuons", "", 50, -3.14, 3.14)
    }
    
    relHistDict = {
        'relDiffRecoGenPhi_muons' : TH1F("relDiffRecoGenPhi_muons", "", 50, -1., 1.),
        'relDiffRecoGenPhi_splitMuons' : TH1F("relDiffRecoGenPhi_splitMuons", "",  50, -1., 1.)
    }
    
    # Load the trees
    for key, tc in treeDict.items():
        for file_match in glob.glob(indir + '/*{}*.root'.format(key)):
            tc.Add(file_match)
        # tc.Print()
        
    # tree loop
    for key, tc in treeDict.items():
        if verbose > 2: print("!!!!! ===== {} ===== ".format(key))
        muon_match_index = recoMuonGenMatch(tc, verbose)
        # Gen - Reco muon collection match
        if verbose > 2: print(len(muon_match_index), tc.GetEntries())
        
        # Fill genPhi
        if key.find("muons")>=0:
            for evt in tc:
                histDict['genPhi'].Fill(evt.gen_phi[1])
        
        for hkey in histDict:
            # event loop
            if hkey.find(key)>=0:
                if verbose > 2: print("found", key, "in", hkey)
                for ind, evt in enumerate(tc):
                    if muon_match_index[ind] != -1: 
                        if verbose > 3: print("plotting: >> match index" , muon_match_index[ind])
                        if verbose > 3: print("plotting: >> muon phi" , evt.muon_phi)
                        histDict[hkey].Fill(evt.muon_phi[muon_match_index[ind]])
        for hkey in relHistDict:
            if hkey.find(key)>=0:
                for ind, evt in enumerate(tc):
                    if muon_match_index[ind] != -1: 
                        # relDiff = (evt.muon_phi[muon_match_index[ind]] - evt.gen_phi[1]) / evt.gen_phi[1]
                        relDiff = (evt.muon_phi[muon_match_index[ind]] + evt.gen_phi[1]) / (-evt.gen_phi[1])
                        relHistDict[hkey].Fill(relDiff)
            
    gStyle.SetOptStat(0)
    SetTDRStyle()
    
    # phi: gen, muons, splitMuons
    can = TCanvas("can", "", 800, 600)
    
    histDict = {k: v for k, v in sorted(histDict.items(), key = lambda item: item[1].GetMaximum(), reverse=True)}
    histList = list(histDict.items())

    histList[0][1].GetXaxis().SetTitle("#phi")
    histList[0][1].GetYaxis().SetTitle("A. U.")
    histList[0][1].SetLineColor(1)
    histList[0][1].Scale(1./histList[0][1].Integral())
    
    histList[1][1].SetLineColor(4)
    histList[1][1].Scale(1./histList[1][1].Integral())
    
    histList[2][1].SetLineColor(2)
    histList[2][1].Scale(1./histList[2][1].Integral())
    # draw
    can.cd()
    histList[0][1].Draw("hist")
    histList[1][1].Draw("same hist")
    histList[2][1].Draw("same hist")
    leg = TLegend(0.25, 0.8, 0.6, 0.9)
    for k, hist in histDict.items():
        leg.AddEntry(hist, k)
    leg.Draw("same")
    can.SaveAs(outdir + "/" + "phi_maxPt_sig.png")
    
    
    # rel diff reco - gen phi
    can2 = TCanvas("can2", "", 800, 600)
    
    relHistDict = {k: v for k, v in sorted(relHistDict.items(), key = lambda item: item[1].GetMaximum(), reverse=True)}
    relHistList = list(relHistDict.items())
    # print(relHistList)
    
    # relHistList[0][1].GetXaxis().SetTitle("(Reco #phi #minus Gen #phi) / Gen #phi")
    relHistList[0][1].GetXaxis().SetTitle("(Reco #phi + Gen #phi) / #minus Gen #phi")
    relHistList[0][1].GetYaxis().SetTitle("A. U.")
    relHistList[0][1].SetLineColor(1)
    relHistList[0][1].Scale(1./relHistList[0][1].Integral())
    
    relHistList[1][1].SetLineColor(4)
    relHistList[1][1].Scale(1./relHistList[1][1].Integral())
    
    can2.cd()
    relHistList[0][1].Draw("hist")
    relHistList[1][1].Draw("same hist")
    leg2 = TLegend(0.25, 0.8, 0.6, 0.9)
    for k, hist in relHistDict.items():
        leg2.AddEntry(hist, k)
    leg2.Draw("same")
    # can2.SaveAs(outdir + "/" + "relDiffRecoGenPhi_maxPt_sig.png")
    # can2.SaveAs(outdir + "/" + "absDiffRecoGenPhi_rvsGenPhi_maxPt_sig.png")
    can2.SaveAs(outdir + "/" + "relDiffRecoGenPhi_rvsGenPhi_maxPt_sig.png")
    
if __name__ == "__main__":
    main()