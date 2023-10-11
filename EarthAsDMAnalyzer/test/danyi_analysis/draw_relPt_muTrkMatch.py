import argparse
import glob
from ROOT import TChain, TMath, TH1F, TCanvas, gStyle, TLegend, TPaveText
import tdrStyle
from draw_pt_muonColl import dR
from draw_pt_muTrkMatch import recoMuonTrackMatch

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
        # 'muons_ctfWithMaterialTracksP5' : TChain(treename),
        # 'muons_standAloneMuons'         : TChain(treename),
        '100GeV_muons_tevMuons'                : TChain(treename),
        '500GeV_muons_tevMuons'                : TChain(treename),
        '1000GeV_muons_tevMuons'               : TChain(treename),
        '3-4TeV_muons_tevMuons'                : TChain(treename),
        # 'splitMuons_tevMuons'           : TChain(treename),
        # 'standAloneMuons' : TChain(treename),
        # 'cosmicMuons'     : TChain(treename),
        # 'cosmicMuons1Leg' : TChain(treename)
    }
    
    # histDict = {
    #     'genPt' : TH1F("genPt", "", 100, 0, 2000),
    #     'recoMuonPt_muons' : TH1F("recoMuonPt_muons", "", 100, 0, 2000),
    #     # 'recoMuonPt_splitMuons' : TH1F("recoMuonPt_splitMuons", "", 100, 0, 2000),
    #     'recoTrackPt_tevMuons' : TH1F("recoMuonPt_tevMuons", "", 100, 0, 2000),
    #     'recoTrackPt_ctfWithMaterialTracksP5' : TH1F("recoMuonPt_ctfWithMaterialTracksP5", "", 100, 0, 2000),
    #     'recoTrackPt_standAloneMuons' : TH1F("recoMuonPt_standAloneMuons", "", 100, 0, 2000),
    # }
    
    relHistDict = {
        # 'relDiffRecoGenPt_muons' : TH1F("relDiffRecoGenPt_muons", "", 50, -1., 1.),
        # 'relDiffRecoGenPt_splitMuons' : TH1F("relDiffRecoGenPt_splitMuons", "",  50, -1., 1.)
        'relDiffRecoGenPt_tevMuons_100GeV' : TH1F("relDiffRecoGenPt_tevMuons_100GeV", "",  50, -1., 1.),
        'relDiffRecoGenPt_tevMuons_500GeV' : TH1F("relDiffRecoGenPt_tevMuons_500GeV", "",  50, -1., 1.),
        'relDiffRecoGenPt_tevMuons_1000GeV' : TH1F("relDiffRecoGenPt_tevMuons_1000GeV", "",  50, -1., 1.),
        'relDiffRecoGenPt_tevMuons_3-4TeV' : TH1F("relDiffRecoGenPt_tevMuons_1000GeV", "",  50, -1., 1.),
    }
    
    # Load the trees
    for key, tc in treeDict.items():
        # # test local ntuples
        # for file_match in glob.glob(indir + '/*{}*.root'.format(key)):
        #     tc.Add(file_match)
        # crab ntuples
        for file_match in glob.glob(indir + '/*{}*/*/*/*.root'.format(key)):
            tc.Add(file_match)
        if key.find('3-4TeV')>=0:
            tc.Add("/home/users/dazhang/works/CMSSW_12_6_5/src/CosmicsAnalyzer/EarthAsDMAnalyzer/test/ntuples/ntuple_MC_RR-91to180Theta-3000to4000GeV_74_muons_tevMuons.root")
        if verbose > 2: print("loaded tree nevents =", tc.GetEntries())
        
    # tree loop
    for key, tc in treeDict.items():
        if verbose > 2: print("!!!!! ===== {} ===== ".format(key))
        
        # Gen - Reco muon collection match
        muon_match_index, track_match_index = recoMuonTrackMatch(tc, verbose, outdir, key)
        
        if verbose > 2: print(len(muon_match_index), tc.GetEntries())
        
        energy_key = key.split("_")[0]
        if verbose > 2: print("energy:~~{}~~".format(energy_key))
                            
        # Fill Reco-Gen track rel diff
        for hkey in relHistDict:
            if energy_key != '' and hkey.find(energy_key)>=0:
                for ind, evt in enumerate(tc):
                    if track_match_index[ind] != -1: 
                        relDiff = (evt.track_pt[track_match_index[ind]] - evt.gen_pt[1]) / evt.gen_pt[1]
                        relHistDict[hkey].Fill(relDiff)
            
    gStyle.SetOptStat(0)
    tdrStyle.SetTDRStyle()
    
    # pT
    can = TCanvas("can", "", 800, 600)
    
    histDict = {k: v for k, v in sorted(relHistDict.items(), key = lambda item: item[1].GetMaximum()/item[1].Integral(), reverse=True)}
    histList = list(histDict.items())
    
    integrals = []
    for i in range(len(histList)):
        hist = histList[i][1]
        if i == 0:
            hist.GetXaxis().SetTitle("(Reco p_{T} #minus Gen p_{T}) / Gen p_{T}")
            hist.GetYaxis().SetTitle("A. U.")
        hist.SetLineColor(tdrStyle.colors['color_comp{}'.format(i+1)])
        integrals.append(hist.Integral())
        hist.Scale(1./hist.Integral())
        hist.SetLineWidth(2)
        
    # Draw
    can.cd()
    # can.SetLogy()
    for key, hist in histDict.items():
        if i == 0:
            # hist.SetMinimum(1e-4)
            hist.Draw("hist")
        else:
            hist.Draw("hist same")
    # legend
    leg = TLegend(0.6, 0.75, 0.94, 0.94)
    for i, (k, hist) in enumerate(histDict.items()):
        leg.AddEntry(hist, k.split("_")[-1] + ": " + str(integrals[i]))
        # leg.AddEntry(hist, k.split("_")[-1] + ": {:.1f} #pm {:.1f}".format(fit_mean_dict[k][0], fit_mean_dict[k][1]))
    leg.Draw("same")
    can.SaveAs(outdir + "/" + "resPt_muTrkMatch_sig.png")
    
    
    # Fit to gaussian
    # fit_mean_dict = {}
    for key, h in histDict.items():
        hist = h.Clone()
        ftr_p = hist.Fit("gaus", "S")  # option S enables retriving the full fit result (TFitResultPtr)
        ftr = ftr_p.Get()  # TFitResult
        mean = ftr.Parameter(0)
        width = ftr.Parameter(2)
        if verbose > 2: print("Gaussian fit mean =", mean, "width =", width)
        # fit_mean_dict[key] = [mean, width]
        fitCan = TCanvas(key, "", 800, 600)
        fitCan.cd()
        hist.SetTitle("Gaussian Fit: tevMuons #theta 91-180, {}".format(key.split("_")[-1]))
        hist.GetXaxis().SetTitle("(Reco p_{T} #minus Gen p_{T}) / Gen p_{T}")
        hist.GetYaxis().SetTitle("A. U.")
        hist.SetMarkerSize(1)
        hist.Draw("E1 func")
        pt = TPaveText(0.2, 0.7, 0.5, 0.94, "NDC")
        pt.SetFillColor(0)
        pt.SetAllWith("", "align",11)
        pt.AddText("tevMuons #theta 91-180, {}".format(key.split("_")[-1]))
        pt.AddText("Gaussian Fit:")
        pt.AddText("mean = {:.3f}".format(mean))
        pt.AddText("width = {:.3f}".format(width))
        # pt.Paint("NDC")
        pt.Draw("same")
        fitCan.SaveAs(outdir + "/" + "resPt_muTrkMatch_fit_{}_sig.png".format(key))
    

    
    
if __name__ == "__main__":
    main()