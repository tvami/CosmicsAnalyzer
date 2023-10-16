import argparse
import glob
from ROOT import TChain, TMath, TH1F, TCanvas, gStyle, TLegend, TPaveText, TF1, TF1Convolution
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
        '100GeV_muons_tevMuons'                : TChain(treename),
        '500GeV_muons_tevMuons'                : TChain(treename),
        '1000GeV_muons_tevMuons'               : TChain(treename),
        '3-4TeV_muons_tevMuons'                : TChain(treename),
    }
    
    ptHistDict = {
        'gen_100GeV' : TH1F("genPt_100GeV", "", 100, 0, 5000),
        'gen_500GeV' : TH1F("genPt_500GeV", "", 100, 0, 5000),
        'gen_1000GeV' : TH1F("genPt_1000GeV", "", 100, 0, 5000),
        'gen_3-4TeV' : TH1F("genPt_3-4TeV", "", 100, 0, 5000),
        '100GeV' : TH1F("Pt_100GeV", "", 100, 0, 5000),
        '500GeV' : TH1F("Pt_500GeV", "", 100, 0, 5000),
        '1000GeV' : TH1F("Pt_1000GeV", "", 100, 0, 5000),
        '3-4TeV' : TH1F("Pt_3-4TGeV", "", 100, 0, 5000),
    }
    
    relHistDict = {
        '100GeV' : TH1F("relDiffRecoGenPt_100GeV", "",  100, -1., 1.),
        '500GeV' : TH1F("relDiffRecoGenPt_500GeV", "",  100, -1., 1.),
        '1000GeV' : TH1F("relDiffRecoGenPt_1000GeV", "",  100, -1., 1.),
        '3-4TeV' : TH1F("relDiffRecoGenPt_3-4TeV", "",  100, -1., 1.),
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
        commonFunctions.vprint( verbose, 3, "loaded tree nevents = {}".format(tc.GetEntries()) )
        
    # tree loop
    for key, tc in treeDict.items():
        commonFunctions.vprint( verbose, 3, "!!!!! ===== {} ===== ".format(key))
        
        # Gen - Reco muon collection match
        muon_match_index = \
                commonFunctions.recoMuonGenMatch(tc, verbose, method='mindR', drawOpt=True, outdir=outdir, key=key+'_sig')
        
        energy_key = key.split("_")[0]
        commonFunctions.vprint( verbose, 3, "energy:~~{}~~".format(energy_key))
                            
        # Fill Reco-Gen pt rel diff
        for hkey in ptHistDict:
            if energy_key != '' and hkey.find(energy_key)>=0:
                for ind, evt in enumerate(tc):
                    if muon_match_index[ind] != -1: 
                        if hkey.find('gen')>=0: 
                            commonFunctions.vprint( verbose, 4, "Filling in gen pt {} hist!!".format(energy_key) )
                            ptHistDict[hkey].Fill(evt.gen_pt[1])
                        else: 
                            commonFunctions.vprint( verbose, 4, "Filling in matched pt {} hist!!".format(energy_key) )
                            ptHistDict[hkey].Fill(evt.muon_fromGenTrack_Pt[muon_match_index[ind]])
                            relDiff = (evt.muon_fromGenTrack_Pt[muon_match_index[ind]] - evt.gen_pt[1]) / evt.gen_pt[1]
                            relHistDict[hkey].Fill(relDiff)
            
    commonFunctions.drawNormHistDict_CompN(ptHistDict, "p_{T}", "pt_compareE_match_sig", outdir, logy=True, maxy=1e1, compN=4)
    commonFunctions.drawNormHistDict(relHistDict, "(Reco p_{T} #minus Gen p_{T}) / Gen p_{T}",\
                                                          "resPt_compareE_match_sig", outdir, maxy=0.4)

    # Fit to gaussian
    # fit_mean_dict = {}
    histDict = relHistDict
    for key, h in histDict.items():
        hist = h.Clone()
        xmin, xmax = -0.25, 0.2
        f1 = TF1("f1", "gaus", xmin, xmax)
        # f_conv = TF1Convolution("gaus", "gaus", -1, 1)
        # f_conv.SetRange(-1, 1)
        # f_conv.SetNofPointsFFT(1000)
        # f1 = TF1("f1", f_conv, -0.5, 0.5, f_conv.GetNpar())
        # f1.SetParLimits(1, -0.15, 0.15)  # mean
        # f1.SetParLimits(2, 0., 0.12)  # width
        
        ftr_p = hist.Fit("f1", "SR")  # option S enables retriving the full fit result (TFitResultPtr)
        ftr = ftr_p.Get()  # TFitResult
        mean = ftr.Parameter(1)
        width = ftr.Parameter(2)
        mean_err = ftr.ParError(1)
        width_err = ftr.ParError(2)
        chi2 = ftr.Chi2()
        commonFunctions.vprint( verbose, 3, "Gaussian fit mean = {}, width = {}".format(mean, width) )
        # fit_mean_dict[key] = [mean, width]
        fitCan = TCanvas(key, "", 800, 600)
        fitCan.cd()
        hist.SetTitle("Gaussian Fit: tevMuons #theta 91-180, {}".format(key.split("_")[-1]))
        hist.GetXaxis().SetTitle("(Reco p_{T} #minus Gen p_{T}) / Gen p_{T}")
        hist.GetYaxis().SetTitle("A. U.")
        hist.SetMarkerSize(1)
        hist.SetMinimum(0)
        hist.Draw("E1 func")
        pt = TPaveText(0.2, 0.7, 0.5, 0.94, "NDC")
        pt.SetFillColor(0)
        pt.SetAllWith("", "align",11)
        pt.AddText("tevMuons #theta 91-180, {}".format(key.split("_")[-1]))
        pt.AddText("Gaussian Fit ({} < x < {}):".format(xmin, xmax))
        pt.AddText("mean = {:.3f} #pm {:.3f}".format(mean, mean_err))
        pt.AddText("width = {:.3f} #pm {:.3f}".format(width, width_err))
        pt.AddText("chi2 = {:.1f}".format(chi2))
        # pt.Paint("NDC")
        pt.Draw("same")
        fitCan.SaveAs(outdir + "/" + "resPt_fit_{}_sig.png".format(key))


    
    
if __name__ == "__main__":
    main()