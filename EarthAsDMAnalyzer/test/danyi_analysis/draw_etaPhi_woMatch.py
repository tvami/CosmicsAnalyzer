import argparse
import glob
from ROOT import TChain, TH1F, TCanvas, gStyle, TLegend
import tdrStyle

def drawNormHistDict(histDict, xtitle, filename, outdir):
    gStyle.SetOptStat(0)
    tdrStyle.SetTDRStyle()
    
    can = TCanvas("can_{}".format(filename), "", 800, 600)
    
    # histDict = {k: v for k, v in sorted(histDict.items(), key = lambda item: item[1].GetMaximum()/item[1].Integral(), reverse=True)}
    histList = list(histDict.items())
    
    integrals = []
    for i in range(len(histList)):
        hist = histList[i][1]
        if i == 0:
            hist.GetXaxis().SetTitle(xtitle)
            hist.GetYaxis().SetTitle("A. U.")
        hist.SetLineColor(tdrStyle.colors['color_comp{}'.format(i+1)])
        integrals.append(hist.Integral())
        hist.Scale(1./hist.Integral())
        hist.SetLineWidth(2)
        hist.SetMaximum(0.15)

    # draw
    can.cd()
    # can.SetLogy()
    for key, hist in histDict.items():
        if i == 0:
            # hist.SetMinimum(1e-4)
            hist.Draw("hist")
        else:
            hist.Draw("hist same")
            
    leg = TLegend(0.25, 0.7, 0.8, 0.94)
    for i, (k, hist) in enumerate(histDict.items()):
        if k.find("gen")>=0:
            leg.AddEntry(hist, "gen : {:.0f}".format(integrals[i]))
    for i, (k, hist) in enumerate(histDict.items()):
        if k.find("muons")>=0:
            leg.AddEntry(hist, "muons: {:.0f}".format(integrals[i]))
    for i, (k, hist) in enumerate(histDict.items()):
        if k.find("DT")>=0:
            leg.AddEntry(hist, "muons DT avg: {:.0f}".format(integrals[i]))
    for i, (k, hist) in enumerate(histDict.items()):
        if not k.find("gen")>=0 and not k.find("muons")>=0 and not k.find("DT")>=0:
            leg.AddEntry(hist, k + ": {:.0f}".format(integrals[i]))

    leg.Draw("same")
    can.SaveAs(outdir + "/" + filename + ".png")

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
        'muons_standAloneMuons'         : TChain(treename),
        'muons_tevMuons'                : TChain(treename),
        'muons_cosmicMuons'             : TChain(treename),
        # 'splitMuons_tevMuons'           : TChain(treename),
        # 'standAloneMuons' : TChain(treename),
        # 'cosmicMuons'     : TChain(treename),
        # 'cosmicMuons1Leg' : TChain(treename)
    }
    
    etaHistDict = {
        'gen' : TH1F("genEta", "", 50, -2., 2.),
        'muons' : TH1F("muonsEta", "", 50, -2., 2.),
        'DT' : TH1F("DTEta", "", 50, -2., 2.),
        'tevMuons' : TH1F("tevMuonsEta", "", 50, -2., 2.),
        'ctfWithMaterialTracksP5' : TH1F("ctfWithMaterialTracksP5Eta", "", 50, -2., 2.),
        'standAloneMuons' : TH1F("standAloneMuonsEta", "", 50, -2., 2.),
        'cosmicMuons' : TH1F("cosmicMuonsEta", "", 50, -2., 2.),
    }
    
    phiHistDict = {
        'gen' : TH1F("genPhi", "", 50, -3.14, 3.14),
        'muons' : TH1F("muonsPhi", "", 50, -3.14, 3.14),
        'DT' : TH1F("DTPhi", "", 50, -3.14, 3.14),
        'tevMuons' : TH1F("tevMuonsPhi", "", 50, -3.14, 3.14),
        'ctfWithMaterialTracksP5' : TH1F("ctfWithMaterialTracksP5Phi", "", 50, -3.14, 3.14),
        'standAloneMuons' : TH1F("standAloneMuonsPhi", "", 50, -3.14, 3.14),
        'cosmicMuons' : TH1F("cosmicMuonsPhi", "", 50, -3.14, 3.14),
    }
    
    # Load the trees
    for key, tc in treeDict.items():
        # # test local ntuples
        for file_match in glob.glob(indir + '/*{}*.root'.format(key)):
            tc.Add(file_match)
        # crab ntuples
        # for file_match in glob.glob(indir + '/*1000GeV_{}*/*/*/*.root'.format(key)):
        #     tc.Add(file_match)
        if verbose > 2: print("loaded tree nevents =", tc.GetEntries())
        
    # tree loop
    for key, tc in treeDict.items():
        if verbose > 2: print("!!!!! ===== {} ===== ".format(key))
        
        muon_key = key.split("_")[0]
        track_key = key.split("_")[1]
        if verbose > 2: print("muon_key:~~{}~~".format(muon_key))
        if verbose > 2: print("track_key:~~{}~~".format(track_key))
        
        # Fill gen pT only for once
        if etaHistDict['gen'].GetEntries() == 0. or phiHistDict['gen'].GetEntries() == 0.:
            if verbose > 2: print("  >> Filling in gen hists")
            for evt in tc:
                etaHistDict['gen'].Fill(evt.gen_eta[1])
                phiHistDict['gen'].Fill(evt.gen_phi[1])
        
        # Fill muons only for once
        if muon_key.find("muons")>=0 and\
           etaHistDict['muons'].GetEntries() == 0. or phiHistDict['muons'].GetEntries() == 0.:
            if verbose > 2: print("  >> Filling in muons hists")
            # for evt in tc:
            #     etaHistDict[muon_key].Fill(evt.muon_eta)
            #     phiHistDict[muon_key].Fill(evt.muon_phi)
            tc.Draw("muon_eta >> {}Eta".format(muon_key), "", "hist")
            tc.Draw("muon_phi >> {}Phi".format(muon_key), "", "hist")
        
        if etaHistDict['DT'].GetEntries() == 0. or phiHistDict['DT'].GetEntries() == 0.:
            tc.Draw("muon_avgEtaFromDTseg >> {}Eta".format("DT"), "", "hist")
            tc.Draw("muon_avgPhiFromDTseg >> {}Phi".format("DT"), "", "hist")
        
        # Fill tracks
        for hkey in etaHistDict:
            if track_key.find(hkey)>=0 and etaHistDict[hkey].GetEntries() == 0.:
                if verbose > 2: print("  >> Filling in {} hists".format(track_key))
                # for evt in tc:
                #     etaHistDict[track_key].Fill(evt.track_eta)
                #     phiHistDict[track_key].Fill(evt.track_phi)
                tc.Draw("track_eta >> {}Eta".format(track_key), "", "hist")
                tc.Draw("track_phi >> {}Phi".format(track_key), "", "hist")
        
    
    drawNormHistDict(etaHistDict, "#eta", "eta_compare_sig", outdir)
    drawNormHistDict(phiHistDict, "#phi", "phi_compare_sig", outdir)
    
    

    
if __name__ == "__main__":
    main()