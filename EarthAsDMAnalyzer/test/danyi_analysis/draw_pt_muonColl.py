import argparse
import glob
from ROOT import TChain, TMath, TH1F, TCanvas, gStyle, TLegend
import tdrStyle

def dR(eta1, phi1, eta2, phi2):
    deta = eta1 - eta2
    dphi = phi1 - phi2
    while dphi > TMath.Pi():
        dphi -= 2 * TMath.Pi()
    while dphi <= -TMath.Pi():
        dphi += 2 * TMath.Pi()
    return TMath.Sqrt(deta * deta + dphi * dphi)

# match reco muon with gen muon (index 1)
def recoMuonGenMatch(tree, verbose):
    match_ind = []
    for ind, evt in enumerate(tree):
        if verbose > 2 and ind % 10000 == 0: print("Processing event", ind)
        if verbose > 3: print("---------------- Event:", ind, "----------------")
        if evt.muon_n > 0:
            genEta = evt.gen_eta[1]
            genPhi = evt.gen_phi[1]
            sum_ind = 0
            min_dR = 9999.
            min_dR_muonInd = -1
            max_pt = 0.
            max_pt_muonInd = -1
            for j in range(evt.muon_n):
                if verbose > 3: print("   >> muon:", j)
                if verbose > 3: print("   >> dtSeg_n:", evt.muon_dtSeg_n[j])
                if verbose > 3: print("   >> dtSeg_eta:", evt.muon_dtSeg_eta[sum_ind:(sum_ind+evt.muon_dtSeg_n[j])])
                if verbose > 3: print("   >> dtSeg_phi:", evt.muon_dtSeg_phi[sum_ind:(sum_ind+evt.muon_dtSeg_n[j])])
                sumEtaFromDTseg = sum(eta for eta in evt.muon_dtSeg_eta[sum_ind:(sum_ind+evt.muon_dtSeg_n[j])] if eta < 9000.)
                sumPhiFromDTseg = sum(phi for phi in evt.muon_dtSeg_phi[sum_ind:(sum_ind+evt.muon_dtSeg_n[j])] if phi < 9000.)
                sum_ind += evt.muon_dtSeg_n[j]
                muon_avgEtaFromDTseg = sumEtaFromDTseg / evt.muon_dtSeg_n[j] if evt.muon_dtSeg_n[j] > 0 else 9999.
                muon_avgPhiFromDTseg = sumPhiFromDTseg / evt.muon_dtSeg_n[j] if evt.muon_dtSeg_n[j] > 0 else 9999.
                muon_dR = dR(genEta, genPhi, muon_avgEtaFromDTseg, muon_avgPhiFromDTseg)
                if verbose > 3: print ("   >> dR:", muon_dR)
                # if muon_dR < 1 and muon_dR < min_dR:
                #     min_dR = muon_dR
                #     min_dR_muonInd = j
                if evt.muon_pt[j] > max_pt:
                    max_pt = evt.muon_pt[j]
                    max_pt_muonInd = j
            # match_ind.append(min_dR_muonInd)
            match_ind.append(max_pt_muonInd)
        else:
            match_ind.append(-1)
        # if ind > 10:
        #     break
    if verbose > 3: print(match_ind)
    return match_ind
    
def recoTrackGenMatch(tree, verbose):
    match_ind = []
    for ind, evt in enumerate(tree):
        if verbose > 2 and ind % 10000 == 0: print("Processing event", ind)
        if verbose > 3: print("---------------- Event:", ind, "----------------")
        if evt.track_n > 0:
            max_pt = 0.
            max_pt_trackInd = -1
            for j in range(evt.track_n):
                if verbose > 3: print("   >> track:", j)
                if verbose > 3: print ("   >> pT:", evt.track_pt[j])
                if evt.track_pt[j] > max_pt:
                    max_pt = evt.track_pt[j]
                    max_pt_trackInd = j
            match_ind.append(max_pt_trackInd)
        else:
            match_ind.append(-1)
        # if ind > 10:
        #     break
    if verbose > 3: print(match_ind)
    return match_ind

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
        '_standAloneMuons'              : TChain(treename),
        'splitMuons_tevMuons'           : TChain(treename),
        # 'standAloneMuons' : TChain(treename),
        # 'cosmicMuons'     : TChain(treename),
        # 'cosmicMuons1Leg' : TChain(treename)
    }
    
    histDict = {
        'genPt' : TH1F("genPt", "", 100, 0, 5000),
        'recoMuonPt_muons' : TH1F("recoMuonPt_muons", "", 100, 0, 5000),
        'recoMuonPt_splitMuons' : TH1F("recoMuonPt_splitMuons", "", 100, 0, 5000),
        'recoTrackPt_tevMuons' : TH1F("recoMuonPt_tevMuons", "", 100, 0, 5000),
        'recoTrackPt_ctfWithMaterialTracksP5' : TH1F("recoMuonPt_ctfWithMaterialTracksP5", "", 100, 0, 5000),
        'recoTrackPt_standAloneMuons' : TH1F("recoMuonPt_standAloneMuons", "", 100, 0, 5000),
    }
    
    relHistDict = {
        'relDiffRecoGenPt_muons' : TH1F("relDiffRecoGenPt_muons", "", 50, -1., 1.),
        'relDiffRecoGenPt_splitMuons' : TH1F("relDiffRecoGenPt_splitMuons", "",  50, -1., 1.)
    }
    
    # Load the trees
    for key, tc in treeDict.items():
        for file_match in glob.glob(indir + '/*{}*.root'.format(key)):
            tc.Add(file_match)
        # tc.Print()
        
    # tree loop
    for key, tc in treeDict.items():
        if verbose > 2: print("!!!!! ===== {} ===== ".format(key))
        
        # Gen - Reco muon collection match
        muon_match_index = recoMuonGenMatch(tc, verbose)
        track_match_index = recoTrackGenMatch(tc, verbose)
        
        if verbose > 2: print(len(muon_match_index), tc.GetEntries())
        
        muon_key = key.split("_")[0]
        track_key = key.split("_")[1]
        if verbose > 2: print("muon_key:~~{}~~".format(muon_key))
        if verbose > 2: print("track_key:~~{}~~".format(track_key))
        
        # Fill genPt
        if muon_key.find("muons")>=0:
            for evt in tc:
                histDict['genPt'].Fill(evt.gen_pt[1])
        
        # Fill reco muon/track pT
        for hkey in histDict:
            if muon_key != '' and hkey.find(muon_key)>=0:
                if verbose > 2: print("found", muon_key, "in", hkey)
                for ind, evt in enumerate(tc):
                    if muon_match_index[ind] != -1: 
                        if verbose > 3: print("plotting: >> match index" , muon_match_index[ind])
                        if verbose > 3: print("plotting: >> muon pt" , evt.muon_pt)
                        histDict[hkey].Fill(evt.muon_pt[muon_match_index[ind]])
            if track_key != '' and hkey.find(track_key)>=0:
                if verbose > 2: print("found", track_key, "in", hkey)
                for ind, evt in enumerate(tc):
                    if track_match_index[ind] != -1: 
                        if verbose > 3: print("plotting: >> match index" , track_match_index[ind])
                        if verbose > 3: print("plotting: >> muon pt" , evt.track_pt)
                        histDict[hkey].Fill(evt.track_pt[track_match_index[ind]])
                        
        # Fill Reco-Gen rel diff
        for hkey in relHistDict:
            if hkey.find(key)>=0:
                for ind, evt in enumerate(tc):
                    if muon_match_index[ind] != -1: 
                        relDiff = (evt.muon_pt[muon_match_index[ind]] - evt.gen_pt[1]) / evt.gen_pt[1]
                        relHistDict[hkey].Fill(relDiff)
            
    gStyle.SetOptStat(0)
    tdrStyle.SetTDRStyle()
    
    # pT
    can = TCanvas("can", "", 800, 600)
    
    histDict = {k: v for k, v in sorted(histDict.items(), key = lambda item: item[1].GetMaximum(), reverse=True)}
    histList = list(histDict.items())
    
    integrals = []
    for i in range(len(histList)):
        hist = histList[i][1]
        if i == 0:
            hist.GetXaxis().SetTitle("p_{T} [GeV]")
            hist.GetYaxis().SetTitle("A. U.")
        hist.SetLineColor(tdrStyle.colors['color_comp{}'.format(i+1)])
        integrals.append(hist.Integral())
        hist.Scale(1./hist.Integral())
        hist.SetLineWidth(2)

    # draw
    can.cd()
    # can.SetLogy()
    for key, hist in histDict.items():
        if i == 0:
            # hist.SetMinimum(1e-3)
            hist.Draw("hist")
        else:
            hist.Draw("hist same")
            
    leg = TLegend(0.25, 0.75, 0.8, 0.94)
    for i, (k, hist) in enumerate(histDict.items()):
        if k.find("gen")>=0:
            leg.AddEntry(hist, "gen : "+str(integrals[i]))
    for i, (k, hist) in enumerate(histDict.items()):
        if not k.find("gen")>=0:
            leg.AddEntry(hist, k.split("_")[1] + ": " + str(integrals[i]))

    leg.Draw("same")
    can.SaveAs(outdir + "/" + "pt_maxPt_sig.png")
    
    '''
    # rel diff reco - gen pT
    can2 = TCanvas("can2", "", 800, 600)
    
    relHistDict = {k: v for k, v in sorted(relHistDict.items(), key = lambda item: item[1].GetMaximum(), reverse=True)}
    relHistList = list(relHistDict.items())
    # print(relHistList)
    
    relHistList[0][1].GetXaxis().SetTitle("(Reco p_{T} #minus Gen p_{T}) / Gen p_{T}")
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
    can2.SaveAs(outdir + "/" + "relDiffRecoGenPt_maxPt_sig.png")
    '''
    
if __name__ == "__main__":
    main()