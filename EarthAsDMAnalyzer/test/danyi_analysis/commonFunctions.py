import copy
from ROOT import TChain, TMath, TH1F, TCanvas, gStyle, TLegend
import tdrStyle

def vprint(verbose, level, msg):
    if verbose >= level:
        print(msg)

def drawNormHistDict(histDict, xtitle, filename, outdir, logy=False, maxy=0.15):
    gStyle.SetOptStat(0)
    tdrStyle.SetTDRStyle()
    
    can = TCanvas("can_{}".format(filename), "", 800, 600)
    
    # histDict = {k: v for k, v in sorted(histDict.items(), key = lambda item: item[1].GetMaximum()/item[1].Integral(), reverse=True)}
    histList = list(histDict.items())
    
    # draw
    can.cd()
    if logy: can.SetLogy()
    
    integrals = []
    for i in range(len(histList)):
        hist = histList[i][1]
        if i == 0:
            hist.GetXaxis().SetTitle(xtitle)
            hist.GetYaxis().SetTitle("A. U.")
        hist.SetLineColor(tdrStyle.colors['color_comp{}'.format(i+1)])
        integrals.append(hist.Integral())
        if hist.Integral() > 0: hist.Scale(1./hist.Integral())
        hist.SetLineWidth(2)
    
    # set range has to be called after Scale()
    for i in range(len(histList)):
        hist = histList[i][1]
        if filename.find("pt")>=0 and logy:
            print("Setting pt range!!")
            hist.GetYaxis().SetRangeUser(1e-3, maxy)
        else:
            hist.SetMaximum(maxy)

    can.Update()
    for key, hist in histDict.items():
        if i == 0:
            # hist.SetMinimum(1e-4)
            hist.Draw("hist")
        else:
            hist.Draw("hist same")
            
    # leg = TLegend(0.25, 0.7, 0.8, 0.94)
    leg = TLegend(0.25, 0.75, 0.85, 0.94)
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
            # leg.AddEntry(hist, k + ": {:.0f}".format(integrals[i]))
            leg.AddEntry(hist, "muon DT match with " + k + ": {:.0f}".format(integrals[i]))

    leg.Draw("same")
    can.SaveAs(outdir + "/" + filename + ".png")

def drawNormHistDict_CompN(histDict, xtitle, filename, outdir, logy=False, maxy=0.15, compN=0):
    gStyle.SetOptStat(0)
    tdrStyle.SetTDRStyle()
    
    can = TCanvas("can_{}".format(filename), "", 800, 600)
    
    # histDict = {k: v for k, v in sorted(histDict.items(), key = lambda item: item[1].GetMaximum()/item[1].Integral(), reverse=True)}
    histList = list(histDict.items())
    
    # draw
    can.cd()
    if logy: can.SetLogy()
    
    integrals = copy.deepcopy(histDict)
    for i in range(compN):
        hist = histList[i][1]
        hist2 = histList[i+compN][1]
        if i == 0:
            hist.GetXaxis().SetTitle(xtitle)
            hist.GetYaxis().SetTitle("A. U.")
        hist.SetLineColor(tdrStyle.colors['color_comp{}'.format(i+1)])
        hist2.SetLineColor(tdrStyle.colors['color_comp{}'.format(i+1)])
        hist.SetLineStyle(2)  # dashed line
        integrals[histList[i][0]] = hist.Integral()
        integrals[histList[i+compN][0]] = hist2.Integral()
        # integrals.append(hist.Integral())
        # integrals2.append(hist2.Integral())
        if hist.Integral() > 0: hist.Scale(1./hist.Integral())
        if hist2.Integral() > 0: hist2.Scale(1./hist2.Integral())
        hist.SetLineWidth(2)
        hist2.SetLineWidth(2)
    
    # set range has to be called after Scale()
    for i in range(len(histList)):
        hist = histList[i][1]
        if filename.find("pt")>=0 and logy:
            print("Setting pt range!!")
            hist.GetYaxis().SetRangeUser(1e-3, maxy)
        else:
            hist.SetMaximum(maxy)

    can.Update()
    for key, hist in histDict.items():
        if i == 0:
            # hist.SetMinimum(1e-4)
            hist.Draw("hist")
        else:
            hist.Draw("hist same")
            
    # leg = TLegend(0.25, 0.7, 0.8, 0.94)
    leg = TLegend(0.25, 0.75, 0.85, 0.94)
    for i, (k, hist) in enumerate(histDict.items()):
        if k.find("gen")>=0:
            leg.AddEntry(hist, "{}: {:.0f}".format(k, integrals[k]))
    for i, (k, hist) in enumerate(histDict.items()):
        if k.find("muons")>=0:
            leg.AddEntry(hist, "muons: {:.0f}".format(integrals[k]))
    for i, (k, hist) in enumerate(histDict.items()):
        if k.find("DT")>=0:
            leg.AddEntry(hist, "muons DT avg: {:.0f}".format(integrals[k]))
    for i, (k, hist) in enumerate(histDict.items()):
        if not k.find("gen")>=0 and not k.find("muons")>=0 and not k.find("DT")>=0:
            # leg.AddEntry(hist, k + ": {:.0f}".format(integrals[i]))
            leg.AddEntry(hist, k + ": {:.0f}".format(integrals[k]))

    leg.Draw("same")
    can.SaveAs(outdir + "/" + filename + ".png")
    
def dR(eta1, phi1, eta2, phi2):
    deta = eta1 - eta2
    dphi = phi1 - phi2
    while dphi > TMath.Pi():
        dphi -= 2 * TMath.Pi()
    while dphi <= -TMath.Pi():
        dphi += 2 * TMath.Pi()
    return TMath.Sqrt(deta * deta + dphi * dphi)


# match reco muon with gen muon (index 1)
def recoMuonGenMatch(tree, verbose, method='mindR', drawOpt=False, outdir=None, key=None):
    ''''
    method: matching method, mindR or maxPt
    drawOpt: if True, will produce dR hist
    '''
    match_ind = []
    h_min_dR = TH1F("h_min_dR", "", 50, 0., 5.)
    for ind, evt in enumerate(tree):
        if verbose > 2 and ind % 10000 == 0: print("Processing event", ind)
        vprint(verbose, 4, "---------------- Event: {} ----------------".format(ind) )
        if evt.muon_n > 0:
            genEta = evt.gen_eta[1]
            genPhi = evt.gen_phi[1]
            sum_ind = 0
            min_dR = 9999.
            min_dR_muonInd = -1
            max_pt = 0.
            max_pt_muonInd = -1
            for j in range(evt.muon_n):
                vprint(verbose, 4, "   >> muon: {}".format(j) )
                vprint(verbose, 4, "   >> dtSeg_n: {}".format(evt.muon_dtSeg_n[j]) )
                # vprint(verbose, 4, "   >> dtSeg_eta:", evt.muon_dtSeg_eta[sum_ind:(sum_ind+evt.muon_dtSeg_n[j])])
                # vprint(verbose, 4, "   >> dtSeg_phi:", evt.muon_dtSeg_phi[sum_ind:(sum_ind+evt.muon_dtSeg_n[j])])
                # sumEtaFromDTseg = sum(eta for eta in evt.muon_dtSeg_eta[sum_ind:(sum_ind+evt.muon_dtSeg_n[j])] if eta < 9000.)
                # sumPhiFromDTseg = sum(phi for phi in evt.muon_dtSeg_phi[sum_ind:(sum_ind+evt.muon_dtSeg_n[j])] if phi < 9000.)
                # sum_ind += evt.muon_dtSeg_n[j]
                # muon_avgEtaFromDTseg = sumEtaFromDTseg / evt.muon_dtSeg_n[j] if evt.muon_dtSeg_n[j] > 0 else 9999.
                # muon_avgPhiFromDTseg = sumPhiFromDTseg / evt.muon_dtSeg_n[j] if evt.muon_dtSeg_n[j] > 0 else 9999.
                muon_dR = dR(genEta, genPhi, evt.muon_avgEtaFromDTseg[j], evt.muon_avgPhiFromDTseg[j])
                if verbose > 3: print ("   >> dR:", muon_dR)
                if muon_dR < 1 and muon_dR < min_dR:
                    min_dR = muon_dR
                    min_dR_muonInd = j
                if evt.muon_pt[j] > max_pt:
                    max_pt = evt.muon_pt[j]
                    max_pt_muonInd = j
            if method == 'mindR':
                match_ind.append(min_dR_muonInd)
                h_min_dR.Fill(min_dR)
            elif method == 'maxPt':
                match_ind.append(max_pt_muonInd)
        else:
            match_ind.append(-1)
        # if ind > 10:
        #     break
    if drawOpt:
        can = TCanvas("canvas_{}".format(key),"",800,600)
        h_min_dR.GetXaxis().SetTitle("min dR")
        h_min_dR.GetYaxis().SetTitle("U.A.")
        h_min_dR.Draw("hist")
        can.SaveAs(outdir + "/" + "mindR_{}.png".format(key))
    vprint(verbose, 4, match_ind)
    
    return match_ind