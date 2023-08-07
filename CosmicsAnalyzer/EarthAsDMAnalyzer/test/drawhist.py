# import numpy as np
# import ROOT as r

# def drawHist(h, xtitle, filename):
#     '''
#         Feed in an array, nbins, low & up range, x-title, filename
#         Plot & return the histgram
#     '''
#     # r.gStyle.SetOptStat(0)

#     can = r.TCanvas()
#     can.Draw()
#     h.GetXaxis().SetTitle(xtitle)
#     h.GetYaxis().SetTitle("Events")
#     h.Draw()
#     can.SaveAs(filename+'.pdf')
#     return h

# def drawHist2d(h, xtitle, ytitle, filename):
#     '''
#         Feed in an array, nbins, low & up range, x-title, filename
#         Plot & return the histgram
#     '''
#     r.gStyle.SetOptStat(0)

#     can = r.TCanvas()
#     can.Draw()
#     h.GetXaxis().SetTitle(xtitle)
#     h.GetYaxis().SetTitle(ytitle)
#     h.Draw("colz")
#     can.SaveAs(filename+'.pdf')
#     return h

# def main():
#     # load tree
#     #tree = r.TChain("muonPhiAnalyzer/MCTree")
#     tree = r.TChain("muonPhiAnalyzer/tree")
#     #tree.Add("MCVariables3.root")
#     tree.Add("TRK-Run3Winter23Reco-00009-10000Events.root")
#     nevents = tree.GetEntries()
#     print("nentries = ", nevents)
    
#     h_muonPhi = r.TH1D( 'h_muonPhi', '', 50, -3.14, 3.14)
#     h_muonPt = r.TH1D( 'h_muonPt', '', 100, 0, 3000)
#     h_muonEta = r.TH1D( 'h_muonEta', '', 50, -1.5, 1.5)
    
#     h_muonPtvsPhi = r.TH2D( 'h_muonPtvsPhi', '', 50, -3.14, 3.14, 100, 0, 3000)
#     h_muonPtvsEta = r.TH2D( 'h_muonPtvsEta', '', 50, -1.5, 1.5, 100, 0, 3000)
    
#     event_count = 0
#     pfreq = 1000
    
#     for event in tree:
        
#         if event_count%pfreq == 0:
#             print('Processing Event: %s'%(event_count))
        
#         for phi in event.muonPhi:
#             h_muonPhi.Fill(phi)
#         for pt in event.muonPt:
#             h_muonPt.Fill(pt)
#         for eta in event.muonEta:
#             h_muonEta.Fill(eta)
        
#         for (pt, phi, eta) in zip(event.muonPt, event.muonPhi, event.muonEta):
#             h_muonPtvsPhi.Fill(phi, pt)
#             h_muonPtvsEta.Fill(eta, pt)
        
#         event_count += 1
    
#     drawHist(h_muonPhi, "#phi", "plots/h_muonPhi")
#     drawHist(h_muonPt, "p_{T} [GeV]", "plots/h_muonPt")
#     drawHist(h_muonEta, "#eta", "plots/h_muonEta")
#     drawHist2d(h_muonPtvsPhi, "#phi", "p_{T} [GeV]", "plots/h_muonPtvsPhi")
#     drawHist2d(h_muonPtvsEta, "#eta", "p_{T} [GeV]", "plots/h_muonPtvsEta")
    
# if __name__=="__main__":
#     main()


# import numpy as np
# import ROOT as r

# def loadTree(file_path, tree_name):
#     tree = r.TChain(tree_name)
#     tree.Add(file_path)
#     return tree

# def drawHist(h, xtitle, filename, color):
#     '''
#         Feed in a histogram, x-title, filename, and color
#         Plot & return the histogram
#     '''
#     r.gStyle.SetOptStat(0)

#     can = r.TCanvas()
#     can.Draw()
#     h.GetXaxis().SetTitle(xtitle)
#     h.GetYaxis().SetTitle("Events")
#     h.SetLineColor(color)
#     h.Draw()
#     can.SaveAs(filename+'.pdf')
#     return h

# def main():
#     # Load first tree
#     tree1 = loadTree("MCntuple0to75-10000Events.root", "muonPhiAnalyzer/tree")
#     nevents1 = tree1.GetEntries()
#     print("nentries1 = ", nevents1)
    
#     # Load second tree
#     tree2 = loadTree("MCntuple91to180-10000Events.root", "muonPhiAnalyzer/tree")
#     nevents2 = tree2.GetEntries()
#     print("nentries2 = ", nevents2)
    
#     h_muonPhi1 = r.TH1D('h_muonPhi1', '', 50, -3.14, 3.14)
#     h_muonPhi2 = r.TH1D('h_muonPhi2', '', 50, -3.14, 3.14)
    
#     event_count = 0
#     pfreq = 1000
    
#     for event in tree1:
#         if event_count % pfreq == 0:
#             print('Processing Event from Tree 1: %s' % (event_count))
#         for phi in event.muon_phi:
#             h_muonPhi1.Fill(phi)
#         event_count += 1

#     event_count = 0
#     for event in tree2:
#         if event_count % pfreq == 0:
#             print('Processing Event from Tree 2: %s' % (event_count))
#         for phi in event.muon_phi:
#             h_muonPhi2.Fill(phi)
#         event_count += 1
    
#     h_muonPhi1 = drawHist(h_muonPhi1, "#phi", "plots/h_muonPhi", r.kRed)
#     h_muonPhi2 = drawHist(h_muonPhi2, "#phi", "plots/h_muonPhi", r.kBlue)

# if __name__ == "__main__":
#     main()

# import numpy as np
# import ROOT as r

# def loadTree(file_path, tree_name):
#     tree = r.TChain(tree_name)
#     tree.Add(file_path)
#     return tree

# def drawHist(h, xtitle, filename, color, label):
#     '''
#         Feed in a histogram, x-title, filename, and color
#         Plot & return the TCanvas
#     '''
#     r.gStyle.SetOptStat(0)

#     can = r.TCanvas()
#     can.Draw()
#     h.GetXaxis().SetTitle(xtitle)
#     h.GetYaxis().SetTitle("Events")
#     h.SetLineColor(color)
#     h.Draw()

#     legend = r.TLegend(0.65, 0.7, 0.85, 0.85)  # Define the legend position
#     legend.SetFillColor(0)
#     legend.SetBorderSize(0)
#     legend.AddEntry(h, label, "l")  # Add entry to the legend for the current histogram
#     legend.Draw()

#     can.SaveAs(filename+'.pdf')
#     return can  # Return the TCanvas

# def main():
#     # Load first tree
#     tree1 = loadTree("MCntuple0to75-10000Events.root", "muonPhiAnalyzer/tree")
#     nevents1 = tree1.GetEntries()
#     print("nentries1 = ", nevents1)
    
#     # Load second tree
#     tree2 = loadTree("MCntuple91to180-10000Events.root", "muonPhiAnalyzer/tree")
#     nevents2 = tree2.GetEntries()
#     print("nentries2 = ", nevents2)
    
#     h_muonPhi1 = r.TH1D('h_muonPhi_0to75 comparison', 'Muon #Phi comparison;#phi;yield;', 80, -3.5, 3.5) #-3.14, 3.14)
#     h_muonPhi2 = r.TH1D('h_muonEta2', '', 80, -3.5, 3.5) #-3.14, 3.14)

#     h_muonPt1 = r.TH1D('h_muonPt_0to75', 'Muon Pt comparison;pt;yield;', 80, -20, 20) #-3.14, 3.14)
#     h_muonPt2 = r.TH1D('h_muonEta2', '', 80, -20, 20) #-3.14, 3.14)

#     h_muonEta1 = r.TH1D('h_muonEta_0to75', 'Muon #Eta comparison;#Eta;yield;', 80, -20, 20) #-3.14, 3.14)
#     h_muonEta2 = r.TH1D('h_muonEta2', '', 80, -20, 20) #-3.14, 3.14)

#     h_muon_timeAtIpOutIn1 = r.TH1D('h_muon_timeAtIpOutIn_0to75', 'Muon timeAtIpOutIn comparison;timeAtIpOutIn;yield;', 80, -20, 20) #-3.14, 3.14)
#     h_muon_timeAtIpOutIn2 = r.TH1D('h_muonEta2', '', 80, -20, 20) #-3.14, 3.14)

#     h_muon_timeAtIpInOut1 = r.TH1D('h_muon_timeAtIpInOut_0to75', 'Muon timeAtIpInOut comparison;timeAtIpInOut;yield;', 80, -20, 20) #-3.14, 3.14)
#     h_muon_timeAtIpInOut2 = r.TH1D('h_muonEta2', '', 80, -20, 20) #-3.14, 3.14)

#     h_muon_freeInvBeta1 = r.TH1D('h_muon_freeInvBeta_0to75', 'Muon freeInvBeta comparison;freeInvBeta;yield;', 80, -20, 20) #-3.14, 3.14)
#     h_muon_freeInvBeta2 = r.TH1D('h_muonEta2', '', 80, -20, 20) #-3.14, 3.14)

#     #Can create the TH1D as part of an array and ask them to fill
#     #Create an intermediate step that takes the root file, creates a new root file that is the TH1D outputs, and then graphs that
    
#     event_count = 0
#     pfreq = 1000
    
#     for event in tree1:
#         if event_count % pfreq == 0:
#             print('Processing Event from Tree 1: %s' % (event_count))
#         for phi in event.muon_phi:
#             if phi > -1000 and phi < 1000:
#                 h_muonPhi1.Fill(phi)
#         for pt in event.muon_pt:
#             #if pt > -1000 and pt < 1000:
#             #    h_muonPt1.Fill(pt)
#         for eta in event.muon_eta:
#             if eta > -1000 and eta < 1000:
#                 h_muonEta1.Fill(eta)
#         for FreeInvBeta in event.muons1leg_comb_freeInvBeta:
#             if FreeInvBeta > -1000 and FreeInvBeta < 1000:
#                 h_muon_freeInvBeta1.Fill(FreeInvBeta)
#         for timeAtIpInOut in event.muons1leg_comb_timeAtIpInOut:
#             if timeAtIpInOut > -1000 and timeAtIpInOut < 1000:
#                 h_muon_timeAtIpInOut1.Fill(timeAtIpInOut)
#         for timeAtIpOutIn in event.muons1leg_comb_timeAtIpOutIn:
#             if timeAtIpOutIn > -1000 and timeAtIpOutIn < 1000:
#                 h_muon_timeAtIpOutIn1.Fill(timeAtIpOutIn)
#             #make sure that the inout is not blown up
#         event_count += 1

#     event_count = 0
#     for event in tree2:
#         if event_count % pfreq == 0:
#             print('Processing Event from Tree 2: %s' % (event_count))
#         for phi in event.muon_phi:
#             if phi > -1000 and phi < 1000:
#                 h_muonPhi2.Fill(phi)
#         for pt in event.muon_pt:
#             #if pt > -1000 and pt < 1000:
#             #    h_muonPhi1.Fill(pt)
#         for eta in event.muon_eta:
#             if eta > -1000 and eta < 1000:
#                 h_muonEta2.Fill(eta)
#         for FreeInvBeta in event.muons1leg_comb_freeInvBeta:
#             if FreeInvBeta > -1000 and FreeInvBeta < 1000:
#                 h_muon_freeInvBeta2.Fill(FreeInvBeta)
#         for timeAtIpInOut in event.muons1leg_comb_timeAtIpInOut:
#             if timeAtIpInOut > -1000 and timeAtIpInOut < 1000:
#                 h_muon_timeAtIpInOut2.Fill(timeAtIpInOut)
#         for timeAtIpOutIn in event.muons1leg_comb_timeAtIpOutIn:
#             if timeAtIpOutIn > -1000 and timeAtIpOutIn < 1000:
#                 h_muon_timeAtIpOutIn2.Fill(timeAtIpOutIn)
#         event_count += 1
    
#     # can = drawHist(h_muonPhi1, "#phi", "plots/h_muonPhi", r.kRed, "0 to 75")
#     # h_muonPhi2.Draw("same")
#     # h_muonPhi2.SetLineColor(r.kBlue)
#     # can.SaveAs("plots/h_muonPhi_both.pdf")

#     can = r.TCanvas()
#     can.Draw()
#     h_muonPhi1.SetLineColor(r.kRed)
#     h_muonPhi1.Draw()

#     h_muonPhi2.SetLineColor(r.kBlue)
#     h_muonPhi2.Draw("same")

#     legend = r.TLegend(0.15, 0.4, 0.35, 0.55)
#     legend.SetFillColor(0)
#     legend.SetBorderSize(1)
#     legend.AddEntry(h_muonPhi1, "0 to 75", "l")
#     #Add + .getmean()
#     legend.AddEntry(h_muonPhi2, "91 to 180", "l")
#     legend.Draw()

#     can.SaveAs("plots/muon_phi_both.pdf")
#     can.SaveAs("plots/muon_phi_both.png")
#     can.close()

#     can2 = r.TCanvas()
#     can2.Draw()
#     h_muonPt1.SetLineColor(r.kRed)
#     h_muonPt1.Draw()

#     h_muonPt2.SetLineColor(r.kBlue)
#     h_muonPt2.Draw("same")

#     legend = r.TLegend(0.15, 0.4, 0.35, 0.55)
#     legend.SetFillColor(0)
#     legend.SetBorderSize(1)
#     legend.AddEntry(h_muonPt1, "0 to 75", "l")
#     #Add + .getmean()
#     legend.AddEntry(h_muonPt1, "91 to 180", "l")
#     legend.Draw()

#     #h_muonPhi1 = drawHist(h_muonPhi1, "#phi", "plots/h_muonPhi", r.kRed,  "0 to 75")
#     #h_muonPhi2 = drawHist(h_muonPhi2, "#phi", "plots/h_muonPhi2", r.kBlue, "91 to 180")


#     can2.SaveAs("plots/muon_pt_both.pdf")
#     can2.SaveAs("plots/muon_pt_both.png")


# if __name__ == "__main__":
#     main()


import ROOT as r
import argparse

def loadTree(file_path, tree_name):
    tree = r.TChain(tree_name)
    tree.Add(file_path)
    return tree

def drawHist(hist, x_label, save_path, color, title, draw_option="HIST", legend_label=""):



    r.gStyle.SetOptStat(0)
    hist.Draw(draw_option)
    hist.SetLineColor(color)
    hist.SetLineWidth(2)
    hist.GetXaxis().SetTitle(x_label)
    hist.GetYaxis().SetTitle("Events")
    hist.SetTitle(title)

    # if legend_label:
    #     legend = r.TLegend(0.7, 0.75, 0.9, 0.9)
    #     legend.AddEntry(hist, legend_label, "l")
    #     legend.Draw()

    return hist

def drawCanvas(hist1, hist2, legendItem1, legendItem2, title,x1,x2,y1,y2):

    can = r.TCanvas()
    can.Draw()
    hist1.SetLineColor(r.kRed)
    hist1.Draw()

    hist2.SetLineColor(r.kBlue)
    hist2.Draw("same")

    legend = r.TLegend(x1,x2,y1,y2)  # Define the legend position
    legend.SetFillColor(0)
    legend.SetBorderSize(1)

    legend.AddEntry(hist1, legendItem1, "l")  # Add entry to the legend for the current histogram
    legend.AddEntry(hist2, legendItem2, "l")
    legend.Draw()

    can.SaveAs(f"plots/test7/{title}_both.pdf")
    can.SaveAs(f"plots/test7/{title}_both.png")

    return can


def draw2DHist(hist1, legendItem1, legendItem2, title,x1,x2,y1,y2):

    can = r.TCanvas()
    can.Draw()
    hist1.SetLineColor(r.kRed)
    hist1.Draw()

    # hist2.SetLineColor(r.kBlue)
    # hist2.Draw("same")

    legend = r.TLegend(x1,x2,y1,y2)  # Define the legend position
    legend.SetFillColor(0)
    legend.SetBorderSize(1)

    legend.AddEntry(hist1, legendItem1, "l")  # Add entry to the legend for the current histogram
    # legend.AddEntry(hist2, legendItem2, "l")
    legend.Draw()

    can.SaveAs(f"plots/2DhistTest/{title}_both.pdf")
    can.SaveAs(f"plots/2DhistTest/{title}_both.png")

    return can

def main():

    # Load first tree
    #tree1 = loadTree("MCntuple0to75-10000Events.root", "muonPhiAnalyzer/tree")
    tree1 = loadTree("ImportantSampleProductions/4-MCAnalysis/MCntuple-0to75-muons.root", "muonPhiAnalyzer/tree")
    nevents1 = tree1.GetEntries()
    print("nentries1 = ", nevents1)
    
    # Load second tree
    #tree2 = loadTree("MCntuple91to180-10000Events.root", "muonPhiAnalyzer/tree")
    tree2 = loadTree("ImportantSampleProductions/4-MCAnalysis/MCntuple-91to180-muons.root", "muonPhiAnalyzer/tree")
    nevents2 = tree2.GetEntries()
    print("nentries2 = ", nevents2)

    h_muonPhi1 = r.TH1F("h_muonPhi1", "Muon #Phi Comparison;#phi;yield", 100, -4, 4)
    h_muonPhi2 = r.TH1F("h_muonPhi2", "Muon #phi - 91 to 180;#phi;yield", 100, -4, 4)
    h_muonPt1 = r.TH1F("h_muonPt1", "Muon pT Comparison;pt;yield", 100, 0, 300)
    h_muonPt2 = r.TH1F("h_muonPt2", "Muon pT Comparison;pt;yield", 100, 0, 300)
    h_muonEta1 = r.TH1F("h_muonEta1", "Muon #eta Comparison;#eta;yield", 100, -2.5, 2.5)
    h_muonEta2 = r.TH1F("h_muonEta2", "Muon #eta Comparison;#eta;yield", 100, -2.5, 2.5)
    h_muon_freeInvBeta1 = r.TH1F("h_muon_freeInvBeta1", "Muon FreeInverseBeta Comparison;freeInvBeta;yield", 100, -50, 50)
    h_muon_freeInvBeta2 = r.TH1F("h_muon_freeInvBeta2", "Muon FreeInverseBeta Comparison;freeInvBeta;yield", 100, -50, 50)
    h_muon_timeAtIpInOut1 = r.TH1F("h_muon_timeAtIpInOut1", "Muon TimeAtIpInOut Comparison;timeAtIpInOut;yield", 100, -500, 500)
    h_muon_timeAtIpInOut2 = r.TH1F("h_muon_timeAtIpInOut2", "Muon TimeAtIpInOut Comparison;timeAtIpInOut;yield", 100, -500, 500)
    h_muon_timeAtIpOutIn1 = r.TH1F("h_muon_timeAtIpOutIn1", "Muon TimeAtIpOutIn Comparison;timeAtIpOutIn;yield", 100, -500, 500)
    h_muon_timeAtIpOutIn2 = r.TH1F("h_muon_timeAtIpOutIn2", "Muon TimeAtIpOutIn Comparison;timeAtIpOutIn;yield", 100, -500, 500)

    h_muonGlobalYvsTime1 = r.TH2D("h_muonGlobalYvsTime1", "Y vs time;#phi;time", 50, -1000, 1000, 50, -200, 200)
    h_muonGlobalYvsTime2 = r.TH2D("h_muonGlobalYvsTime2", "Muon #phi - 91 to 180;#phi;yield", 50, -4, 4, 50, 0, 3000)


    # h_muon1LegPhi1 = r.TH1F("h_muon1LegPhi1", "muon1Leg #phi comparison;#phi;yield", 100, -4, 4)
    # h_muon1LegPhi2 = r.TH1F("h_muon1LegPhi2", "muon1Leg #phi comparison;#phi;yield", 100, -4, 4)
    # h_muon1LegPt1 = r.TH1F("h_muon1LegPt1", "muon1Leg pT comparison;pt;yield", 100, 0, 300)
    # h_muon1LegPt2 = r.TH1F("h_muon1LegPt2", "muon1Leg pT comparison;pt;yield", 100, 0, 300)
    # h_muon1LegEta1 = r.TH1F("h_muon1LegEta1", "muon1Leg #eta comparison;#eta;yield", 100, -2.5, 2.5)
    # h_muon1LegEta2 = r.TH1F("h_muon1LegEta2", "muon1Leg #eta comparison;#eta;yield", 100, -2.5, 2.5)
    # h_muon1Leg_freeInvBeta1 = r.TH1F("h_muon1Leg_freeInvBeta1", "muon1Leg FreeInverseBeta comparison;freeInvBeta;yield", 100, -50, 50)
    # h_muon1Leg_freeInvBeta2 = r.TH1F("h_muon1Leg_freeInvBeta2", "muon1Leg FreeInverseBeta comparison;freeInvBeta;yield", 100, -50, 50)
    # h_muon1Leg_timeAtIpInOut1 = r.TH1F("h_muon1Leg_timeAtIpInOut1", "muon1Leg TimeAtIpInOut comparison;timeAtIpInOut;yield", 100, -500, 500) #0 to 75
    # h_muon1Leg_timeAtIpInOut2 = r.TH1F("h_muon1Leg_timeAtIpInOut2", "muon1Leg TimeAtIpInOut comparison;timeAtIpInOut;yield", 100, -500, 500) #91 to 180
    # h_muon1Leg_timeAtIpOutIn1 = r.TH1F("h_muon1Leg_timeAtIpOutIn1", "muon1Leg TimeAtIpOutIn comparison;timeAtIpOutIn;yield", 100, -500, 500)
    # h_muon1Leg_timeAtIpOutIn2 = r.TH1F("h_muon1Leg_timeAtIpOutIn2", "muon1Leg TimeAtIpOutIn comparison;timeAtIpOutIn;yield", 100, -500, 500)



    # h_splitmuonPhi1 = r.TH1F("h_splitmuonPhi1", "splitMuon #Phi Comparison;#phi;yield", 100, -4, 4)
    # h_splitmuonPhi2 = r.TH1F("h_splitmuonPhi2", "#phi - 91 to 180;#phi;yield", 100, -4, 4)
    # h_splitmuonPt1 = r.TH1F("h_splitmuonPt1", "splitMuon pT Comparison;pt;yield", 100, 0, 300)
    # h_splitmuonPt2 = r.TH1F("h_splitmuonPt2", "splitMuon pT Comparison;pt;yield", 100, 0, 300)
    # h_splitmuonEta1 = r.TH1F("h_splitmuonEta1", "splitmuon #eta Comparison;#eta;yield", 100, -2.5, 2.5)
    # h_splitmuonEta2 = r.TH1F("h_splitmuonEta2", "splitmuon #eta Comparison;#eta;yield", 100, -2.5, 2.5)
    # h_splitmuon_freeInvBeta1 = r.TH1F("h_splitmuon_freeInvBeta1", "splitmuon FreeInverseBeta Comparison;freeInvBeta;yield", 100, -50, 50)
    # h_splitmuon_freeInvBeta2 = r.TH1F("h_splitmuon_freeInvBeta2", "splitmuon FreeInverseBeta Comparison;freeInvBeta;yield", 100, -50, 50)
    # h_splitmuon_timeAtIpInOut1 = r.TH1F("h_splitmuon_timeAtIpInOut1", "splitmuon TimeAtIpInOut Comparison;timeAtIpInOut;yield", 100, -500, 500)
    # h_splitmuon_timeAtIpInOut2 = r.TH1F("h_splitmuon_timeAtIpInOut2", "splitmuon TimeAtIpInOut Comparison;timeAtIpInOut;yield", 100, -500, 500)
    # h_splitmuon_timeAtIpOutIn1 = r.TH1F("h_splitmuon_timeAtIpOutIn1", "splitmuon TimeAtIpOutIn Comparison;timeAtIpOutIn;yield", 100, -500, 500)
    # h_splitmuon_timeAtIpOutIn2 = r.TH1F("h_splitmuon_timeAtIpOutIn2", "splitmuon TimeAtIpOutIn Comparison;timeAtIpOutIn;yield", 100, -500, 500)


    # h_lhcSTAMuonPhi1 = r.TH1F("h_lhcSTAMuonsPhi1", "lhcSTAMuons #Phi Comparison;#phi;yield", 100, -4, 4)
    # h_lhcSTAMuonPhi2 = r.TH1F("h_lhcSTAMuonsPhi2", "lhcSTAMuons #phi - 91 to 180;#phi;yield", 100, -4, 4)
    # h_lhcSTAMuonPt1 = r.TH1F("h_lhcSTAMuonPt1", "lhcSTAMuons pT Comparison;pt;yield", 100, 0, 300)
    # h_lhcSTAMuonPt2 = r.TH1F("h_lhcSTAMuonPt2", "lhcSTAMuons pT Comparison;pt;yield", 100, 0, 300)
    # h_lhcSTAMuonEta1 = r.TH1F("h_lhcSTAMuonEta1", "lhcSTAMuons #eta Comparison;#eta;yield", 100, -2.5, 2.5)
    # h_lhcSTAMuonEta2 = r.TH1F("h_lhcSTAMuonEta2", "lhcSTAMuons #eta Comparison;#eta;yield", 100, -2.5, 2.5)
    # h_lhcSTAMuon_freeInvBeta1 = r.TH1F("h_lhcSTAMuon_freeInvBeta1", "lhcSTAMuons FreeInverseBeta Comparison;freeInvBeta;yield", 100, -50, 50)
    # h_lhcSTAMuon_freeInvBeta2 = r.TH1F("h_lhcSTAMuon_freeInvBeta2", "lhcSTAMuons FreeInverseBeta Comparison;freeInvBeta;yield", 100, -50, 50)
    # h_lhcSTAMuon_timeAtIpInOut1 = r.TH1F("h_lhcSTAMuon_timeAtIpInOut1", "lhcSTAMuons TimeAtIpInOut Comparison;timeAtIpInOut;yield", 100, -500, 500)
    # h_lhcSTAMuon_timeAtIpInOut2 = r.TH1F("h_lhcSTAMuon_timeAtIpInOut2", "lhcSTAMuons TimeAtIpInOut Comparison;timeAtIpInOut;yield", 100, -500, 500)
    # h_lhcSTAMuon_timeAtIpOutIn1 = r.TH1F("h_lhcSTAMuon_timeAtIpOutIn1", "lhcSTAMuons TimeAtIpOutIn Comparison;timeAtIpOutIn;yield", 100, -500, 500)
    # h_lhcSTAMuon_timeAtIpOutIn2 = r.TH1F("h_lhcSTAMuon_timeAtIpOutIn2", "lhcSTAMuons TimeAtIpOutIn Comparison;timeAtIpOutIn;yield", 100, -500, 500)




    event_count = 0
    pfreq = 1000

    for tree in [tree1, tree2]:
        for event in tree:
            #if event_count % pfreq == 0:
                #print('Processing Event from Tree %d: %d' % (tree, event_count))

            for phi in event.muon_phi:
                if phi > -1000 and phi < 1000: # and event.muon_tofMap_found == 1:
                    if tree == tree1:
                        h_muonPhi1.Fill(phi)
                    else:
                        h_muonPhi2.Fill(phi)

            for pt in event.muon_pt:
                if pt > -1000 and pt < 1000: # and event.muon_tofMap_found == 1:
                    if tree == tree1:
                        h_muonPt1.Fill(pt)
                    else:
                        h_muonPt2.Fill(pt)

            for eta in event.muon_eta:
                if eta > -1000 and eta < 1000: # and event.muon_tofMap_found == 1:
                    if tree == tree1:
                        h_muonEta1.Fill(eta)
                    else:
                        h_muonEta2.Fill(eta)

            for FreeInvBeta in event.muon_comb_freeInvBeta:
                if FreeInvBeta > -1000 and FreeInvBeta < 1000 and FreeInvBeta != 0: #event.muon_tofMap_found == 1:  #and FreeInvBeta !=0:
                        if tree == tree1:
                            h_muon_freeInvBeta1.Fill(FreeInvBeta)
                        else:
                            h_muon_freeInvBeta2.Fill(FreeInvBeta)

            for timeAtIpInOut in event.muon_comb_timeAtIpInOut:
                if timeAtIpInOut > -1000 and timeAtIpInOut < 1000 and timeAtIpInOut != 0: #event.muon_tofMap_found == 1: #and timeAtIpInOut != 0:
                    if tree == tree1:
                        h_muon_timeAtIpInOut1.Fill(timeAtIpInOut)
                    else:
                        h_muon_timeAtIpInOut2.Fill(timeAtIpInOut)

            for timeAtIpOutIn in event.muon_comb_timeAtIpOutIn:
                if timeAtIpOutIn > -1000 and timeAtIpOutIn < 1000 and timeAtIpOutIn != 0: #event.muon_tofMap_found == 1:
                    if tree == tree1:
                        h_muon_timeAtIpOutIn1.Fill(timeAtIpOutIn)
                    else:
                        h_muon_timeAtIpOutIn2.Fill(timeAtIpOutIn)

            # for dtSeg_globY in event.muon_dtSeg_globY:
            #     if dtSeg_globY > -1000 and dtSeg_globY < 1000 and dtSeg_globY != 0: #event.muon_tofMap_found == 1:
            #         if tree == tree1:
            #             h_muonGlobalYvsTime1.Fill(dtSeg_globY)
            #         else:
            #             h_muonGlobalYvsTime1.Fill(dtSeg_globY)

            for (dtSeg_t0timing, dtSeg_globY) in zip (event.muon_dtSeg_t0timing, event.muon_dtSeg_globY):
                if dtSeg_t0timing > -1000 and dtSeg_t0timing < 1000 and dtSeg_t0timing != 0: #event.muon_tofMap_found == 1:
                    if tree == tree1:
                        h_muonGlobalYvsTime1.Fill(dtSeg_globY, dtSeg_t0timing)
                    else:
                        h_muonGlobalYvsTime2.Fill(dtSeg_t0timing, dtSeg_globY)

#         for (pt, phi, eta) in zip(event.muonPt, event.muonPhi, event.muonEta):
#             h_muonPtvsPhi.Fill(phi, pt)
#             h_muonPtvsEta.Fill(eta, pt)



#             #1 Leg Values
#             for phi in event.muons1leg_phi:
#                 if phi > -1000 and phi < 1000: # and event.muons1leg_tofMap_found == 1:
#                     if tree == tree1:
#                         h_muon1LegPhi1.Fill(phi)
#                     else:
#                         h_muon1LegPhi2.Fill(phi)

#             for pt in event.muons1leg_pt:
#                 if pt > -1000 and pt < 1000: # and event.muons1leg_tofMap_found == 1:
#                 #    h_muonPt1.Fill(pt)
#                     if tree == tree1:
#                         h_muon1LegPt1.Fill(pt)
#                     else:
#                         h_muon1LegPt2.Fill(pt)

#             for eta in event.muons1leg_eta:
#                 if eta > -1000 and eta < 1000: # and event.muons1leg_tofMap_found == 1:
#                     if tree == tree1:
#                         h_muon1LegEta1.Fill(eta)
#                     else:
#                         h_muon1LegEta2.Fill(eta)

#             for FreeInvBeta in event.muons1leg_comb_freeInvBeta:
#                 if FreeInvBeta > -1000 and FreeInvBeta < 1000 and FreeInvBeta == 0: # event.muons1leg_tofMap_found == 1:
#                     if tree == tree1:
#                         h_muon1Leg_freeInvBeta1.Fill(FreeInvBeta)
#                     else:
#                         h_muon1Leg_freeInvBeta2.Fill(FreeInvBeta)

#             for timeAtIpInOut in event.muons1leg_comb_timeAtIpInOut:
#                 if timeAtIpInOut > -1000 and timeAtIpInOut < 1000 and timeAtIpInOut == 0: # event.muons1leg_tofMap_found == 1:
#                     if tree == tree1:
#                         h_muon1Leg_timeAtIpInOut1.Fill(timeAtIpInOut)
#                     else:
#                         h_muon1Leg_timeAtIpInOut2.Fill(timeAtIpInOut)

#             for timeAtIpOutIn in event.muons1leg_comb_timeAtIpOutIn:
#                 if timeAtIpOutIn > -1000 and timeAtIpOutIn < 1000 and timeAtIpOutIn == 0: #event.muons1leg_tofMap_found == 1:
#                     if tree == tree1:
#                         h_muon1Leg_timeAtIpOutIn1.Fill(timeAtIpOutIn)
#                     else:
#                         h_muon1Leg_timeAtIpOutIn2.Fill(timeAtIpOutIn)



#             #Splitmuons
#             for phi in event.splitMuons_phi:
#                 if phi > -1000 and phi < 1000: # and event.splitMuons_tofMap_found == 1:
#                     if tree == tree1:
#                         h_splitmuonPhi1.Fill(phi)
#                     else:
#                         h_splitmuonPhi2.Fill(phi)

#             for pt in event.splitMuons_pt:
#                 if pt > -1000 and pt < 1000: # and event.splitMuons_tofMap_found == 1:
#                 #    h_muonPt1.Fill(pt)
#                     if tree == tree1:
#                         h_splitmuonPt1.Fill(pt)
#                     else:
#                         h_splitmuonPt2.Fill(pt)

#             for eta in event.splitMuons_eta:
#                 if eta > -1000 and eta < 1000: # and event.splitMuons_tofMap_found == 1:
#                     if tree == tree1:
#                         h_splitmuonEta1.Fill(eta)
#                     else:
#                         h_splitmuonEta2.Fill(eta)

#             for FreeInvBeta in event.splitMuons_comb_freeInvBeta:
#                 if FreeInvBeta > -1000 and FreeInvBeta < 1000 and FreeInvBeta == 0: # event.splitMuons_tofMap_found == 1:
#                     if tree == tree1:
#                         h_lhcSTAMuon_freeInvBeta1.Fill(FreeInvBeta)
#                     else:
#                         h_lhcSTAMuon_freeInvBeta2.Fill(FreeInvBeta)

#             for timeAtIpInOut in event.splitMuons_comb_timeAtIpInOut:
#                 if timeAtIpInOut > -1000 and timeAtIpInOut < 1000 and timeAtIpInOut == 0: #event.splitMuons_tofMap_found == 1:
#                     if tree == tree1:
#                         h_splitmuon_timeAtIpInOut1.Fill(timeAtIpInOut)
#                     else:
#                         h_splitmuon_timeAtIpInOut2.Fill(timeAtIpInOut)

#             for timeAtIpOutIn in event.splitMuons_comb_timeAtIpOutIn:
#                 if timeAtIpOutIn > -1000 and timeAtIpOutIn < 1000 and timeAtIpOutIn == 0: #event.splitMuons_tofMap_found == 1:
#                     if tree == tree1:
#                         h_splitmuon_timeAtIpOutIn1.Fill(timeAtIpOutIn)
#                     else:
#                         h_splitmuon_timeAtIpOutIn2.Fill(timeAtIpOutIn)




# #lhcSTAMuons
#             for phi in event.lhcSTAMuons_phi:
#                 if phi > -1000 and phi < 1000: #and event.lhcSTAMuons_tofMap_found == 1:
#                     if tree == tree1:
#                         h_lhcSTAMuonPhi1.Fill(phi)
#                     else:
#                         h_lhcSTAMuonPhi2.Fill(phi)

#             for pt in event.lhcSTAMuons_pt:
#                 if pt > -1000 and pt < 1000: #and event.lhcSTAMuons_tofMap_found == 1:
#                     if tree == tree1:
#                         h_lhcSTAMuonPt1.Fill(pt)
#                     else:
#                         h_lhcSTAMuonPt2.Fill(pt)

#             for eta in event.lhcSTAMuons_eta:
#                 if eta > -1000 and eta < 1000: #event.lhcSTAMuons_tofMap_found == 1:
#                     if tree == tree1:
#                         h_lhcSTAMuonEta1.Fill(eta)
#                     else:
#                         h_lhcSTAMuonEta2.Fill(eta)

#             for FreeInvBeta in event.lhcSTAMuons_comb_freeInvBeta:
#                 if FreeInvBeta > -1000 and FreeInvBeta < 1000 and FreeInvBeta != 0: # event.lhcSTAMuons_tofMap_found == 1:
#                     if tree == tree1:
#                         h_splitmuon_freeInvBeta1.Fill(FreeInvBeta)
#                     else:
#                         h_splitmuon_freeInvBeta2.Fill(FreeInvBeta)

#             for timeAtIpInOut in event.lhcSTAMuons_comb_timeAtIpInOut and event.lhcSTAMuons_tofMap_found:
#                 if timeAtIpInOut > -1000 and timeAtIpInOut < 1000 and lhcSTAMuons_tofMap_found == 1:
#                     if event.lhcSTAMuons_tofMap_found == 1:
#                         if tree == tree1:
#                             h_lhcSTAMuon_timeAtIpInOut1.Fill(timeAtIpInOut)
#                         else:
#                             h_lhcSTAMuon_timeAtIpInOut2.Fill(timeAtIpInOut)

#             for timeAtIpOutIn in event.lhcSTAMuons_comb_timeAtIpOutIn:
#                 if timeAtIpOutIn > -1000 and timeAtIpOutIn < 1000 and event.lhcSTAMuons_tofMap_found == 1:
#                     if tree == tree1:
#                         h_lhcSTAMuon_timeAtIpOutIn1.Fill(timeAtIpOutIn)
#                     else:
#                         h_lhcSTAMuon_timeAtIpOutIn2.Fill(timeAtIpOutIn)

            event_count += 1

    # can = r.TCanvas()
    # can.Draw()
    # h_muonPhi1 = drawHist(h_muonPhi1, "#phi", "plots/h_muonPhi", r.kRed, "0 to 75")
    # h_muonPhi2.Draw("same")
    # h_muonPhi2.SetLineColor(r.kBlue)
    # h_muonPhi1 = drawHist(h_muonPhi1, "#phi", "plots/h_muonPhi", r.kRed,  "0 to 75")
    # h_muonPhi2 = drawHist(h_muonPhi2, "#phi", "plots/h_muonPhi2", r.kBlue, "91 to 180")

    # legend = r.TLegend(0.15, 0.4, 0.35, 0.55)
    # legend.SetFillColor(0)
    # legend.SetBorderSize(1)
    # legend.AddEntry(h_muonPhi1, "0 to 75", "l")
    # legend.AddEntry(h_muonPhi2, "91 to 180", "l")
    # legend.Draw()

    # can = r.TCanvas()
    # can.Draw()
    # h_muonPhi1.SetLineColor(r.kRed)
    # h_muonPhi1.Draw()

    # h_muonPhi2.SetLineColor(r.kBlue)
    # h_muonPhi2.Draw("same")


    # can.SaveAs("plots/muon_phi_both.pdf")
    # can.SaveAs("plots/muon_phi_both.png")
    # can.Close()

    # drawHist(h_muonPhi1, "#phi", "plots/h_muonPhi1_test", kRed)
    # drawHist(h_muonPhi2, "#phi", "plots/h_muonPhi2_test", kblue)

    # # Draw the first histogram
    # h_muonPhi1 = drawHist(h_muonPhi1, "X-axis Label", "plots/h_muonPhi1_test", r.kBlue, "Histogram 1")

    # # Draw the second histogram on top of the first one
    # drawHist(h_muonPhi2, "X-axis Label", "save_path_2.png", r.kRed, "Histogram 2", "same")

    # drawCanvas(h_muonPhi1,h_muonPhi2, r'$0^\circ$ - $75^\circ$', r'$91^\circ$ - $180^\circ$', 'muonPhi', 0.75, 0.4, 0.95, 0.55)
    # drawCanvas(h_muonPt1,h_muonPt2, r'$0^\circ$ - $75^\circ$', r'$91^\circ$ - $180^\circ$', 'muonPt', 0.75, 0.4, 0.95, 0.55)
    # drawCanvas(h_muonEta1,h_muonEta2, r'$0^\circ$ - $75^\circ$', r'$91^\circ$ - $180^\circ$', 'muonEta', 0.75, 0.4, 0.95, 0.55)
    # drawCanvas(h_muon_freeInvBeta1,h_muon_freeInvBeta2, r'$0^\circ$ - $75^\circ$', r'$91^\circ$ - $180^\circ$', 'muonfreeInvBeta', 0.75, 0.4, 0.95, 0.55)
    # drawCanvas(h_muon_timeAtIpInOut1,h_muon_timeAtIpInOut2, r'$0^\circ$ - $75^\circ$', r'$91^\circ$ - $180^\circ$', 'muontimeAtIpInOut', 0.75, 0.4, 0.95, 0.55)
    # drawCanvas(h_muon_timeAtIpOutIn1,h_muon_timeAtIpOutIn2, r'$0^\circ$ - $75^\circ$', r'$91^\circ$ - $180^\circ$', 'muontimeAtIpOutIn', 0.75, 0.4, 0.95, 0.55)



    draw2DHist(h_muonGlobalYvsTime1, r'$0^\circ$ - $75^\circ$', r'$91^\circ$ - $180^\circ$', '2Dhisttest', 0.75, 0.4, 0.95, 0.55)

    # drawCanvas(h_muon1LegPhi1,h_muon1LegPhi2, r'$0^\circ$ - $75^\circ$', r'$91^\circ$ - $180^\circ$', 'muon1LegPhi', 0.75, 0.4, 0.95, 0.55)
    # drawCanvas(h_muon1LegPt1,h_muon1LegPt2, r'$0^\circ$ - $75^\circ$', r'$91^\circ$ - $180^\circ$', 'muon1LegPt', 0.75, 0.4, 0.95, 0.55)
    # drawCanvas(h_muon1LegEta1,h_muon1LegEta2, r'$0^\circ$ - $75^\circ$', r'$91^\circ$ - $180^\circ$', 'muon1LegEta', 0.75, 0.4, 0.95, 0.55)
    # drawCanvas(h_muon1Leg_freeInvBeta1,h_muon1Leg_freeInvBeta2, r'$0^\circ$ - $75^\circ$', r'$91^\circ$ - $180^\circ$', 'muon1LegfreeInvBeta', 0.75, 0.4, 0.95, 0.55)
    # drawCanvas(h_muon1Leg_timeAtIpInOut1,h_muon1Leg_timeAtIpInOut2, r'$0^\circ$ - $75^\circ$', r'$91^\circ$ - $180^\circ$', 'muon1LegtimeAtIpInOut', 0.75, 0.4, 0.95, 0.55)
    # drawCanvas(h_muon1Leg_timeAtIpOutIn1,h_muon1Leg_timeAtIpOutIn2, r'$0^\circ$ - $75^\circ$', r'$91^\circ$ - $180^\circ$', 'muon1LegtimeAtIpOutIn', 0.75, 0.4, 0.95, 0.55)

    # drawCanvas(h_splitmuonPhi1,h_splitmuonPhi2, r'$0^\circ$ - $75^\circ$', r'$91^\circ$ - $180^\circ$', 'splitmuonPhi', 0.75, 0.4, 0.95, 0.55)
    # drawCanvas(h_splitmuonPt1,h_splitmuonPt2, r'$0^\circ$ - $75^\circ$', r'$91^\circ$ - $180^\circ$', 'splitmuonPt', 0.75, 0.4, 0.95, 0.55)
    # drawCanvas(h_splitmuonEta1,h_splitmuonEta2, r'$0^\circ$ - $75^\circ$', r'$91^\circ$ - $180^\circ$', 'splitmuonEta', 0.75, 0.4, 0.95, 0.55)
    # drawCanvas(h_splitmuon_freeInvBeta1,h_splitmuon_freeInvBeta2, r'$0^\circ$ - $75^\circ$', r'$91^\circ$ - $180^\circ$', 'splitmuon_freeInvBeta', 0.75, 0.4, 0.95, 0.55)
    # drawCanvas(h_splitmuon_timeAtIpInOut1,h_splitmuon_timeAtIpInOut2, r'$0^\circ$ - $75^\circ$', r'$91^\circ$ - $180^\circ$', 'splitmuon_timeAtIpInOut', 0.75, 0.4, 0.95, 0.55)
    # drawCanvas(h_splitmuon_timeAtIpOutIn1,h_splitmuon_timeAtIpOutIn2, r'$0^\circ$ - $75^\circ$', r'$91^\circ$ - $180^\circ$', 'splitmuon_timeAtIpOutIn', 0.75, 0.4, 0.95, 0.55)

    # drawCanvas(h_lhcSTAMuonPhi1,h_lhcSTAMuonPhi2, r'$0^\circ$ - $75^\circ$', r'$91^\circ$ - $180^\circ$', 'lhcSTAMuonPhi', 0.75, 0.4, 0.95, 0.55)
    # drawCanvas(h_lhcSTAMuonPt1,h_lhcSTAMuonPt2, r'$0^\circ$ - $75^\circ$', r'$91^\circ$ - $180^\circ$', 'lhcSTAMuonPt', 0.75, 0.4, 0.95, 0.55)
    # drawCanvas(h_lhcSTAMuonEta1,h_lhcSTAMuonEta2, r'$0^\circ$ - $75^\circ$', r'$91^\circ$ - $180^\circ$', 'lhcSTAMuonEta', 0.75, 0.4, 0.95, 0.55)
 #   drawCanvas(h_lhcSTAMuon_freeInvBeta1,h_lhcSTAMuon_freeInvBeta2, r'$0^\circ$ - $75^\circ$', r'$91^\circ$ - $180^\circ$', 'lhcSTAMuon_freeInvBeta', 0.75, 0.4, 0.95, 0.55)
 #   drawCanvas(h_lhcSTAMuon_timeAtIpInOut1,h_lhcSTAMuon_timeAtIpInOut2, r'$0^\circ$ - $75^\circ$', r'$91^\circ$ - $180^\circ$', 'lhcSTAMuon_timeAtIpInOut', 0.75, 0.4, 0.95, 0.55)
 #   drawCanvas(h_lhcSTAMuon_timeAtIpOutIn1,h_lhcSTAMuon_timeAtIpOutIn2, r'$0^\circ$ - $75^\circ$', r'$91^\circ$ - $180^\circ$', 'lhcSTAMuon_timeAtIpOutIn', 0.75, 0.4, 0.95, 0.55)

    # Draw the first histogram with a label "Phi - 0 to 75"
    # h_muonPhi1 = drawHist(h_muonPhi1, "#phi", "plots/h_muonPhi1_test.png", r.kBlue, "Phi - 0 to 75", legend_label="Phi - 0 to 75")

    # # Draw the second histogram with a label "Phi - 91 to 180" and overlay it on top of the first one
    # h_muonPhi2 = drawHist(h_muonPhi2, "#phi", "plots/BothMuons.png", r.kRed, "Phi - 91 to 180", "same", legend_label="Phi - 91 to 180")

    # can2 = r.TCanvas()
    # can2.Draw()

    # h_muonPt1.SetLineColor(r.kRed)
    # h_muonPt1.Draw()
    # h_muonPt2.SetLineColor(r.kBlue)
    # h_muonPt2.Draw("same")

    # can2.SaveAs("plots/muon_pt_both.pdf")
    # can2.SaveAs("plots/muon_pt_both.png")





if __name__ == "__main__":
    main()