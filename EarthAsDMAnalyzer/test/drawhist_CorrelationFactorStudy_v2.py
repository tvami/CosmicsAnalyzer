import ROOT as r
import os
import array
import numpy as np
import argparse
import os.path
import matplotlib.pyplot as plt

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

def normalizeHistogram(hist):
    totalEntries = hist.GetEntries()
    if totalEntries != 0:
        scale = 1.0 / totalEntries
        hist.Scale(scale)


def drawCanvas(hist1, hist2, legendItem1, legendItem2, title,file,x1,x2,y1,y2):

    can = r.TCanvas()
    can.Draw()

    hist1.SetLineColor(r.kRed)
    hist1.Draw()

    hist2.SetLineColor(r.kBlue)

    integral_hist1 = hist1.Integral()
    integral_hist2 = hist2.Integral()

    # Scale hist2 by the ratio of integrals to normalize it to hist1
    if integral_hist1 != 0 and integral_hist2 != 0:
        scaling_factor = integral_hist1 / integral_hist2
        hist2.Scale(scaling_factor)


    hist2.Draw("same")

    legend = r.TLegend(x1,x2,y1,y2)  # Define the legend position
    legend.SetFillColor(0)
    legend.SetBorderSize(1)

    legend.AddEntry(hist1, legendItem1, "l")  # Add entry to the legend for the current histogram
    legend.AddEntry(hist2, legendItem2, "l")
    legend.Draw()

    # can.SaveAs(f"plots/{file}/{title}_both.pdf")
    can.SaveAs(f"plots/{file}/{title}_both.png")

    return can


def drawCanvasAndCreateRootFile(hist1, hist2, legendItem1, legendItem2, title,file, save, RootName, x1,x2,y1,y2):

    can = r.TCanvas()
    # can.Draw()

    output_file = r.TFile(f"{save}/{RootName}.root", "RECREATE")
    # print("This is output_file", output_file)


    hist1.SetLineColor(r.kRed)
    hist1.Draw()
    hist1.Write()

    hist2.SetLineColor(r.kBlue)

    integral_hist1 = hist1.Integral()
    integral_hist2 = hist2.Integral()

    # Scale hist2 by the ratio of integrals to normalize it to hist1
    if integral_hist1 != 0 and integral_hist2 != 0:
        scaling_factor = integral_hist1 / integral_hist2
        hist2.Scale(scaling_factor)


    hist2.Draw("same")
    hist2.Write()

    legend = r.TLegend(x1,x2,y1,y2)  # Define the legend position
    legend.SetFillColor(0)
    legend.SetBorderSize(1)

    legend.AddEntry(hist1, legendItem1, "l")  # Add entry to the legend for the current histogram
    legend.AddEntry(hist2, legendItem2, "l")
    legend.Draw()

    # can.SaveAs(f"plots/{file}/{title}_both.pdf")
    # can.SaveAs(f"{file}/{title}_both.png")

    output_file.Close()

    # print(f"save: {save}")
    # print(f"RootName: {RootName}")


    # print(f"Attempting to open file: plots/{save}/{title}_both.root")
    # try:
    #     output_file = r.TFile(f"plots/{save}/{title}_both.root", "RECREATE")
    #     print("File opened successfully.")
    # except OSError as e:
    #     print(f"Error: {e}")


    return can

def Create_RootFile_WithAllHistograms(input_directory, Output_Directory, Output_Name, histogram_name):

    # output_histogram = r.TCanvas()
    # # can.Draw()

    # for i in input_directory

    # output_file = r.TFile(f"{Output_Directory}/{Output_Name}.root", "RECREATE")

    os.makedirs(Output_Directory, exist_ok=True)

    # Create the output ROOT file
    output_file = r.TFile(f"{Output_Directory}/{Output_Name}.root", "RECREATE")

    # Iterate through files in the input_directory
    for root_file in os.listdir(input_directory):
        if root_file.endswith(".root"):
            # Open each ROOT file in the input_directory
            input_file = r.TFile(os.path.join(input_directory, root_file))

            # Assuming you have histograms named "hist" in your input files
            histogram = input_file.Get(histogram_name)

            if histogram:
                # Clone the histogram to the output file
                output_histogram = histogram.Clone()

                # Write the histogram to the output file
                output_histogram.Write()

            # Close the input file
            input_file.Close()

    # Close the output file
    output_file.Close()








def combine_ROOT_codes(ROOT_codes, input_directory, output_directory,output_name, histogram_name):
    # Check if the input list is empty
    if not ROOT_codes:
        return None

    print("This is ROOT_codes",ROOT_codes)

    for code in os.listdir(ROOT_codes):
        # Load the ROOT file with the given code (replace this with your actual loading logic)
        input_file = r.TFile(f"{input_directory}/{code}")

        if not input_file or input_file.IsZombie():
            print(f"Failed to load input file: {code}. Skipping.")
            continue


        # Assuming you have a histogram named "hist" in each input file (replace with your actual histogram)
        hist = input_file.Get(histogram_name)

        if not hist:
            print(f"Histogram not found in input file: {code}. Skipping.")
            input_file.Close()
            continue


        # Clone the histogram and add it to the output file
        hist_clone = hist.Clone()
        hist_clone.SetName(code)
        print("This is hist_clone", hist_clone)
        hist_clone.Write()
        
        # Close the input file
        input_file.Close()



    can = r.TCanvas()

    # Create and open the output ROOT file for writing

    output_file_path = f"{output_directory}/{output_name}.root"
    print(f"Attempting to create ROOT file: {output_file_path}")
    output_file = r.TFile(output_file_path, "RECREATE")

    # output_file = r.TFile(f"{output_directory}/{output_name}.root", "RECREATE")


    output_file.Close()


    
    return output_file







def draw2DHist(hist1, legendItem1, legendItem2, title): #x1,x2,y1,y2):

    can = r.TCanvas()
    can.Draw()
    hist1.SetLineColor(r.kRed)
    hist1.Draw("COLZ")

    # hist2.SetLineColor(r.kBlue)
    # hist2.Draw("same")

    legend = r.TLegend(x1,x2,y1,y2)  # Define the legend position
    legend.SetFillColor(0)
    legend.SetBorderSize(1)

    legend.AddEntry(hist1, legendItem1, "l")  # Add entry to the legend for the current histogram
    # legend.AddEntry(hist2, legendItem2, "l")
    legend.Draw()

    can.SaveAs(f"plots/Trigger1Leg/{title}_both.pdf")
    can.SaveAs(f"plots/Trigger1Leg/{title}_both.png")

    return can

def draw2DHistNoLegend(hist1, title, save): #x1,x2,y1,y2):

    can = r.TCanvas()
    can.Draw()
    hist1.SetLineColor(r.kRed)
    #hist1.SetStats(0);
    hist1.Draw("COLZ")

    output_file = r.TFile(f"plots/{save}/{title}_both.root", "RECREATE")
    hist1.Write()


    can.SaveAs(f"plots/{save}/{title}_both.pdf")
    can.SaveAs(f"plots/{save}/{title}_both.png")

    output_file.Close()

    return can

def testdraw2DHistNoLegend(hist1, title): #x1,x2,y1,y2):
    print("test")

def GetCorrelationPlot_old(hist1, title, type): #x1,x2,y1,y2):
    # Create a TProfile to display the correlation factor values
    profile = r.TProfile("Correlation Profile", "Correlation Factor vs Event Number_{type}", total_events, 0, total_events)

    # Fill the TProfile with the correlation factor values
    for i, factor in enumerate(correlation_factors):
        profile.Fill(i, factor)

    # Create a canvas to display the TProfile
    canvasCorrelationFactor = r.TCanvas("canvas", "Correlation Factor Profile", 800, 600)
    profile.Draw("P")  # "P" option indicates to draw points

    # Label axes and set the title
    profile.GetXaxis().SetTitle("Event Number")
    profile.GetYaxis().SetTitle("Correlation Factor")

    # Display the canvas
    canvasCorrelationFactor.Draw()
    canvasCorrelationFactor.SaveAs(f"plots/BugHuntCorrelationFactor_v6/Signal_CorrelationFactor.png")



def GetCorrelationPlot(name, list1, numberEvents, title, typeof, firstfolder, secondfolder, RootSave): #x1,x2,y1,y2):
    name = r.TCanvas("canvas", "Correlation Factor Profile", 800, 600)
    name.Draw()
    output_file = r.TFile(f"plots/{firstfolder}/{RootSave}/{title}.root", "RECREATE")

    graph = r.TGraph(len(numberEvents), array.array('d',numberEvents), array.array('d',list1))
    graph.SetTitle(f"Correlation Factor vs {title}")
    graph.GetXaxis().SetTitle("Event")
    graph.GetYaxis().SetTitle(f"{typeof} Correlation Factor")
    graph.Draw("AP")  # "APL" option indicates to draw points and lines

    graph.Write()

    # name.Draw()
    name.SaveAs(f"plots/{firstfolder}/{secondfolder}/{title}.png")

    output_file.Close()

def GetCorrelationPlotOverlay(hist1,hist2,hist3,hist4,hist5):
    name = r.TCanvas("canvas", "Correlation Factor Profile", 800, 600)
    name.Draw()
    output_file = r.TFile(f"plots/{firstfolder}/{RootSave}/{title}.root", "RECREATE")



    print("yes")

def overlay_root_files(input_directory, output_file_name,output_root_file_name):
    # Create a canvas to draw histograms
    # Create a canvas to draw histograms
    canvas = r.TCanvas("canvas", "Overlayed Histograms", 800, 600)
    
    # Create a legend to label histograms
    legend = r.TLegend(0.7, 0.7, 0.9, 0.9)
    
    # Create a list to hold histograms
    histograms = []

    # Create a list to hold graphs
    graphs = []
    
    for root_file in os.listdir(input_directory):
        if not root_file.endswith(".root"):
            continue  # Skip non-ROOT files
        
        file_name = os.path.join(input_directory, root_file)
        
        # Open the ROOT file
        file = r.TFile.Open(file_name)
        
        if not file or file.IsZombie():
            print(f"Failed to open file: {file_name}")
            continue
        
        # Assuming the histogram you want to overlay is named "hist"
        hist = file.Get("Graph")
        
        if not hist:
            print(f"Failed to get histogram from file: {file_name}")
            file.Close()
            continue
        
        # Set the line color, style, and marker style for the histogram
        hist.SetLineColor(len(histograms) + 1)
        hist.SetLineStyle(len(histograms) + 1)
        hist.SetMarkerStyle(20 + len(histograms))
        
        # Add the histogram to the legend
        legend.AddEntry(hist, root_file, "l")
        
        # Add the histogram to the list
        histograms.append(hist)
        
        file.Close()
    
    # Draw the histograms with option "same" to overlay them
    for i, hist in enumerate(histograms):
        if i == 0:
            hist.Draw()
        else:
            hist.Draw("same")
    
    # Draw the legend
    legend.Draw()

    for hist in histograms:
        hist.SetMinimum(-1)
        hist.SetMaximum(1)

    
    # Create a new TMultiGraph to hold and overlay the graphs
    multi_graph = r.TMultiGraph()
    
    for i in range(1, len(histograms)):
        multi_graph.Add(histograms[i])

    
    # Draw the overlayed graphs
    multi_graph.Draw("A")
    
    # Draw the legend
    legend.Draw()
    
    # Set the y-axis range for the overlayed graph
    multi_graph.SetMinimum(-1)
    multi_graph.SetMaximum(1)
    
    # Create a new ROOT file to save the overlayed graph
    output_root_file = r.TFile(output_root_file_name, "RECREATE")
    
    # Write the overlayed graph to the new ROOT file
    multi_graph.Write()
    
    # Close the new ROOT file
    output_root_file.Close()


def createAndAppendToRootFile(output_directory, output_file_name, histograms_to_append):

    # Create the full path for the output ROOT file
    # output_file_path = os.path.join(output_directory, output_file_name)

    os.makedirs(output_directory, exist_ok=True)

    # Create a new ROOT file in "RECREATE" mode
    # output_file = r.TFile(output_directory, "RECREATE")
    output_file = r.TFile(f"{output_directory}/{output_file_name}.root", "RECREATE")


    # Loop through the list of histograms to append
    for histogram in histograms_to_append:
        if histogram:
            # Clone the histogram to the output file
            output_histogram = histogram.Clone()

            # Write the histogram to the output file
            output_histogram.Write()

    # Close the output ROOT file
    output_file.Close()

    return output_file


def DoesItemPassTrigger(hist1, title, MuonType): #x1,x2,y1,y2):

    can = r.TCanvas()
    can.Draw()
    hist1.SetLineColor(r.kRed)
    #hist1.SetStats(0);
    hist1.Draw("COLZ")

    return can


def create_and_save_tprofile(hist2D, title, x_axis_title, y_axis_title):
    # Creating a TProfile
    num_x_bins = hist2D.GetNbinsX()
    x_min = hist2D.GetXaxis().GetXmin()
    x_max = hist2D.GetXaxis().GetXmax()
    can2 = r.TProfile(f"profile_{title}", title, num_x_bins, x_min, x_max)

    # Setting axis titles
    can2.GetXaxis().SetTitle(x_axis_title)
    can2.GetYaxis().SetTitle(y_axis_title)

    # Loop through the bins and fill the TProfile with mean values
    for i in range(1, num_x_bins + 1):
        bin_content_0 = hist2D.GetBinContent(i, 1)
        bin_content_1 = hist2D.GetBinContent(i, 2)
        mean_value = (bin_content_0 + bin_content_1) / 2.0
        can2.Fill(hist2D.GetXaxis().GetBinCenter(i), mean_value)

    # Drawing the TProfile
    can2.Draw()

    # Saving the TProfile as PDF and PNG files
    can2.SaveAs(f"plots/test/{title}_both.pdf")
    can2.SaveAs(f"plots/test/{title}_both.png")

    return can2




def main():

    # Load first tree
    #tree1 = loadTree("MCntuple0to75-10000Events.root", "muonPhiAnalyzer/tree")
    #tree1 = loadTree("MCntuple-0to75-splitMuons-3Triggers-WithEnergy.root", "muonPhiAnalyzer/tree")
    

    # # path = 
    # ###Background path
    path = "/ceph/cms/store/user/lbrennan/EarthAsDM/CMSSW_13_0_10/src/Cosmics/Cosmics/crab_Ntuplizer-0to75Theta-4to3000GeV-126X_mcRun3_2022cosmics_realistic_deco_v1_v4/230907_052916/0000/Ntuplizer-0to75Theta-4to3000GeV_v1_40.root"
    tree1 = loadTree(os.path.abspath(path), "muonPhiAnalyzer/tree")
        #("/ceph/cms/store/user/lbrennan/EarthAsDM/Cosmics/crab_RAWtoReco-91to180Theta-3000to4000GeV-126X_mcRun3_2022cosmics_realistic_deco_v1_v4/230811_215442/0000/3RR-91to180Theta-3000to4000GeV_Total.root"), "muonPhiAnalyzer/tree")
        #/home/users/lbrennan/work/MCNtuplizer/CMSSW_13_0_10/src/CosmicsAnalyzer/EarthAsDMAnalyzer/test/Ntuplizer-91to180Theta-3000to4000GeV_OnlyTheFirstEvent_v7.root"), "muonPhiAnalyzer/tree")
        #"/ceph/cms/store/user/lbrennan/EarthAsDM/CMSSW_13_0_10/src/Cosmics/Cosmics/crab_Ntuplizer-0to75Theta-4to3000GeV-126X_mcRun3_2022cosmics_realistic_deco_v1_v4/230907_052916/0000/Ntuplizer-0to75Theta-4to3000GeV_v1_10.root"), "muonPhiAnalyzer/tree")
    
    nevents1 = tree1.GetEntries()
    print("nentries1 = ", nevents1)
    
    # # Load second tree
    # ###91 to 180
    # ###Signal
    # tree2 = loadTree(os.path.abspath("/ceph/cms/store/user/lbrennan/EarthAsDM/CMSSW_13_0_10/src/Cosmics/Cosmics/crab_Ntuplizer-91to180Theta-3000to4000GeV-126X_mcRun3_2022cosmics_realistic_deco_v1_v6/230901_183652/0000/Ntuplizer-91to180Theta-3000to4000GeV_v2_10.root"), "muonPhiAnalyzer/tree")
    # #(os.path.abspath("/ceph/cms/store/user/lbrennan/EarthAsDM/CMSSW_13_0_10/src/Cosmics/Cosmics/crab_Ntuplizer-91to180Theta-3000to4000GeV-126X_mcRun3_2022cosmics_realistic_deco_v1_v6/230901_183652/0000/Ntuplizer-91to180Theta-3000to4000GeV_v2_10.root"), "muonPhiAnalyzer/tree")



    # # tree2 = loadTree("MCntuple-91to180-splitMuons-3Triggers-WithEnergy.root", "muonPhiAnalyzer/tree") 
    # nevents2 = tree2.GetEntries()


    # print("nentries2 = ", nevents2)


    h_muonPhi2 = r.TH1F("h_muonPhi2", "Muon #phi - 91 to 180;#phi;yield", 100, -4, 4)
    h_muonPt2 = r.TH1F("h_muonPt2", "Muon pT Comparison;pt;yield", 100, 0, 300)
    h_muonEta2 = r.TH1F("h_muonEta2", "Muon #eta Comparison;#eta;yield", 100, -2.5, 2.5)

    h_muon_freeInvBeta1 = r.TH1F("h_muon_freeInvBeta1", "Muon FreeInverseBeta Comparison;freeInvBeta;yield", 100, -50, 50)
    h_muon_freeInvBeta2 = r.TH1F("h_muon_freeInvBeta2", "Muon FreeInverseBeta Comparison;freeInvBeta;yield", 100, -50, 50)
    h_muon_timeAtIpInOut1 = r.TH1F("h_muon_timeAtIpInOut1", "Muon TimeAtIpInOut Comparison;timeAtIpInOut;yield", 100, -500, 500)
    h_muon_timeAtIpInOut2 = r.TH1F("h_muon_timeAtIpInOut2", "Muon TimeAtIpInOut Comparison;timeAtIpInOut;yield", 100, -500, 500)
    h_muon_timeAtIpOutIn1 = r.TH1F("h_muon_timeAtIpOutIn1", "Muon TimeAtIpOutIn Comparison;timeAtIpOutIn;yield", 100, -500, 500)
    h_muon_timeAtIpOutIn2 = r.TH1F("h_muon_timeAtIpOutIn2", "Muon TimeAtIpOutIn Comparison;timeAtIpOutIn;yield", 100, -500, 500)

    h_muonEnergy1 = r.TH1F("h_muonEnergy1", "Muon Energy Comparison;energy;yield", 100, -100, 3000)
    h_muonEnergy2 = r.TH1F("h_muonEnergy2", "Muon Energy Comparison;energy;yield", 100, -100, 3000)

    h_muonGlobalYvsTime1 = r.TH2D("h_muonGlobalYvsTime1", "Background MC GlobalY vs time;time[ns];GlobalY", 100, -100, 100, 50, -800, 800)
    h_muonGlobalYvsTime2 = r.TH2D("h_muonGlobalYvsTime2", "Signal MC GlobalY vs time;time[ns];GlobalY", 100, -100, 100, 50, -800, 800)

    h_muonGlobalYvsTime1CorrelationFactor = r.TH1F("h_muonGlobalYvsTime1CorrelationFactor", "Muon correlationFactor Comparison;#Event;correlationFactor", 100, -1, 1)
    #h_muonGlobalYvsTime1CorrelationFactor2 = r.TH1F

    # h_muonGlobalYvsTime1CorrelationFactor = r.TH1F("h_muonGlobalYvsTime1CorrelationFactor", "Muon correlationFactor Comparison;#Event;correlationFactor", 100, -1, 1)
    h_muon_n = r.TH1F("h_muon_n", "Muon #Phi Comparison;#phi;yield", 100, -4, 4)

    h_TriggervsEnergy1 = r.TH2D("h_TriggervsEnergy1", "Background MC Cosmic Muon trigger vs energy;0 or 1;energy", 2, -0.5, 1.5, 100, 0, 1000)
    h_TriggervsEnergy2 = r.TH2D("h_TriggervsEnergy2", "Signal MC Cosmic Muon trigger vs energy;0 or 1;energy", 2, -0.5, 1.5, 101, 0.0, 1000.0)


#the ndof and dt segment plots
    # h_TriggervsEnergy2.SetOption("overflow")


    h_muonPhi_PositiveCorrelationFactor1 = r.TH1F("h_muonPhi_PositiveCorrelationFactor1", "Muon #Phi Comparison;#phi;yield", 50, -4, 4)
    h_muonPhi_NegativeCorrelationFactor1 = r.TH1F("h_muonPhi_NegativeCorrelationFactor1", "Muon #Phi Comparison;#phi;yield", 50, -4, 4)    

    h_muonEta_PositiveCorrelationFactor1 = r.TH1F("h_muonEta_PositiveCorrelationFactor1", "Muon #eta Comparison;#eta;yield", 50, -2.5, 2.5)
    h_muonEta_NegativeCorrelationFactor1 = r.TH1F("h_muonEta_NegativeCorrelationFactor1", "Muon #eta Comparison;#eta;yield", 50, -2.5, 2.5)


    h_muonPt_PositiveCorrelationFactor1 = r.TH1F("h_muonPt_PositiveCorrelationFactor1", "Muon pT Comparison;pt;yield", 50, 0, 300)
    h_muonPt_NegativeCorrelationFactor1 = r.TH1F("h_muonPt_NegativeCorrelationFactor1", "Muon pT Comparison;pt;yield", 50, 0, 300)



    h_muon_comb_ndof_PositiveCorrelationFactor1 = r.TH1F("h_muon_ndof_PositiveCorrelationFactor1", "Muon ndof Comparison;ndof;yield", 50, 0, 40)
    h_muon_comb_ndof_NegativeCorrelationFactor1 = r.TH1F("h_muon_ndof_NegativeCorrelationFactor1", "Muon ndof Comparison;ndof;yield", 50, 0, 40)




    h_muon_dtSeg_globY_PositiveCorrelationFactor1 = r.TH1F("h_muon_dtSeg_globY_PositiveCorrelationFactor1", "Muon_dtSeg_globY_Comparison;dtSeg_globY;yield", 50, -1000, 1000)
    h_muon_dtSeg_globY_NegativeCorrelationFactor1 = r.TH1F("h_muon_dtSeg_globY_NegativeCorrelationFactor1", "Muon_dtSeg_globY Comparison;dtSeg_globY;yield", 50, -1000, 1000)

    # h_muon_dtSeg_globY_PositiveCorrelationFactor2 = r.TH1F("h_muon_dtSeg_globY_PositiveCorrelationFactor2", "Muon dtSeg_globY Comparison;dtSeg_globY;yield", 50, -2, 2)
    # h_muon_dtSeg_globY_NegativeCorrelationFactor2 = r.TH1F("h_muon_dtSeg_globY_NegativeCorrelationFactor2", "Muon dtSeg_globY Comparison;dtSeg_globY;yield", 50, -2, 2)


    h_muon_dtSeg_globZ_PositiveCorrelationFactor1 = r.TH1F("h_muon_dtSeg_globZ_PositiveCorrelationFactor1", "Muon dtSeg_globZ Comparison;dtSeg_globY;yield", 50, -1000, 1000)
    h_muon_dtSeg_globZ_NegativeCorrelationFactor1 = r.TH1F("h_muon_dtSeg_globZ_NegativeCorrelationFactor1", "Muon dtSeg_globZ Comparison;dtSeg_globY;yield", 50, -1000, 1000)

    # h_muon_dtSeg_globZ_PositiveCorrelationFactor2 = r.TH1F("h_muon_dtSeg_globZ_PositiveCorrelationFactor2", "Muon dtSeg_globZ Comparison;dtSeg_globY;yield", 100, -2, 2)
    # h_muon_dtSeg_globZ_NegativeCorrelationFactor2 = r.TH1F("h_muon_dtSeg_globZ_NegativeCorrelationFactor2", "Muon dtSeg_globZ Comparison;dtSeg_globY;yield", 100, -2, 2)



    h_TriggervsPhi1 = r.TH2D("h_TriggervsPhi1", "Background MC Cosmic Muon trigger vs Phi;0 or 1;Phi", 2, -0.5, 1.5, 100, -4, 4)
    h_TriggervsPhi2 = r.TH2D("h_TriggervsPhi2", "Signal MC Cosmic Muon trigger vs Phi;0 or 1;Phi", 2, -0.5, 1.5, 100, -4, 4)

    h_TriggervsPt1 = r.TH2D("h_TriggervsPt1", "Background MC Cosmic Muon vs pT;0 or 1;pT", 2, -0.5, 1.5, 101, 0, 1000)
    h_TriggervsPt2 = r.TH2D("h_TriggervsPt2", "Signal MC Cosmic Muon trigger vs pT;0 or 1;pT", 2, -0.5, 1.5, 101, 0, 2000) #100, 0, 2, 2, 0, 2000)

    h_TriggervsEta1 = r.TH2D("h_TriggervsEta1", "Background MC Cosmic Muon vs Eta;0 or 1;Eta", 2, -0.5, 1.5, 100, -4, 4)
    h_TriggervsEta2 = r.TH2D("h_TriggervsEta2", "Signal MC Cosmic Muon trigger vs Eta;0 or 1;Eta", 2, -0.5, 1.5, 100, -4, 4)

    h_GenVsRecoTracks1 = r.TH2D("h_GenVsRecoTracks1", "Signal MC GenVsReco Tracks;0 or 1;Eta", 2, -0.5, 1.5, 100, 0, 1000)
    h_GenVsRecoTracks2 = r.TH2D("h_GenVsRecoTracks2", "Signal MC trigger vs Eta;0 or 1;Eta", 100, 0, 2, 2, 0, 2000)


    event_count = 0
    pfreq = 1000
    accepted_printed = False
    rejected_printed = False
    correlation_factors = []
    Event_Number = []
    total_events = 0
    # graph = r.TGraph()
    muon_n = []
    muon_phi_test = []
    muon_dtSeg_globY = []
    muon_dtSeg_n =[]


    muon_phi_average =0
    muon_dtSeg_globY_average = 0
    muon_dtSeg_n_average = 0


    # for tree in [tree1, tree2]:

    for tree in [tree1]:
        for event in tree:

            total_events += 1
            Event_Number.append(total_events)
            muon_n.append(event.muon_n)
            for (dtSeg_t0timing, dtSeg_globY, dtSeg_found) in zip (event.muon_dtSeg_t0timing, event.muon_dtSeg_globY, event.muon_dtSeg_found):
                if dtSeg_t0timing > -900 and dtSeg_t0timing < 900 and dtSeg_t0timing != 0 and dtSeg_found == 1: #event.muon_tofMap_found == 1:
                    if tree == tree1:
                        h_muonGlobalYvsTime1.Fill(dtSeg_t0timing, dtSeg_globY)
                        correlationFactor = h_muonGlobalYvsTime1.GetCorrelationFactor()
                        h_muonGlobalYvsTime1CorrelationFactor.Fill(correlationFactor)
                        correlation_factors.append(correlationFactor)
                        if correlationFactor > 0:
                            for phi in event.muon_phi:
                                lengthPhi = len(event.muon_phi)
                                muon_phi_average += phi
                                h_muonPhi_PositiveCorrelationFactor1.Fill(phi)
                            for pt in event.muon_pt:
                                h_muonPt_PositiveCorrelationFactor1.Fill(pt)   
                            for eta in event.muon_eta:
                                h_muonEta_PositiveCorrelationFactor1.Fill(eta)
                            for ndof in event.muon_comb_ndof:
                                lengthPhi = len(event.muon_phi)
                                muon_phi_average += phi
                                h_muon_comb_ndof_PositiveCorrelationFactor1.Fill(ndof)
                            for dtSeg_globY in event.muon_dtSeg_globY:
                                if dtSeg_globY > 10 or dtSeg_globY < -10:
                                    h_muon_dtSeg_globY_PositiveCorrelationFactor1.Fill(dtSeg_globY)
                            for dtSeg_globZ in event.muon_dtSeg_globZ:
                                if dtSeg_globZ > 10 or dtSeg_globZ < -10:
                                    h_muon_dtSeg_globZ_PositiveCorrelationFactor1.Fill(dtSeg_globZ)


                        if correlationFactor <= 0:
                            for phi in event.muon_phi:
                                lengthPhi = len(event.muon_phi)
                                muon_phi_average += phi
                                h_muonPhi_NegativeCorrelationFactor1.Fill(phi)
                            for pt in event.muon_pt:
                                h_muonPt_NegativeCorrelationFactor1.Fill(pt)

                            for eta in event.muon_eta:
                                h_muonEta_NegativeCorrelationFactor1.Fill(eta)

                            for ndof in event.muon_comb_ndof:
                                lengthPhi = len(event.muon_phi)
                                muon_phi_average += phi
                                h_muon_comb_ndof_NegativeCorrelationFactor1.Fill(ndof)
                            for dtSeg_globY in event.muon_dtSeg_globY:
                                if dtSeg_globY > 10 or dtSeg_globY < -10:
                                    h_muon_dtSeg_globY_NegativeCorrelationFactor1.Fill(dtSeg_globY)

                            for dtSeg_globZ in event.muon_dtSeg_globZ:
                                if dtSeg_globZ > 10 or dtSeg_globZ < -10:
                                    h_muon_dtSeg_globZ_NegativeCorrelationFactor1.Fill(dtSeg_globZ)


            event_count += 1


    RootFilesPath = 'plots/BackgroundPlots/test_CorrelationFactorVsBackground_v1/RootFiles'
    pngFilePath = 'plots/BackgroundPlots/test_CorrelationFactorVsBackground_v1/pngs'

    #histogram_name = "h_muon_dtSeg_globY_PositiveCorrelationFactor1"
    # RootSavedPath_inputDirectory = 'plots/BackgroundPlots/CorrelationFactorvsBackgroundEventParameters_RootFileCreationTest_v1/RootFiles'
    RootSavedPath_outputDirectory = 'plots/BackgroundPlots/test_CorrelationFactorVsBackground_v1/ProofItWorks'

    histograms_to_append = [h_muon_dtSeg_globY_PositiveCorrelationFactor1, h_muon_dtSeg_globZ_PositiveCorrelationFactor1, h_muon_comb_ndof_PositiveCorrelationFactor1, h_muonPhi_PositiveCorrelationFactor1, h_muonPt_PositiveCorrelationFactor1, h_muonEta_PositiveCorrelationFactor1,
    h_muonPhi_NegativeCorrelationFactor1]
    createAndAppendToRootFile(RootSavedPath_outputDirectory, 'test2', histograms_to_append)




if __name__ == "__main__":
    main()