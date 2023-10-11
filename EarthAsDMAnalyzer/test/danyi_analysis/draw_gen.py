import argparse
import glob
from ROOT import TChain, TMath, TH1F, TCanvas, gStyle, TLegend, Math
import tdrStyle
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
    
    # load tree and assign process type
    prc = 'sig' if indir.find('91to180') >= 0 else 'bkg'
    tc = TChain(treename)
    tc.Add(indir)
    
    # create hist lists
    pt_list = []
    eta_list = []
    phi_list = []
    energy_list = []
    for i in range(4):
        pt_list.append(TH1F("pt_gen{}".format(i), "", 100, 0., 5000.))
        eta_list.append(TH1F("eta_gen{}".format(i), "", 50, 0., 3.))
        phi_list.append(TH1F("phi_gen{}".format(i), "", 50, -3.14, 3.14))
        energy_list.append(TH1F("energy_gen{}".format(i), "", 100, 0., 5000.))
    
    # event loop
    for evt in tc:
        gen_v = Math.PtEtaPhiMVector()
        for i in range(4):
            gen_v.SetCoordinates(evt.gen_pt[i], evt.gen_eta[i], evt.gen_phi[i], evt.gen_mass[i])
            pt_list[i].Fill(evt.gen_pt[i])
            eta_list[i].Fill(evt.gen_eta[i])
            phi_list[i].Fill(evt.gen_phi[i])
            energy_list[i].Fill(gen_v.E())
    
    # draw
    gStyle.SetOptStat(0)
    tdrStyle.SetTDRStyle()
    
    hist_dict = {
        'pt' : pt_list,
        'eta' : eta_list,
        'phi' : phi_list,
        'energy' : energy_list
    }
    
    xtitle_dict = {
        'pt' : 'p_{T} GeV',
        'eta' : '#eta',
        'phi' : '#phi',
        'energy' : 'E [GeV]'
    }
    
    for key, hist_list in hist_dict.items():
        can = TCanvas("mycan_{}".format(key), "", 800, 600)
        for i in range(4):
            if i == 0:
                hist_list[i].GetXaxis().SetTitle(xtitle_dict[key])
                hist_list[i].GetYaxis().SetTitle("A. U.")
            hist_list[i].SetLineColor(tdrStyle.colors['color_comp{}'.format(i+1)])
            hist_list[i].Scale(1./hist_list[i].Integral())
            hist_list[i].SetLineWidth(2)
        can.cd()
        for i, hist in enumerate(hist_list):
            if i == 0:
                # hist.SetMinimum(1e-4)
                hist.Draw("hist")
            else:
                hist.Draw("hist same")
                # if i >=1:
                #     break
        leg = TLegend(0.35, 0.75, 0.6, 0.94)
        for i, hist in enumerate(hist_list):
            leg.AddEntry(hist, "gen{}".format(i))
        leg.Draw("same")
        can.SaveAs(outdir + "/" + "gen_{}_{}.png".format(key,prc))
            
    

    
if __name__ == "__main__":
    main()