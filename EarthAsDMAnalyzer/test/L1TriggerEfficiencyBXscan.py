#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BX / quality scan of the per-hemisphere L1 DT trigger efficiency, data vs MC.

Companion to L1TriggerEfficiency.py. It runs the same tag-and-probe but, in a
single pass over each tree, stores for every probed leg the list of (|bx|, quality)
of the DT trigger primitives matched to it. The full (|bx|<=cut, quality>=cut) grid
is then evaluated without re-looping, and the inclusive efficiencies and the
data/MC scale factor (SF = eff_data / eff_MC) are tabulated.

This separates the two sources of any data/MC difference:
  - timing      : how the efficiency depends on the allowed BX window
  - intrinsic   : how it depends on the DT trigger quality code

It also writes a plot of efficiency vs BX window (data vs MC) at a fixed quality.

Usage:
  python3 L1TriggerEfficiencyBXscan.py -d Ntuplizer-Data-Run2025A.root \
                                       -m Ntuplizer-MC-CosmicToMu.root \
                                       -n L1DTTrigEff_BXscan
"""
import argparse, math, os
import ROOT

def dphi(a, b):
    d = a - b
    while d > math.pi:  d -= 2 * math.pi
    while d < -math.pi: d += 2 * math.pi
    return d

def dR(e1, p1, e2, p2):
    return math.hypot(e1 - e2, dphi(p1, p2))

def circ_mean(angles):
    s = sum(math.sin(a) for a in angles)
    c = sum(math.cos(a) for a in angles)
    return math.atan2(s, c)

def collect(path, tree_name, drcut, minseg):
    """One pass: per both-hemisphere event, store (up_prims, low_prims) where each
    is the list of (abs(bx), quality) of DT primitives matched (dR<cut) to the leg."""
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise SystemExit("Cannot open " + path)
    t = f.Get(tree_name)
    if not t:
        raise SystemExit("No tree '%s' in %s" % (tree_name, path))
    evts = []
    for ev in t:
        if not any(0 < p < 1e4 for p in ev.muon_tuneP_Pt):
            continue
        up, low = [], []
        for j in range(len(ev.muon_dtSeg_globY)):
            y = ev.muon_dtSeg_globY[j]
            if abs(y) > 9000:
                continue
            (up if y > 0 else low).append((ev.muon_dtSeg_eta[j], ev.muon_dtSeg_phi[j]))
        if len(up) < minseg or len(low) < minseg:
            continue
        uleg = (sum(e for e, p in up) / len(up),  circ_mean([p for e, p in up]))
        lleg = (sum(e for e, p in low) / len(low), circ_mean([p for e, p in low]))

        def matched(leg, want_up):
            le, lp = leg
            out = []
            for k in range(ev.dtTrigPh_n):
                y = ev.dtTrigPh_globY[k]
                if abs(y) > 9000 or (y > 0) != want_up:
                    continue
                if dR(le, lp, ev.dtTrigPh_globEta[k], ev.dtTrigPh_globPhi[k]) < drcut:
                    out.append((abs(ev.dtTrigPh_bx[k]), ev.dtTrigPh_quality[k]))
            return out

        evts.append((matched(uleg, True), matched(lleg, False)))
    f.Close()
    return evts

def efficiency(evts, bxcut, qcut):
    tag = pas = 0
    for up, low in evts:
        ufire = any(b <= bxcut and q >= qcut for b, q in up)
        lfire = any(b <= bxcut and q >= qcut for b, q in low)
        if lfire:
            tag += 1; pas += 1 if ufire else 0
        if ufire:
            tag += 1; pas += 1 if lfire else 0
    return (pas / tag if tag else 0.0)

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("-d", "--data", default=os.path.join(here, "Ntuplizer-Data-Run2025A.root"))
    ap.add_argument("-m", "--mc",   default=os.path.join(here, "Ntuplizer-MC-CosmicToMu.root"))
    ap.add_argument("-t", "--tree", default="muonPhiAnalyzer/tree")
    ap.add_argument("-o", "--outdir", default=here)
    ap.add_argument("-n", "--name", default="L1DTTrigEff_BXscan")
    ap.add_argument("--dr", type=float, default=0.4)
    ap.add_argument("--minseg", type=int, default=1)
    ap.add_argument("--bxlist",   default="0,1,2,3,4", help="BX windows to scan")
    ap.add_argument("--quallist", default="0,2,3,4,5,6", help="quality thresholds to scan")
    ap.add_argument("--plotqual", type=int, default=4, help="quality for the eff-vs-BX plot")
    args = ap.parse_args()

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)

    bxs   = [int(x) for x in args.bxlist.split(",")]
    quals = [int(x) for x in args.quallist.split(",")]

    D = collect(args.data, args.tree, args.dr, args.minseg)
    M = collect(args.mc,   args.tree, args.dr, args.minseg)
    print("both-hemisphere events: data=%d  mc=%d\n" % (len(D), len(M)))

    for bxcut in bxs:
        print("================ |bx| <= %d ================" % bxcut)
        print("%7s %10s %10s %12s" % ("qual>=", "eff_data", "eff_MC", "SF=data/MC"))
        for qcut in quals:
            ed = efficiency(D, bxcut, qcut)
            em = efficiency(M, bxcut, qcut)
            print("%7d %10.3f %10.3f %12.3f" % (qcut, ed, em, (ed / em if em else 0)))
        print()

    # plot: efficiency vs BX window, data vs MC, at fixed quality
    gd, gm = ROOT.TGraph(), ROOT.TGraph()
    for i, b in enumerate(bxs):
        gd.SetPoint(i, b, efficiency(D, b, args.plotqual))
        gm.SetPoint(i, b, efficiency(M, b, args.plotqual))
    c = ROOT.TCanvas("c", "c", 800, 600); c.SetGridy(True)
    gd.SetTitle("L1 DT efficiency vs BX window (quality #geq %d)"
                ";max |BX| of DT primitive;per-hemisphere L1 DT efficiency" % args.plotqual)
    gd.SetMarkerStyle(20); gd.SetMarkerColor(ROOT.kBlack); gd.SetLineColor(ROOT.kBlack); gd.SetLineWidth(2)
    gm.SetMarkerStyle(24); gm.SetMarkerColor(ROOT.kRed);   gm.SetLineColor(ROOT.kRed);   gm.SetLineWidth(2)
    gd.GetYaxis().SetRangeUser(0, 1.05); gd.GetXaxis().SetLimits(bxs[0] - 0.3, bxs[-1] + 0.3)
    gd.Draw("ALP"); gm.Draw("LP SAME")
    leg = ROOT.TLegend(0.45, 0.20, 0.88, 0.38); leg.SetBorderSize(0); leg.SetFillStyle(0)
    leg.AddEntry(gd, "Data", "lp"); leg.AddEntry(gm, "MC", "lp"); leg.Draw()
    base = os.path.join(args.outdir, args.name)
    c.SaveAs(base + ".png"); c.SaveAs(base + ".pdf")
    print("Wrote", base + ".png/.pdf")

if __name__ == "__main__":
    main()
