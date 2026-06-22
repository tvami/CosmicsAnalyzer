#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Per-hemisphere L1 (DT Local Trigger) efficiency on cosmic muons.

This reproduces the tag-and-probe method of the CMS cosmic commissioning paper
(JINST 5 (2010) T03002, Fig. 12) using the dtTrigPh_* branches of the v5
ntuplizer.

Idea:
  - Select a cosmic muon that crosses BOTH hemispheres (offline DT segments with
    global Y > 0 and < 0). This is an unbiased reference: a real through-going
    muon traversed the top and the bottom of the detector.
  - TAG  : the DT Local Trigger fired in one hemisphere (>=1 in-time, good-quality
           DT trigger primitive matched to the offline leg in that hemisphere).
  - PROBE: ask whether the DT Local Trigger also fired in the OTHER hemisphere.
  - Efficiency( probe pT ) = N(tag & probe) / N(tag).

Because a through-going cosmic has a single momentum (both legs measure the same
pT), the efficiency is binned in the offline muon pT. The measurement is filled
symmetrically (tag=lower/probe=upper AND tag=upper/probe=lower).

--acceptance applies fiducial cuts inspired by the paper (offline pT > 5 GeV,
>= 2 DT stations on the probe leg, probe leg away from DT sector phi-cracks),
which lift the plateau the same way the paper's "with acceptance cuts" curve does.

Usage:
  # single sample
  python3 L1TriggerEfficiency.py -i ntuple.root -n L1DTTriggerEfficiency
  # data vs MC overlay
  python3 L1TriggerEfficiency.py -i data.root --mc mc.root -n L1DTTriggerEfficiency_DataMC
  # with acceptance cuts
  python3 L1TriggerEfficiency.py -i data.root --mc mc.root --acceptance -n ..._acc
"""
import argparse, math, array, os
import ROOT

# ---------------------------------------------------------------------------
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

PTBINS = array.array('d', [2, 3, 5, 7, 10, 15, 20, 25, 30, 40, 50,
                           70, 100, 150, 200, 300, 500, 1000])

# ---------------------------------------------------------------------------
def passes_acceptance(leg_phi, n_stations, pt):
    """Fiducial cuts approximating the paper's 'with acceptance cuts' series."""
    if pt < 5.0:
        return False
    if n_stations < 2:
        return False
    # DT sectors are 30 deg wide; reject probe legs within 5 deg of a sector edge.
    pos = (math.degrees(leg_phi)) % 30.0      # 0..30, sector centre at 0/30, edge at 15
    if abs(pos - 15.0) < 5.0:
        return False
    return True

# ---------------------------------------------------------------------------
def compute(tree, qual, bxmax, drcut, minseg, acceptance):
    """Run the tag-and-probe over a tree; return (eff_comb, eff_up, eff_low, stats)."""
    suffix = str(id(tree))
    eff_comb = ROOT.TEfficiency("eff_comb_" + suffix, "", len(PTBINS) - 1, PTBINS)
    eff_up   = ROOT.TEfficiency("eff_up_"   + suffix, "", len(PTBINS) - 1, PTBINS)
    eff_low  = ROOT.TEfficiency("eff_low_"  + suffix, "", len(PTBINS) - 1, PTBINS)
    for e in (eff_comb, eff_up, eff_low):
        e.SetStatisticOption(ROOT.TEfficiency.kFCP)  # Clopper-Pearson

    n_total = n_both = n_tag = n_pass = 0

    for ev in tree:
        n_total += 1
        pts = [p for p in ev.muon_tuneP_Pt if 0 < p < 1e4]
        if not pts:
            continue
        pt = max(pts)

        # offline legs from DT segments (position-based hemispheres)
        up, low = [], []   # lists of (eta, phi, station)
        for j in range(len(ev.muon_dtSeg_globY)):
            y = ev.muon_dtSeg_globY[j]
            if abs(y) > 9000:
                continue
            rec = (ev.muon_dtSeg_eta[j], ev.muon_dtSeg_phi[j], ev.muon_dtSeg_Station_[j])
            (up if y > 0 else low).append(rec)

        if len(up) < minseg or len(low) < minseg:
            continue
        n_both += 1

        upper_leg = (sum(r[0] for r in up) / len(up),  circ_mean([r[1] for r in up]))
        lower_leg = (sum(r[0] for r in low) / len(low), circ_mean([r[1] for r in low]))
        up_nst  = len(set(r[2] for r in up))
        low_nst = len(set(r[2] for r in low))

        def fired(leg, want_upper):
            le, lp = leg
            for k in range(ev.dtTrigPh_n):
                y = ev.dtTrigPh_globY[k]
                if abs(y) > 9000 or (y > 0) != want_upper:
                    continue
                if abs(ev.dtTrigPh_bx[k]) > bxmax or ev.dtTrigPh_quality[k] < qual:
                    continue
                if dR(le, lp, ev.dtTrigPh_globEta[k], ev.dtTrigPh_globPhi[k]) < drcut:
                    return True
            return False

        up_fired  = fired(upper_leg, True)
        low_fired = fired(lower_leg, False)

        # tag = lower, probe = upper
        if low_fired and (not acceptance or passes_acceptance(upper_leg[1], up_nst, pt)):
            n_tag += 1
            eff_up.Fill(up_fired, pt); eff_comb.Fill(up_fired, pt)
            if up_fired: n_pass += 1
        # tag = upper, probe = lower
        if up_fired and (not acceptance or passes_acceptance(lower_leg[1], low_nst, pt)):
            n_tag += 1
            eff_low.Fill(low_fired, pt); eff_comb.Fill(low_fired, pt)
            if low_fired: n_pass += 1

    return eff_comb, eff_up, eff_low, (n_total, n_both, n_tag, n_pass)

# ---------------------------------------------------------------------------
def report(label, stats, qual, bxmax, drcut, minseg, acceptance):
    n_total, n_both, n_tag, n_pass = stats
    print("=" * 64)
    print("Sample                               :", label)
    print("Total events                         :", n_total)
    print("Cosmics crossing both hemispheres    :", n_both)
    print("Tag-and-probe pairs (denominator)    :", n_tag)
    print("Pairs where probe also fired (numer.):", n_pass)
    if n_tag:
        p  = n_pass / n_tag
        lo = ROOT.TEfficiency.ClopperPearson(n_tag, n_pass, 0.6827, False)
        hi = ROOT.TEfficiency.ClopperPearson(n_tag, n_pass, 0.6827, True)
        print("Inclusive per-hemisphere efficiency  : %.4f  [%.4f, %.4f]" % (p, lo, hi))
    print("Cuts: quality>=%d |bx|<=%d dR<%.2f minseg=%d acceptance=%s"
          % (qual, bxmax, drcut, minseg, acceptance))
    print("=" * 64)

def style(eff, color, marker):
    eff.SetLineColor(color); eff.SetMarkerColor(color); eff.SetMarkerStyle(marker)
    eff.SetMarkerSize(1.0)

def frame_and_draw(first, ytitle="L1 DT efficiency"):
    first.Draw("AP")
    ROOT.gPad.Update()
    g = first.GetPaintedGraph()
    if g:
        g.GetYaxis().SetRangeUser(0.0, 1.05)
        g.GetXaxis().SetLimits(PTBINS[0], PTBINS[-1])
    ROOT.gPad.Update()

# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("-i", "--input",  default=os.path.join(here, "Ntuplizer-Cosmics-Data.root"))
    ap.add_argument("--mc",           default=None, help="second (MC) ntuple to overlay")
    ap.add_argument("-t", "--tree",   default="muonPhiAnalyzer/tree")
    ap.add_argument("-o", "--outdir", default=here)
    ap.add_argument("-n", "--name",   default="L1DTTriggerEfficiency")
    ap.add_argument("--qual", type=int,   default=4)
    ap.add_argument("--bx",   type=int,   default=1)
    ap.add_argument("--dr",   type=float, default=0.4)
    ap.add_argument("--minseg", type=int, default=1)
    ap.add_argument("--acceptance", action="store_true",
                    help="apply paper-like fiducial cuts (pT>5, >=2 stations, away from phi-cracks)")
    args = ap.parse_args()

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)

    def load(path):
        f = ROOT.TFile.Open(path)
        if not f or f.IsZombie():
            raise SystemExit("Cannot open " + path)
        t = f.Get(args.tree)
        if not t:
            raise SystemExit("No tree '%s' in %s" % (args.tree, path))
        return f, t

    fdat, tdat = load(args.input)
    eff_c, eff_u, eff_l, st = compute(tdat, args.qual, args.bx, args.dr, args.minseg, args.acceptance)
    lbl0 = ("DATA: " if args.mc else "") + os.path.basename(args.input)
    report(lbl0, st, args.qual, args.bx, args.dr, args.minseg, args.acceptance)

    c = ROOT.TCanvas("c", "c", 800, 600)
    c.SetLogx(True); c.SetGridy(True)
    leg = ROOT.TLegend(0.50, 0.16, 0.88, 0.40)
    leg.SetBorderSize(0); leg.SetFillStyle(0)

    title = ("Per-hemisphere L1 DT trigger efficiency" + (" (acceptance cuts)" if args.acceptance else "")
             + ";offline muon p_{T} [GeV];L1 DT efficiency")

    if args.mc:
        # data vs MC: overlay the combined efficiency
        fmc, tmc = load(args.mc)
        effm_c, effm_u, effm_l, stm = compute(tmc, args.qual, args.bx, args.dr, args.minseg, args.acceptance)
        report("MC:   " + os.path.basename(args.mc), stm, args.qual, args.bx, args.dr, args.minseg, args.acceptance)

        style(eff_c,  ROOT.kBlack, 20)
        style(effm_c, ROOT.kRed,   24)
        eff_c.SetTitle(title)
        frame_and_draw(eff_c)
        effm_c.Draw("P SAME")
        leg.AddEntry(eff_c,  "Data (both hemispheres)", "lp")
        leg.AddEntry(effm_c, "MC (both hemispheres)",   "lp")
        leg.Draw()
        keep = [effm_c, effm_u, effm_l]
    else:
        style(eff_c, ROOT.kBlack, 20)
        style(eff_u, ROOT.kRed,   24)
        style(eff_l, ROOT.kBlue,  25)
        eff_c.SetTitle(title)
        frame_and_draw(eff_c)
        eff_u.Draw("P SAME"); eff_l.Draw("P SAME")
        leg.AddEntry(eff_c, "Both hemispheres (combined)", "lp")
        leg.AddEntry(eff_u, "Probe = upper", "lp")
        leg.AddEntry(eff_l, "Probe = lower", "lp")
        leg.Draw()
        keep = []

    base = os.path.join(args.outdir, args.name + "_vs_pt")
    c.SaveAs(base + ".png"); c.SaveAs(base + ".pdf")

    fout = ROOT.TFile(os.path.join(args.outdir, args.name + ".root"), "RECREATE")
    eff_c.Write("eff_comb"); eff_u.Write("eff_probeUpper"); eff_l.Write("eff_probeLower")
    if args.mc:
        keep[0].Write("eff_comb_mc"); keep[1].Write("eff_probeUpper_mc"); keep[2].Write("eff_probeLower_mc")
    fout.Close()
    print("Wrote:", base + ".png/.pdf and", args.name + ".root")

if __name__ == "__main__":
    main()
