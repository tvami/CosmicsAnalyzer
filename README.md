# CosmicsAnalyzer

## Instructions to run the code on lxplus

```
cmsrel CMSSW_13_0_10
cd CMSSW_13_0_10/src/
cmsenv
git clone -b main git@github.com:tvami/CosmicsAnalyzer.git
scram b -j
cd CosmicsAnalyzer/EarthAsDMAnalyzer/test/
// just an example file
xrdcp root://cms-xrd-global.cern.ch//store/data/Run2023D/Cosmics/RAW-RECO/CosmicSP-PromptReco-v1/000/369/811/00000/2ad63d9f-234b-4c68-8c18-49d485d42bc7.root .
cmsRun muon_analyzer_cfg.py
```

## Ntuplizer versions

The ntuplizer is tagged in git (`git tag`). Each version adds branches on top of the previous one; older ntuples remain readable, they simply lack the newer branches.

| Version | Date | Main additions |
| --- | --- | --- |
| v2 | 2025-11-11 | More general-track info, and the `HLT_Random` trigger (PR #17) |
| v3 | 2025-11-12 | Extra info on the muon-to-general-track match (`muon_fromGenTrack_*`), removed the R-squared variable (PR #18) |
| v4 | 2026-02-13 | Number of valid tracker hits on the muon object (`muon_numberOfValidHits`); `track_n` made `uint` (PR #19) |
| v5 | 2026-06-22 | L1 trigger info for the per-hemisphere L1 trigger-efficiency measurement on cosmics: uGMT muon candidates (`L1muon_*`) and L1 DT Local Trigger primitives (`dtTrigPh_*`, one per DT chamber, with chamber global position for hemisphere assignment). Also the magnetic field `bField` [T]: in data read from the online DCS magnet current (so B-off runs show ~0), in MC the nominal 3.8 T (set `isData=0` for MC, the DCS current is dummy there) |

## Per-hemisphere L1 DT trigger efficiency (`L1TriggerEfficiency.py`)

`EarthAsDMAnalyzer/test/L1TriggerEfficiency.py` measures the L1 DT Local Trigger
efficiency with a tag-and-probe on cosmics that cross both hemispheres, following
the CMS cosmic commissioning paper (JINST 5 (2010) T03002, Fig. 12): tag = the DT
trigger fired in one half, probe = did it also fire in the other half, as a
function of the offline muon pT. It uses the `dtTrigPh_*` branches (v5+).

### 1. Produce the ntuples (the `dtTrigPh_*` branches must be filled)

The DT trigger primitives are **not** stored in standard files, so the ntuplizer
config has to provide them. The relevant config snippets are already set up in
`test/4N-Data_cfg.py`, `test/CRAB/4N-Data_cfg.py` and `test/CRAB/4N-CentralMC_cfg.py`:

- **Data (RAW-RECO)** — re-unpack from the raw FED data:
  ```python
  process.load('EventFilter.L1TRawToDigi.bmtfDigis_cfi')
  process.muonPhiAnalyzer.dtTrigPhiCollection = cms.InputTag("bmtfDigis")
  process.p = cms.Path(process.bmtfDigis * process.muonPhiAnalyzer)
  ```
- **MC (GEN-SIM-RECO / with `muonDTDigis`)** — re-emulate from the DT digis:
  ```python
  process.load('L1Trigger.DTTrigger.dtTriggerPrimitiveDigis_cfi')
  process.dtTriggerPrimitiveDigis.digiTag = cms.InputTag("muonDTDigis")
  process.muonPhiAnalyzer.dtTrigPhiCollection = cms.InputTag("dtTriggerPrimitiveDigis")
  process.p = cms.Path(process.dtTriggerPrimitiveDigis * process.muonPhiAnalyzer)
  ```

### 2. Run the measurement

```
cd CosmicsAnalyzer/EarthAsDMAnalyzer/test/

# single sample (data or MC): combined + per-hemisphere (upper/lower) curves
python3 L1TriggerEfficiency.py -i Ntuplizer-Data-Run2025A.root -n L1DTTrigEff_Data

# data vs MC overlay of the combined efficiency
python3 L1TriggerEfficiency.py -i Ntuplizer-Data-Run2025A.root \
                               --mc Ntuplizer-MC-CosmicToMu.root \
                               -n L1DTTrigEff_DataMC

# with paper-like fiducial cuts (pT > 5, >= 2 DT stations, away from phi-cracks)
python3 L1TriggerEfficiency.py -i data.root --mc mc.root --acceptance -n L1DTTrigEff_DataMC_acc
```

Each run prints the inclusive efficiency (Clopper-Pearson interval) and writes
`<name>_vs_pt.png`, `<name>_vs_pt.pdf` and `<name>.root` (holding the `TEfficiency`
objects `eff_comb`, `eff_probeUpper`, `eff_probeLower`).

Options: `-i` input ntuple, `--mc` second ntuple to overlay, `-n` output basename,
`-o` output dir, `-t` tree name (default `muonPhiAnalyzer/tree`), `--qual` min DT
primitive quality code (default 4), `--bx` max |bx| (default 1), `--dr` max dR of
the primitive-to-leg match (default 0.4), `--minseg` min DT segments per hemisphere
(default 1), `--acceptance` apply the fiducial cuts.

> Note: on an `el8`/`el9` mismatched host, run inside `cmssw-el8`. If the input
> ntuples live on `/ceph` (not mounted in the container), bind it first, e.g.
> `export SINGULARITY_BINDPATH="/ceph:/pnfs"` and refer to the files under `/pnfs`.

### 3. BX / quality scan (timing vs intrinsic, `L1TriggerEfficiencyBXscan.py`)

`L1TriggerEfficiencyBXscan.py` tabulates the inclusive efficiency and the data/MC
scale factor (`SF = eff_data / eff_MC`) over a grid of BX windows and quality
thresholds in a single pass, and writes a plot of efficiency vs BX window:

```
python3 L1TriggerEfficiencyBXscan.py -d Ntuplizer-Data-Run2025A.root \
                                     -m Ntuplizer-MC-CosmicToMu.root \
                                     -n L1DTTrigEff_BXscan
```

This disentangles the two sources of any data/MC difference: the BX window probes
**timing**, the quality threshold probes the **intrinsic** local-trigger quality.
For the Run2025A vs CosmicToMu test samples the intrinsic efficiency agrees
(`SF ~ 1.0` at `|bx|<=2`, quality `>= 4`); the apparent gap at the nominal
`|bx|<=1` is a BX/timing-synchronisation effect (a cosmic crosses the two
hemispheres ~1-2 BX apart, which the MC emulation models too tightly).
Options: `-d` data ntuple, `-m` MC ntuple, `--bxlist`, `--quallist`, `--plotqual`,
`--dr`, `--minseg`, `-n`/`-o` output name/dir.

