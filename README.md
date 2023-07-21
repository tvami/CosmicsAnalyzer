# CosmicsAnalyzer

## Instructions to run the code on lxplus

```
cmsrel CMSSW_13_0_10
cd CMSSW_13_0_10/src/
cmsenv
mkdir CosmicsAnalyzer
git clone -b main git@github.com:tvami/CosmicsAnalyzer.git CosmicsAnalyzer/MyAnalyzer
scram b -j
cd CosmicsAnalyzer/MyAnalyzer/test/
// just an example file
xrdcp root://cms-xrd-global.cern.ch//store/data/Run2023D/Cosmics/RAW-RECO/CosmicSP-PromptReco-v1/000/369/811/00000/2ad63d9f-234b-4c68-8c18-49d485d42bc7.root .
cmsRun muon_analyzer_cfg.py
```
