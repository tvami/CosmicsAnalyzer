import os
import glob

from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.requestName = 'Ntuplizer-0to75Theta-4to3000GeV-126X_mcRun3_2022cosmics_realistic_deco_v1_v5' #Change this

config.section_('JobType')
config.JobType.allowUndistributedCMSSW = True
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '4Analysis_muon_analyzer_cfg.py' #'4HLT2AOD.py' #Change this
config.JobType.outputFiles = ['Ntuplizer-0to75Theta-4to3000GeV_v1.root']   #['HSCP_AOD.root'] #Change This
config.JobType.disableAutomaticOutputCollection = True
config.JobType.maxMemoryMB = 4000
config.JobType.numCores = 1

config.section_('Data')
config.Data.inputDBS = 'phys03' #phys03 gave an issue
config.Data.inputDataset = '/Cosmics/lbrennan-crab_RAWtoReco-0to75Theta-4to3000GeV-126X_mcRun3_2022cosmics_realistic_deco_v1_v4-6f8577e9b9f3e6c5d2a25bde23678f33/USER'  #'/HSCP_tauPrimeCharge2e_Pythia8_TuneCP2_Mass1200_v4_ZPrimeMass6000/tvami-crab_HSCP_tauPrimeCharge2e_Pythia8_TuneCP2_Mass1200_v4_ZPrimeMass6000_HLT_2018_13TeV_v1-b403a189a2d057e62e59ed092120c7f4/USER' #Change This
config.Data.outLFNDirBase = '/store/user/lbrennan/EarthAsDM/CMSSW_13_0_10/src/Cosmics'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
#config.Data.totalUnits = 1
config.Data.ignoreLocality = True
config.Data.publication = True

config.section_('Site')
config.Site.storageSite = 'T2_US_UCSD'
config.Site.whitelist = ['T2_DE_DESY','T2_CH_CERN','T2_IT_Bari','T1_IT_*','T2_US_*', 'T3_US_FNALLPC','T2_HU_Budapest','T2_FR_*', 'T2_UK_London_IC']