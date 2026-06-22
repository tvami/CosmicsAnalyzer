import sys, os, time, re
import numpy as np
from threading import Thread

datasetList = [
    "/Cosmics/Run2022B-CosmicTP-PromptReco-v1/RAW-RECO",
    "/Cosmics/Run2022C-CosmicTP-PromptReco-v1/RAW-RECO",
    "/Cosmics/Run2022D-CosmicTP-PromptReco-v1/RAW-RECO",
    "/Cosmics/Run2022D-CosmicTP-PromptReco-v2/RAW-RECO",
    "/Cosmics/Run2022D-CosmicTP-PromptReco-v3/RAW-RECO",
    "/Cosmics/Run2022E-CosmicTP-PromptReco-v1/RAW-RECO",
    "/Cosmics/Run2022F-CosmicTP-PromptReco-v1/RAW-RECO",
    "/Cosmics/Run2022G-CosmicTP-PromptReco-v1/RAW-RECO",
    "/Cosmics/Run2023A-CosmicTP-PromptReco-v1/RAW-RECO",
    "/Cosmics/Run2023A-CosmicTP-PromptReco-v2/RAW-RECO",
    "/Cosmics/Run2023B-CosmicTP-PromptReco-v1/RAW-RECO",
    "/Cosmics/Run2023C-CosmicTP-PromptReco-v1/RAW-RECO",
    "/Cosmics/Run2023C-CosmicTP-PromptReco-v2/RAW-RECO",
    "/Cosmics/Run2023C-CosmicTP-PromptReco-v3/RAW-RECO",
    "/Cosmics/Run2023C-CosmicTP-PromptReco-v4/RAW-RECO",
    "/Cosmics/Run2023D-CosmicTP-PromptReco-v1/RAW-RECO",
    "/Cosmics/Run2023D-CosmicTP-PromptReco-v2/RAW-RECO",
    "/Cosmics/Run2023E-CosmicTP-PromptReco-v1/RAW-RECO",
    "/Cosmics/Run2023F-CosmicTP-PromptReco-v1/RAW-RECO",
    "/Cosmics/Run2024A-CosmicTP-PromptReco-v1/RAW-RECO",
    "/Cosmics/Run2024B-CosmicTP-PromptReco-v1/RAW-RECO",
    "/Cosmics/Run2024C-CosmicTP-PromptReco-v1/RAW-RECO",
    "/Cosmics/Run2024D-CosmicTP-PromptReco-v1/RAW-RECO",
    "/Cosmics/Run2024E-CosmicTP-PromptReco-v1/RAW-RECO",
    "/Cosmics/Run2024E-CosmicTP-PromptReco-v2/RAW-RECO",
    "/Cosmics/Run2024F-CosmicTP-PromptReco-v1/RAW-RECO",
    "/Cosmics/Run2024G-CosmicTP-PromptReco-v1/RAW-RECO",
    "/Cosmics/Run2024H-CosmicTP-PromptReco-v1/RAW-RECO",
    "/Cosmics/Run2024I-CosmicTP-PromptReco-v1/RAW-RECO",
    "/Cosmics/Run2024I-CosmicTP-PromptReco-v2/RAW-RECO",
    "/Cosmics/Run2024J-CosmicTP-PromptReco-v1/RAW-RECO",
    "/Cosmics/Run2025A-CosmicTP-PromptReco-v1/RAW-RECO",
    "/Cosmics/Run2025A-CosmicTP-PromptReco-v2/RAW-RECO",
    "/Cosmics/Run2025B-CosmicTP-PromptReco-v1/RAW-RECO",
    "/Cosmics/Run2025C-CosmicTP-PromptReco-v1/RAW-RECO",
    "/Cosmics/Run2025C-CosmicTP-PromptReco-v2/RAW-RECO",
    "/Cosmics/Run2025D-CosmicTP-PromptReco-v1/RAW-RECO",
    "/Cosmics/Run2025E-CosmicTP-PromptReco-v1/RAW-RECO",
    "/Cosmics/Run2025F-CosmicTP-PromptReco-v1/RAW-RECO",
    "/Cosmics/Run2025F-CosmicTP-PromptReco-v2/RAW-RECO",
    "/Cosmics/Run2025G-CosmicTP-PromptReco-v1/RAW-RECO",
    "/Cosmics/Run2026A-CosmicTP-PromptReco-v1/RAW-RECO",
    "/Cosmics/Run2026B-CosmicTP-PromptReco-v1/RAW-RECO",
    "/Cosmics/Run2026C-CosmicTP-PromptReco-v1/RAW-RECO",
    "/Cosmics/Run2026D-CosmicTP-PromptReco-v1/RAW-RECO"
]

if not os.path.exists("submittedConfigs"): os.makedirs("submittedConfigs")

if not os.path.exists("4crab_Template.py"):
  TEMPLATE = '''
import os
import glob

from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.workArea = 'crab_projects'
config.General.requestName = 'Ntuplizer-ROVIDMINTA_v5' #Change this

config.section_('JobType')
config.JobType.allowUndistributedCMSSW = True
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '4N-Data_cfg.py'
config.JobType.outputFiles = ['Ntuplizer-Cosmics-Data.root']
config.JobType.disableAutomaticOutputCollection = True
config.JobType.numCores = 1

config.section_('Data')
config.Data.inputDBS = 'global'
config.Data.inputDataset = 'MINTA'
config.Data.outLFNDirBase = '/ceph/cms/store/user/tvami/EarthAsDM/Ntuples/Ntuples_v5.0.0/' #Change this
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.ignoreLocality = True
config.Data.publication = False

config.section_('Site')
config.Site.storageSite = 'T2_US_UCSD'
config.Site.whitelist = ['T2_DE_DESY','T2_CH_CERN','T2_IT_Bari','T1_IT_*','T2_US_*', 'T3_US_FNALLPC','T2_HU_Budapest','T2_FR_*', 'T2_UK_London_IC']
  '''

  with open("4crab_Template.py", "w") as text_file:
      text_file.write(TEMPLATE)

def task(i):
  # print("Submit for sample "+i)
  os.system("cp 4crab_Template.py 4crab_toSubmit"+str(i.replace("/","_"))+".py")
  shortSampleName = i[1:(i.find('RAW-RECO'))-1].replace("/","_")
  replaceROVIDMINTA = "sed -i 's/ROVIDMINTA/"+shortSampleName+"/g' 4crab_toSubmit"+str(i.replace("/","_"))+".py"
  os.system(replaceROVIDMINTA)
  replaceMINTA = "sed -i 's/MINTA/"+i.replace("/","\/")+"/g' 4crab_toSubmit"+str(i.replace("/","_"))+".py"
  os.system(replaceMINTA)
  os.system("crab submit -c 4crab_toSubmit"+str(i.replace("/","_"))+".py")
  os.system("mv 4crab_toSubmit"+str(i.replace("/","_"))+".py submittedConfigs/.")

threads = []
for dataset in datasetList:
  t = Thread(target=task, args=(dataset,))
  threads.append(t)
  t.start()

for t in threads:
    t.join()

os.system("rm 4crab_Template.py")
