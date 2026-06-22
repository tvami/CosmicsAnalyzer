import sys, os, time, re
import numpy as np
from threading import Thread

datasetList = []

datasetListSignal = [
  "/Signal/smasanam-crab_RAWtoReco-91to180Theta-3to4000GeV-126X_mcRun3_2022cosmics_realistic_deco_v4_s3-6f8577e9b9f3e6c5d2a25bde23678f33/USER"
]

datasetListBkg = [
  "/Cosmics/lbrennan-crab_RAWtoReco-0to75Theta-4to3000GeV-126X_mcRun3_2022cosmics_realistic_deco_SimHits_v1_v9-6f8577e9b9f3e6c5d2a25bde23678f33/USER"
]

datasetList.extend(datasetListBkg)

if not os.path.exists("submittedConfigs"): os.makedirs("submittedConfigs")

if not os.path.exists("4crab_Template.py"):
  TEMPLATE = '''
import os
import glob

from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.workArea = 'crab_projects'
config.General.transferLogs = True
config.General.requestName = 'Ntuplizer-ROVIDMINTA_v3' #Change this

config.section_('JobType')
config.JobType.allowUndistributedCMSSW = True
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '4N-CentralMC_cfg.py'
config.JobType.outputFiles = ['Ntuplizer-Cosmics-MC.root']
config.JobType.disableAutomaticOutputCollection = True
#config.JobType.maxMemoryMB = 4500
config.JobType.numCores = 1

config.section_('Data')
config.Data.inputDBS = 'phys03'
config.Data.inputDataset = 'MINTA'
config.Data.outLFNDirBase = '/store/user/tvami/EarthAsDM/'
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
  # Extract a shorter sample name (e.g., "RAWtoReco-91to180Theta-3to4000GeV")
  temp = i.split('crab_')[1] if 'crab_' in i else i
  shortSampleName = temp.split('-126X')[0].replace('/', '_')
  # Ensure it's under 100 chars when combined with prefix
  if len("Ntuplizer-" + shortSampleName + "_v3") > 100:
    shortSampleName = shortSampleName[:85]  # Leave room for prefix and suffix
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