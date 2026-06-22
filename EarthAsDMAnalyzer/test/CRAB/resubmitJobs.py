import sys, os
from optparse import OptionParser
from threading import Thread

parser = OptionParser(usage="Usage: python3 %prog codeVersion")
(opt, args) = parser.parse_args()

datasetList = []

codeVersion = sys.argv[1]

for fname in os.listdir("./"):
  if codeVersion in fname:
    datasetList.append("./" + fname)

def resubmit_job(dataset):
  outTask = f"crab resubmit -d {dataset} --maxmemory 4500 --maxjobruntime 3000"
  os.system(outTask)

threads = []
for dataset in datasetList:
  t = Thread(target=resubmit_job, args=(dataset,))
  t.start()
  threads.append(t)

for t in threads:
  t.join()