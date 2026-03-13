#!/usr/bin/env python
# coding: utf-8

import os
import re

###remember to make dictionary os.mkdir("Concat")
demulti = "/nesi/project/uoo03773/AUS/Demulti/" ##Demultiplexed files
pop = "/nesi/project/uoo03773/AUS/AUSpop.txt" ##population info for individual samples
out = "/nesi/project/uoo03773/AUS/Concat/" ##output dictionary

allfiles_to_concat = os.listdir(demulti)
alltokeep = []

with open(pop) as f:
  for line in f:
    samplename = line.split()[0] #This assumes the sample name is always the first column.
    print(samplename)
    checkfiles = [filename for filename in allfiles_to_concat if samplename+".r" in filename and filename.startswith(samplename)] #the filename contains and starts with sample name

    if len(checkfiles) != 2: #check each sample read has one normal and one rem (2 forward and 2 reverse for each sample)
      raise Exception
    else:
      nonrem = [x.replace("rem.", "") for x in checkfiles] #rename the .rem files back to normal file name (in processing memory not the actual files), in order to concatenate them into one
    print(nonrem)
    
    for x in nonrem: #concatenate forward + reverse + forward.rem + reverse.rem into one file for each sample
        assert os.path.exists(demulti + x)
        os.system("zcat " + demulti+nonrem[0] + " " + demulti+nonrem[1] + " " + demulti+checkfiles[0] + " " + demulti+checkfiles[1] + " | gzip -c > " + out + samplename + ".fq.gz")
