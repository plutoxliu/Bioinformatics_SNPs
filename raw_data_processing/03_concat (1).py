#!/usr/bin/env python
# coding: utf-8

# In[ ]:

#python
import os
import re
###remember to make dictionary os.mkdir("Concat")
allfiles_to_concat = os.listdir("/nesi/project/uoo03773/AUS/Demulti/")
alltokeep=[]
with open("/nesi/project/uoo03773/AUS/AUSpop.txt") as f:
  for line in f:
    #print(line)
    samplename=line.split()[0]
    print (samplename)
    checkfiles=[filename for filename in allfiles_to_concat if samplename+".r" in filename and filename.startswith(samplename)]
    #print(checkfiles)
    if len (checkfiles)!=2: # that was weirdly complicated because some sample name are contained in others different ways, but the vcheck above solve it uysing the rem file
      #print(line)
      raise Exception
    else:
      nonrem = [x.replace("rem.","")for x in checkfiles] # the 2 normal files at once
    print(nonrem)
    for x in nonrem:
        assert os.path.exists("/nesi/project/uoo03773/AUS/Demulti/"+x)
        os.system("zcat "+" /nesi/project/uoo03773/AUS/Demulti/"+nonrem[0]+" /nesi/project/uoo03773/AUS/Demulti/"+nonrem[1]+" /nesi/project/uoo03773/AUS/Demulti/"+checkfiles[0]+" /nesi/project/uoo03773/AUS/Demulti/"+checkfiles[0]+ " | gzip -c > /nesi/project/uoo03773/AUS/Concat/"+samplename+".fq.gz" )


# In[ ]:




