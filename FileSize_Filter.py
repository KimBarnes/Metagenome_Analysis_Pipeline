# -*- coding: utf-8 -*-
"""
Created on Fri Dec 29 14:36:50 2017

@author: Kim
"""

import os
import glob
import shutil as sh
import sys

path=sys.path[0]
for line in open("%s/Parameters.txt" %path, "r"):
    if line.startswith("Directory"):
        directory=line.split(":")[1].strip()
    if line.startswith("File Size"):
        filesize=int(line.split(":")[1].strip())*1000
        
os.mkdir("%s/ANALYSIS" %directory)

for file in glob.glob(r"%s/*.fasta" %directory):
    if os.path.getsize("%s" %file) > filesize:
        sh.copy2("%s" %file, "%s/ANALYSIS" %directory)