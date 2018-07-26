# -*- coding: utf-8 -*-
"""
Created on Sun Feb 11 14:54:07 2018

@author: Kim
"""

import sys
import os
import glob
import shutil as sh

path=sys.path[0]
for line in open("%s/Parameters.txt" %path, "r"):
    if line.startswith("Directory"):
        directory=line.split(":")[1].strip()
    if line.startswith("Completeness_threshold"):
        completeness_threshold=int((line.split(":")[1]).strip())
    if line.startswith("Contamination_threshold"):
        contamination_threshold=int((line.split(":")[1]).strip())
        
#IDENTIFICATION OF CLUSTERS + GENERATING DIRECTORIES FOR ANALYSIS
os.chdir(r'%s/ANALYSIS' %directory)
with open("Metagenome_CheckM_Results.txt", "a") as output:
    output.write("OTU Name\tGenome Size\tTaxonomy\tCompleteness (%)\tContamination (%)\r\n")
os.chdir(r'%s/ANALYSIS/CheckM/storage' %directory)
for line in open("bin_stats_ext.tsv"):
    genome_bin=line.split("\t")[0]+".fasta"
    genome_size=line.split("'Genome size': ")[1]
    genome_size=float(genome_size.split(",")[0])
    completeness=line.split("'Completeness': ")[1]
    completeness=float(completeness.split(",")[0])
    contamination=line.split("'Contamination': ")[1]
    contamination=float(contamination.split(",")[0])
    kingdom=line.split("'marker lineage': ")[1]
    kingdom=kingdom.split(",")[0]
    if kingdom != "'root'":
        kingdom=kingdom.split("__")[1]
        kingdom=kingdom.split("'")[0]
    os.chdir(r'%s/ANALYSIS' %directory)
    with open("Metagenome_CheckM_Results.txt", "a") as output:
        output.write("%s\t%s\t%s\t%s\t%s\r\n" %(genome_bin,genome_size,kingdom,completeness,contamination))

os.chdir(r'%s/ANALYSIS' %directory)
os.chdir(r'%s/ANALYSIS' %directory)
for line in open("Metagenome_CheckM_Results.txt"):
    if not line.startswith("OTU"):
        line=line.strip()
        genome_bin=line.split("\t")[0]
        kingdom=line.split("\t")[2]
        completeness=float(line.split("\t")[3])
        contamination=float(line.split("\t")[4])
        if contamination < contamination_threshold and completeness > completeness_threshold:
            if not os.path.exists('%s/ANALYSIS/Genome_Annotation' %directory):
                os.mkdir(r'%s/ANALYSIS/Genome_Annotation'%directory)
            if not os.path.exists('%s/ANALYSIS/Genome_Annotation/Bacteria' %directory):
                os.mkdir(r'%s/ANALYSIS/Genome_Annotation/Bacteria' %directory)
            if not os.path.exists('%s/ANALYSIS/Genome_Annotation/Archaea' %directory):
                os.mkdir(r'%s/ANALYSIS/Genome_Annotation/Archaea' %directory)
            for file in glob.glob(r"%s/ANALYSIS/%s" %(directory,genome_bin)):
                if kingdom != "'root'":
                    if "rcha" in kingdom:
                        sh.copy2("%s" %file, "%s/ANALYSIS/Genome_Annotation/Archaea" %directory)
                    else:
                        sh.copy2("%s" %file, "%s/ANALYSIS/Genome_Annotation/Bacteria" %directory)
                        