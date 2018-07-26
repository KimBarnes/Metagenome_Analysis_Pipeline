# -*- coding: utf-8 -*-#
"""
Created on Thu Nov 23 15:04:51 2017

@author: Kim
"""
from __future__ import division
import matplotlib
matplotlib.use("Agg")
import glob
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import os
import sys

#Only need to change the name of the directory
path=sys.path[0]
for line in open("%s/Parameters.txt" %path, "r"):
    if line.startswith("Directory"):
        directory=line.split(":")[1].strip()

os.mkdir("%s/ANALYSIS/Metabolic_Pathways" %directory)

#Provides overview of the %completeness of pathways in each cluster
with open("%s/ANALYSIS/Metabolic_Pathways/Metagenome_Metabolic_Pathway_Analysis.tsv" %directory, "a") as output:
    output.write("Cluster Name")


for file in glob.glob(r'%s/Metabolic_Pathways/*/*.txt' %path):
    metabolic_filename=file.split("/")[-1]
    metabolic_pathway=metabolic_filename.split(".")[0]
    with open ("%s/ANALYSIS/Metabolic_Pathways/Metagenome_Metabolic_Pathway_Analysis.tsv" %directory, "a") as output:
        output.write("\t%s" %metabolic_pathway)

for file in glob.glob(r'%s/ANALYSIS/Genome_Annotation/*/Prokka/*/*.gbf'%directory):
     sequence_filename=file.split("/")[-1]
     sequence_name=sequence_filename.split(".")[0] 
     with open("%s/ANALYSIS/Metabolic_Pathways/%s_pathway_genes.txt" %(directory,sequence_name),"a") as p_output:
         p_output.write("%s" %sequence_name)
         
     with open("%s/ANALYSIS/Metabolic_Pathways/Metagenome_Metabolic_Pathway_Analysis.tsv" %directory, "a") as output:
         output.write("\n%s" %sequence_name)
     for metabolic_file in glob.glob(r'%s/Metabolic_Pathways/*/*.txt' %path):
        score=0
        text_list=[]
        metabolic_filename=metabolic_file.split("/")[-1]
        metabolic_pathway=metabolic_filename.split(".")[0]
        metabolic_pathway_length=sum(1 for line in open('%s' %metabolic_file) if not line.startswith("-"))
        for metabolic_line in open("%s" %metabolic_file):
            if not metabolic_line.startswith("-"):
                Result=False
                gene_list=[]
                metabolic_line=metabolic_line.strip()
                step_name=metabolic_line.split("!")
                for i in step_name:
                    redundant_enzymes=i.split("|")[1]
                    multiple_genes=redundant_enzymes.split(",")
                    for genes in multiple_genes:
                            gene='"'+genes+'"'
                            gene_='"'+genes+'_'
                            gene_list.append(gene)
                            gene_list.append(gene_)
                for sequence_line in open("%s" %file):
                    for gene in gene_list:
                        if gene in sequence_line:
                            if score>=metabolic_pathway_length:
                                break
                            Result=True
                            score+=1
                            text_list.append(sequence_line)
                    if Result==True:
                            break
        pathway_completeness=score/metabolic_pathway_length
        with open("%s/ANALYSIS/Metabolic_Pathways/Metagenome_Metabolic_Pathway_Analysis.tsv" %directory, "a") as m_output:
            m_output.write("\t%.2f" %pathway_completeness)
        if pathway_completeness > 0:
            with open("%s/ANALYSIS/Metabolic_Pathways/%s_pathway_genes.txt" %(directory,sequence_name),"a") as p_output:
                p_output.write("\r\n\r\n\r\n%s\r\n\r\n" %metabolic_pathway)
            for metabolic_line in open("%s" %metabolic_file):
                if metabolic_line.startswith("-"):
                    with open("%s/ANALYSIS/Metabolic_Pathways/%s_pathway_genes.txt" %(directory,sequence_name),"a") as p_output:
                        p_output.write("\r\n%s" %metabolic_line)                
                if not metabolic_line.startswith("-"):
                    gene_list=[]
                    metabolic_line=metabolic_line.strip()
                    step_name=metabolic_line.split("!")
                    for i in step_name:
                        redundant_enzymes=i.split("|")[1]
                        enzyme_name=i.split("|")[0]
                        multiple_genes=redundant_enzymes.split(",")
                        for genes in multiple_genes:
                                gene='"'+genes+'"'
                                gene_='"'+genes+'_'
                                gene_list.append((gene,enzyme_name))
                                gene_list.append((gene_,enzyme_name))
                    for sequence_line in open("%s" %file):
                        for i in gene_list:
                            if str(i[0]) in sequence_line:
                                gene=sequence_line.split('"')[1]
                                with open("%s/ANALYSIS/Metabolic_Pathways/%s_pathway_genes.txt" %(directory,sequence_name),"a") as p_output:
                                    p_output.write("\t-%s-\t%s\t" %(gene,str(i[1])))
                        
os.chdir("%s/ANALYSIS/Metabolic_Pathways" %directory)
df = pd.read_table("Metagenome_Metabolic_Pathway_Analysis.tsv", index_col=0 )
sns.set(color_codes=True)
sns.set(font_scale=0.2)
cg = sns.clustermap(df, col_cluster=False, row_cluster=False, cmap="Blues")
plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0, fontsize="small")
plt.setp(cg.ax_heatmap.yaxis.set_label_text("Clusters"), fontsize="x-large")
plt.setp(cg.ax_heatmap.xaxis.get_majorticklabels(), rotation=40, fontsize="small")
plt.setp(cg.ax_heatmap.xaxis.set_label_text("Metabolic Pathways"), fontsize="x-large")
plt.setp(cg.ax_heatmap.set_title("Metagenome Metabolic Pathway Analysis"), fontsize="xx-large")
plt.savefig("Metagenome Metabolic Pathway Analysis.pdf")
#plt.show()
#plt.clf()

