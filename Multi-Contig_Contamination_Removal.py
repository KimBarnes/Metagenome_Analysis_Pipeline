# -*- coding: utf-8 -*-
"""
Created on Thu Dec 28 22:48:05 2017

@author: Kim
"""

#This script uses the longest contig within a given cluster for phylogenetic analysis using blastn, and the assumption is made that all contigs within the cluster should belong to the phylum assigned by the longest contig.
#Potential outliers are identified first by %GC content and then tetranucleotide frequency (faster)
#Each outlier is analysed for phylogenetic similarity using a blastn search, and either accepted or rejected based upon the results and similarity to the phylum assigned by the longest contig.



from __future__ import division
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from Bio.SeqUtils import GC
from Bio import SeqIO 
import glob
import pylab
import os
import seaborn as sns
import pandas as pd
from Bio.Blast.Applications import NcbiblastnCommandline
import string
import sys
import math

pylab.ioff()

def GC_calculator():
    with open ("%s_GC.txt" %name, "w") as GC_output:
         GC_output.write("Contig\tLength\tG+C%\r\n")
         for rec in SeqIO.parse("%s" %filename, "fasta"):
             GC_output.write("%s\t%i\t%s\r\n" %(rec.id, len(rec.seq), format(GC(rec.seq),".2f")))
    print ("%s GC Content Calculated" %name)
    
def GC_scattergraph():
    pylab.clf()
    Max_y=0
  #  Max_x=0
    for rec in SeqIO.parse("%s" %filename, "fasta"):
         x=GC(rec.seq)
         y=len(rec.seq)
         pylab.scatter(x,y, color = "midnightblue")
         if y>Max_y:
             Max_y=y
             #Max_x=x
    pylab.xlabel("GC%")
    pylab.ylabel("Contig Length(bp)")
    #pylab.axvline(Max_x-5, dashes=[5,5])
    #pylab.axvline(Max_x+5, dashes=[5,5])
    pylab.title("%s\n%i Contigs\nGC%% %0.1f to %0.1f" %(name, len(list(SeqIO.parse("%s" %filename, "fasta"))), min(sorted(GC(rec.seq) for rec in SeqIO.parse("%s" %filename, "fasta"))), max(sorted(GC(rec.seq) for rec in SeqIO.parse("%s" %filename, "fasta")))))
    pylab.savefig("%s_GC_Scatter.pdf" %name, dpi=300)
    pylab.close()
    print ("%s GC Graph Created" %name)
    
def BLAST():
    global input_file
    global output_file
    blast_results = NcbiblastnCommandline(query="%s" %input_file, db="/biol/programs/ncbiBLAST/db/nt", outfmt='"6 qseqid stitle evalue"', out="%s" %output_file, negative_gilist="%s/Environmental_samples.gi" %path, word_size=11, num_alignments=10, max_hsps=1, num_threads=70 )
    stdout, stderr = blast_results()
    if stdout != '' or stderr != '':
        print("BLAST returned unexpected output")
        print(stdout, stderr)
    
def Phylum_assignment():
    global BlastGenus
    global ContigPhylum
    global input_file
    ContigPhylum="No/Missing Phylum"
    if BlastGenus != "-Candidatus-":
        for file in glob.glob(r'%s' %Phylum_files):
            with open(file) as f:
                contents = f.read()
                if BlastGenus in contents:
                    Phylum=file.split("/")[-1]
                    ContigPhylum=Phylum.split(".")[0]
    elif BlastGenus == "-Candidatus-":
        for line in open("%s" %input_file):
            line=line.strip()
            GenusList=line.split()[1:3]
            Genus="-"
            for i in GenusList:
                Genus+=str(i)+" "
            BlastGenus=str(Genus[:-1])+"-"
        for file in glob.glob(r'%s' %Phylum_files):
            with open(file) as f:
                contents = f.read()
                if BlastGenus in contents:
                    Phylum=file.split("/")[-1]
                    ContigPhylum=Phylum.split(".")[0]
        if ContigPhylum=="No/Missing Phylum":
            for line in open("%s" %input_file):
                line=line.strip()
                BlastGenus="-"+line.split()[2]+"-"
                for file in glob.glob(r'%s' %Phylum_files):
                    with open(file) as f:
                        contents = f.read()
                        if BlastGenus in contents:
                            Phylum=file.split("/")[-1]
                            ContigPhylum=Phylum.split(".")[0]
    if ContigPhylum == "No/Missing Phylum":
        print("Genus not found:%s" %BlastGenus)
    return ContigPhylum

def Tetranucleotide_frequency_calculator():
    output_file = open('%s_TetranucleotideFrequency.txt' %name,'w')
    output_file.write(">Contig")
    tetralist = []
    bases = ['A', 'C', 'G', 'T']
    for b1 in bases:
        for b2 in bases:
            for b3 in bases:
                for b4 in bases:
                    tetramer = b1+b2+b3+b4
                    tetralist.append(tetramer)
    for i in tetralist:
        output_file.write("\t%s" %i)
    for rec in SeqIO.parse("%s" %filename, "fasta"):    
        gene_length = len(rec.seq)
        output_file.write("\r\n%s" %rec.id)
        for tetramer in tetralist:
            tetracount = rec.seq.count(tetramer)/gene_length*100
            output_file.write("\t%f" %tetracount)
    output_file.close()
    print("%s Tetranucleotide Frequency Calculated" %name)
    
def Tetranucleotide_frequency_heatmap():
    global number_of_contigs
    df = pd.read_table('%s_TetranucleotideFrequency.txt' %name, index_col=0 )
    cg = sns.clustermap(df, method="average", col_cluster=False)
    sns.set(font_scale=0.2)
    plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0, fontsize="small")
    plt.setp(cg.ax_heatmap.yaxis.set_label_text("Contigs"), fontsize="xx-large")
    plt.setp(cg.ax_heatmap.xaxis.get_majorticklabels(), rotation=90, fontsize="small")
    plt.setp(cg.ax_heatmap.xaxis.set_label_text("Tetranucleotides"), fontsize="xx-large")
    plt.setp(cg.ax_heatmap.set_title("%s Tetranucleotide Frequencies" %name), fontsize="xx-large")
    sns.plt.savefig("%s_TetranucleotideHeatmap.pdf" %name)
    plt.clf()
    clustered_list= cg.dendrogram_row.reordered_ind
    plt.close()
    print("%s Tetranucleotide Frequency Heatmap Created" %name)
    return clustered_list
           
#Only need to change the name of the directory
path=sys.path[0]
for line in open("%s/Parameters.txt" %path, "r"):
    if line.startswith("Directory"):
        directory=line.split(":")[1].strip()
    if line.startswith("Completeness_threshold"):
        completeness_threshold=int((line.split(":")[1]).strip())
    if line.startswith("Contamination_threshold"):
        contamination_threshold=int((line.split(":")[1]).strip())
#Phylum files created from the NCBI Taxonomy Browser Nov 2017
Phylum_files="%s/Phylum/*.txt" %path

#IDENTIFICATION OF CLUSTERS + GENERATING DIRECTORIES FOR ANALYSIS
FileList=[]
os.chdir(r'%s/ANALYSIS' %directory)
for line in open("Metagenome_CheckM_Results.txt"):
    if not line.startswith("OTU"):
        line=line.strip()
        genome_bin=line.split("\t")[0]
        completeness=float(line.split("\t")[2])
        contamination=float(line.split("\t")[3])
        if contamination > contamination_threshold and completeness > completeness_threshold:
            FileList.append(genome_bin)            

#Creating directories for each fasta file to put analysis in
os.chdir(r'%s/ANALYSIS' %directory)
for i in FileList:
    os.chdir(r'%s/ANALYSIS' %directory)
    if not os.path.exists('%s/ANALYSIS/Sequence_Composition' %directory):
        os.mkdir('%s/ANALYSIS/Sequence_Composition' %directory)
    clustername=i.split(".")[0]
    cluster_filename=i
    os.mkdir("%s/ANALYSIS/Sequence_Composition/%s" %(directory,clustername))
    for rec in SeqIO.parse("%s" %cluster_filename, "fasta"):
        os.chdir(r"%s/ANALYSIS/Sequence_Composition/%s" %(directory,clustername))
        with open("%s" %cluster_filename, "a") as file:
            file.write(">%s\r\n%s\r\n" %(rec.id,str(rec.seq)))
    os.chdir(r"%s/ANALYSIS/Sequence_Composition/%s" %(directory,clustername))

        #Assigning variables
    reject_file=True
    letter_list=sorted(string.ascii_uppercase)
    name=clustername
    filename=cluster_filename
    
    #Generating sequence composition statistics for original cluster
    GC_calculator()
    GC_scattergraph()
    Tetranucleotide_frequency_calculator()
    Tetranucleotide_frequency_heatmap()
                    
#LONGEST CONTIG + %GC CONTENT ANALYSIS
    while reject_file==True:
      Rejected_contigs=[]
      number_of_contigs=0
      ReasonsToRepeat=True
      for rec in SeqIO.parse("%s" %filename, "fasta"):
        number_of_contigs+=1
      if number_of_contigs>1:

       #Calculating the largest contig                
       Evalue=True
       LengthList=[]
       contiglist=[]
       for rec in SeqIO.parse("%s" %filename, "fasta"):
            LengthList.append((len(rec.seq),rec.id))
       while Evalue==True:
            for i in range(1,number_of_contigs+1):
             if Evalue==True:
                x=sorted(LengthList)[-i]
                BiggestContig=x[1]
            
                with open("BiggestContig.txt", "w") as output:
                    for rec in SeqIO.parse("%s" %filename, "fasta"):
                        if rec.id == BiggestContig:
                            BiggestContigLength = len(rec.seq)
                            BiggestContigGC= GC(rec.seq)
                            output.write(">%s\r\n%s" %(rec.id, str(rec.seq)))
                            Rec_id=rec.id
        
            #BLAST longest contig
                blast_results = NcbiblastnCommandline(query="BiggestContig.txt", db="/biol/programs/ncbiBLAST/db/nt", outfmt='"6 qseqid stitle evalue"', out="%s_LargestContigBlastN.txt" %clustername, negative_gilist="%s/Environmental_samples.gi" %path, word_size=11, num_alignments=1, max_hsps=1, num_threads=60 )
                stdout, stderr = blast_results()
                if stdout != '' or stderr != '':
                    print("BLAST returned unexpected output")
                    print(stdout, stderr)
#num_alignments=how many species to return
#max_hsps= Maximum number of HSPs (alignments) to keep for any single query-subject pair

            #Keep all the Blast Results in the same file
                with open("%s_BlastResults.txt" %clustername, "a") as BlastFile:
                    with open("%s_LargestContigBlastN.txt" %clustername, "r") as LongestContigBlast:
                        for line in LongestContigBlast:
                            line=line.strip()
                            BlastFile.write("\r\n\r\n%s Blast Results\r\n%s_A Blast Results\r\nQuery\tGenus\tEvalue\tPhylum\r\n\r\nLongestContig\r\n%s::\t" %(clustername, name,line))

    #Extracting BLAST results to get Genus name and Evalue + assigning Phylum to the cluster
                with open("%s_LargestContigBlastN.txt" %clustername) as BlastResults:
                    for line in BlastResults:
                        Species = line.split()[1]
                        BlastGenus= "-"+ Species.split()[0] + "-"
                        LargestContigEvalue = line.split()[-1]
                    if LargestContigEvalue=="0.0":
                        input_file="%s_LargestContigBlastN.txt" %clustername
                        LongestContigPhylum=Phylum_assignment()
                        Evalue=False
                        LongestContig=True
                    if LargestContigEvalue!="0.0":
                        if Rec_id not in contiglist:
                            contiglist.append(Rec_id)
                        if len(contiglist) == number_of_contigs:
                            reject_file=False
                            LongestContig=False
                            Evalue=False
            os.remove("BiggestContig.txt")
            os.remove("%s_LargestContigBlastN.txt" %clustername)
       if LongestContig==True:                            
        genome_length=0
        for rec in SeqIO.parse("%s" %filename, "fasta"):
            genome_length+=len(rec.seq)
 
        #Printing overview + BLAST results to file
        with open("%s_BlastResults.txt" %clustername, "a") as BlastFile:
            BlastFile.write("%s\r\n\r\n" %LongestContigPhylum)
        with open("%s_Overview.txt" %clustername, "a") as overview:
            overview.write("\r\n\r\n%s Overview\r\n\r\nNumber of Contigs:\r\n%s\r\nTotal Genome Length\r\n%ibp\r\nGC Percentage Range:\r\n%0.1f to %0.1f\r\n\r\nLongest Contig:\r\nContig Name\tLength\tGC Percentage\r\n%s\t%i\t%s\r\nGenus:\r\n%s\r\nEvalue:\r\n%s\r\nPhylum:\r\n%s\r\n\r\n" %(name,len(list(SeqIO.parse("%s" %filename, "fasta"))), genome_length, min(sorted(GC(rec.seq) for rec in SeqIO.parse("%s" %filename, "fasta"))), max(sorted(GC(rec.seq) for rec in SeqIO.parse("%s" %filename, "fasta"))),BiggestContig,BiggestContigLength,BiggestContigGC,Species,LargestContigEvalue,LongestContigPhylum))

        print ("%s Analysing GC Content Outliers" %name)
 
        #Upper and lower bounds of %GC outliers calculated
        GC_range=((max(sorted(GC(rec.seq) for rec in SeqIO.parse("%s" %filename, "fasta")))) - (min(sorted(GC(rec.seq) for rec in SeqIO.parse("%s" %filename, "fasta")))))
        GC_point=0.04*math.exp(0.14*GC_range)
        Lower_bound=float(BiggestContigGC-GC_point)
        Upper_bound=float(BiggestContigGC+GC_point)
        with open("%s_Overview.txt" %clustername, "a") as overview:
            overview.write("%s Analysis\r\nGC Percentage Analysis\r\nUpper bound\t%s\r\nLower bound\t%s\r\n\r\nGC Percentage Outliers\r\n" %(name,Upper_bound, Lower_bound))

        #Identifying outlying contigs based on %GC content
        GC_outliers=[]
        for rec in SeqIO.parse("%s" %filename, "fasta"):
            if GC(rec.seq) > Upper_bound:
                GC_outliers.append(rec.id)
            if GC(rec.seq) < Lower_bound:
                GC_outliers.append(rec.id)
                
        #Analysing each outlying contig
        if len(GC_outliers) > 0:
            with open("%s_BlastResults.txt" %clustername, "a") as BlastFile:
                BlastFile.write("%GC Outliers\r\n") 
            for rec in SeqIO.parse("%s" %filename, "fasta"):
                List=[]
                if rec.id in GC_outliers:
                        
                    #Check to see if it's already been BLASTED, should save some time                             
                    if rec.id+"\t" in open("%s_BlastResults.txt" %clustername, "r+").read():
                        BlastList=[]
                        for line in open("%s_BlastResults.txt" %clustername, "r"):
                            if rec.id+"\t" in line:
                                Line=line.strip()
                                if len(BlastList)<10:
                                    BlastList.append(Line)
                                    contig_phylum=Line.split("::")[1]
                                    contig_phylum=contig_phylum.strip()
                                    List.append(contig_phylum)
                        with open("%s_BlastResults.txt" %clustername, "a") as BlastResults:
                            for i in BlastList:
                                    BlastResults.write("\r\n%s" %i)
                    else:
                        with open("%s_TempGCOutliers.txt" %clustername, "w") as output:
                            output.write(">%s\r\n%s" %(rec.id, str(rec.seq)))
                        #BLAST GC outliers. Takes top 10 hits since contigs could be quite small
                        input_file="%s_TempGCOutliers.txt" %clustername
                        output_file="%s_TempGCOutliersBlastN.txt" %clustername
                        BLAST()

                    #Extracting BLAST results to get Genus name + assigning Phylum to the contig. Contig either accepted or rejected based on this Phylum
                        with open("%s_TempGCOutliersBlastN.txt" %clustername, "r") as BlastResults:
                            for line in BlastResults:
                                line=line.strip()
                                l=line.split()[1]
                                BlastGenus="-"+l.split()[0]+"-"
                            
                                input_file="%s_TempGCOutliersBlastN.txt" %clustername
                                contig_phylum=Phylum_assignment()
                            
                                List.append(contig_phylum)
                                with open("%s_BlastResults.txt" %clustername, "a") as BlastFile:                
                                            BlastFile.write("\r\n%s::%s" %(line, contig_phylum)) 
                    if LongestContigPhylum in List:
                            Result="Accepted"
                    else:
                            Result="Rejected"
                    if Result == "Rejected":
                            #If a contig is rejected, it's added to a list of rejected contigs
                            Rejected_contigs.append(rec.id)

                    with open("%s_Overview.txt" %clustername, "a") as overview:                
                        overview.write("%s:\t%s\r\n" %(rec.id, Result))
      
        #Write contigs that have not been rejected to "*_A.fasta"
        name="%s-%s" %(clustername,letter_list[0])
        del letter_list[0]
        with open("%s.fasta" %name, "w") as AcceptedContigs:
            for rec in SeqIO.parse("%s" %filename, "fasta"):
                 if rec.id not in Rejected_contigs:
                     AcceptedContigs.write(">%s\r\n%s\r\n" %(rec.id, str(rec.seq)))
     
        if os.path.exists("%s_TempGCOutliers.txt" %clustername):
            os.remove("%s_TempGCOutliers.txt" %clustername)
        if os.path.exists("%s_TempGCOutliersBlastN.txt" %clustername):
            os.remove("%s_TempGCOutliersBlastN.txt" %clustername)
        print ("%s GC Outliers Analysed" %name)
        
        with open("%s_Overview.txt" %clustername, "a") as overview:
            overview.write("\r\n\r\nTetranucleotide Frequency Outliers")
            
#TETRANUCLEOTIDE FREQUENCY ANALYSIS
        while ReasonsToRepeat==True:
            filename="%s.fasta" %name
            number_of_contigs=0
            for rec in SeqIO.parse("%s" %filename, "fasta"):
                number_of_contigs+=1
            if number_of_contigs>1:
              ten_contigs=[]
        
              #Revised version of "*_A.fasta"
              with open("%s" %filename, "w") as AcceptedContigs:
                    for rec in SeqIO.parse("%s" %cluster_filename, "fasta"):
                        if rec.id not in Rejected_contigs:
                            AcceptedContigs.write(">%s\r\n%s\r\n" %(rec.id, str(rec.seq)))
        
                #Recalculating tetranucleotide frequencies and heatmap
              Tetranucleotide_frequency_calculator()
              number_of_contigs=0
              for rec in SeqIO.parse("%s" %filename, "fasta"):
                number_of_contigs+=1
              if number_of_contigs>1:
                HeatmapList=Tetranucleotide_frequency_heatmap()
                if HeatmapList==0:
                    number_of_contigs=0
                    ReasonsToRepeat=False
                with open('%s_TetranucleotideFrequency.txt' %name) as file:
                    ContigList=[]
                    for line in file:
                        if not line.startswith(">Contig"):
                            line=line.split()[0]
                            ContigList.append(line)
                ContigList=list(enumerate(ContigList))
            
                #Extracting ID's of 10 most outlying contigs from the tetranucleotide frequency heatmap
                ReorderedContigsList=[]
                for i in HeatmapList:
                    for item in ContigList:
                        if item[0] == i:
                            ReorderedContigs= item[1]
                            ReorderedContigsList.append(ReorderedContigs)
        
                TNF_Outliers=[]
                for i in ReorderedContigsList[0:5]:
                    TNF_Outliers.append(i)
                for i in ReorderedContigsList[-6:-1]:
                    TNF_Outliers.append(i)
                       
                #Extract sequences of TNF_Outliers + BLAST
                for rec in SeqIO.parse("%s" %filename, "fasta"):
                    List=[]
                    if rec.id in TNF_Outliers:
                        if rec.id+"\t" in open("%s_BlastResults.txt" %clustername, "r+").read():
                            BlastList=[]
                            for line in open("%s_BlastResults.txt" %clustername, "r"):
                                if rec.id+"\t" in line:
                                    Line=line.strip()
                                    if len(BlastList)<10:
                                        BlastList.append(Line)
                                        contig_phylum=Line.split("::")[1]
                                        contig_phylum=contig_phylum.strip()
                                        List.append(contig_phylum)
                            with open("%s_BlastResults.txt" %clustername, "a") as BlastResults:
                                    for i in BlastList:
                                        BlastResults.write("\r\n%s" %i)
                        
                        else:
                            with open("%s_TempTNFOutliers.txt" %clustername, "w") as TNFContig_output:
                                TNFContig_output.write(">%s\r\n%s" %(rec.id, str(rec.seq)))  
                            
                            input_file="%s_TempTNFOutliers.txt" %clustername
                            output_file="%s_TempTNFOutliersBlastN.txt" %clustername
                            BLAST()
        
                        #Extracting BLAST results to get Genus name + assigning Phylum to the contig. Contig either accepted or rejected based on this Phylum
                            with open("%s_TempTNFOutliersBlastN.txt" %clustername) as BlastResults:
                                for line in BlastResults:
                                    line=line.strip()
                                    l=line.split()[1]
                                    BlastGenus="-"+l.split()[0]+"-"
                                    input_file="%s_TempTNFOutliersBlastN.txt" %clustername
                                    contig_phylum=Phylum_assignment()
                                    List.append(contig_phylum)
                                    with open("%s_BlastResults.txt" %clustername, "a") as BlastFile:                
                                        BlastFile.write("\r\n%s::\t%s" %(line, contig_phylum))    
                                        
                        if LongestContigPhylum in List:
                            Result="Accepted"
                        else:
                            Result="Rejected"
                        #If it finds a "Rejected" contig in the outlying 10 contigs: removes rejected sequence and reclusters in a looping fashion untill all 10 are accepted
                        if Result == "Rejected":
                            ten_contigs.append("R")
                            Rejected_contigs.append(rec.id)
                        if not Result == "Rejected":
                            ten_contigs.append("A")
                            
                        with open("%s_Overview.txt" %clustername, "a") as overview:                
                            overview.write("\r\n%s:\t%s" %(rec.id, Result))
                                                        
                if "R" not in ten_contigs:
                    ReasonsToRepeat=False
                else:
                    with open("%s_Overview.txt" %clustername, "a") as overview:
                        overview.write("\r\n")
            else:
                ReasonsToRepeat=False
        #Once all 10 contigs are accepted, the loop is broken and the final version of "*_A.fasta" written
        with open("%s" %filename, "w") as AcceptedContigs:
            for rec in SeqIO.parse("%s" %cluster_filename, "fasta"):
                if rec.id not in Rejected_contigs:
                    AcceptedContigs.write(">%s\r\n%s\r\n" %(rec.id, str(rec.seq)))
                    
        genome_length=0
        for rec in SeqIO.parse("%s" %filename, "fasta"):
            genome_length+=len(rec.seq)
            
        with open("%s_Overview.txt" %clustername, "a") as overview:
            overview.write("\r\n\r\n%s Overview\r\n\r\nNumber of Contigs:\r\n%s\r\nTotal Genome Length\r\n%ibp\r\nGC Percentage Range:\r\n%0.1f to %0.1f\r\n\r\nLongest Contig:\r\nContig Name\tLength\tGC Percentage\r\n%s\t%i\t%s\r\nGenus:\r\n%s\r\nEvalue:\r\n%s\r\nPhylum:\r\n%s\r\n\r\n" %(name,len(list(SeqIO.parse("%s" %filename, "fasta"))), genome_length, min(sorted(GC(rec.seq) for rec in SeqIO.parse("%s" %filename, "fasta"))), max(sorted(GC(rec.seq) for rec in SeqIO.parse("%s" %filename, "fasta"))),BiggestContig,BiggestContigLength,BiggestContigGC,Species,LargestContigEvalue,LongestContigPhylum))
            
        if os.path.exists("%s_TempTNFOutliers.txt" %clustername):                                
            os.remove("%s_TempTNFOutliers.txt" %clustername)
        if os.path.exists("%s_TempTNFOutliersBlastN.txt" %clustername):
            os.remove("%s_TempTNFOutliersBlastN.txt" %clustername)
        print("%s Tetranucleotide Frequency Analysis Completed" %name)
            
        #Calculating %GC content, contig length + %GC scatter graph for final version of "*_A.fasta". Purely visual 
        GC_calculator()
        GC_scattergraph()

        #Writing all rejected contigs from [Rejected_contigs] to a file
        with open("%s_RejectedContigs.fasta" %clustername, "w") as RejectedContigs:
                for rec in SeqIO.parse("%s.fasta" %clustername, "fasta"):
                    if rec.id in Rejected_contigs:
                        RejectedContigs.write(">%s\r\n%s\r\n" %(rec.id, str(rec.seq)))

        #If "*_RejectedContigs.fasta" contains anything, a loop analysing the contigs within it is started. Allows the script to analyse clusters with more than two genomes/phyla
        cluster_filename="%s_RejectedContigs.fasta" %clustername
        name="%s_RejectedContigs" %clustername
        filename=cluster_filename
        
        number_of_contigs=0
        for rec in SeqIO.parse("%s_RejectedContigs.fasta" %clustername, "fasta"):
            number_of_contigs+=1
        if number_of_contigs<1:
                reject_file=False
                os.remove("%s_RejectedContigs.fasta" %clustername)
        else:
                reject_file=True
      else:
        with open("%s_UnassignedContigs.fasta" %clustername, "a") as output:
                  for rec in SeqIO.parse("%s" %filename, "fasta"):
                      output.write(">%s\r\n%s\r\n" %(rec.id, str(rec.seq)))
        reject_file=False
        if os.path.exists("%s_RejectedContigs.fasta" %clustername):
            os.remove("%s_RejectedContigs.fasta" %clustername)