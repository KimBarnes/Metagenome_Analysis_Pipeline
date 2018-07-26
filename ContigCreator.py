# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 15:26:21 2017

@author: Kim
"""

#Creates contigs from whole sequences
#if genome isn't whole it's a problem and you'll need to remove the extra > lines


from random import randint
import glob
from Bio.SeqUtils import GC
from Bio import SeqIO
import pylab
pylab.ioff()



for file in glob.glob(r"C:/Users/Kim/Documents/Masters/Python In Progress/Control Sequences/Whole Sequences/*.fasta"):
 filename=file.split("\\")[1]
 name=filename.split(".")[0]
 with open("%s" %filename, "r") as f:
    line=f.readlines()
    lineNumber= 1
    x=25
    fragmentLength=25
    fragmentNumber=10000
    for i in list(range(1,1000000)):
            fragmentNumber+=1
            if (i % 30)==0:
                fragmentLength=fragmentLength+round(0.1*fragmentLength)
            else:
                fragmentLength=fragmentLength
            contigStart=lineNumber
            contigEnd=lineNumber+fragmentLength
            fragment=line[contigStart:contigEnd]
            if len(fragment) > 24:
              with open("%s_pieces.txt" %name, "a") as output_file:
                output_file.write(">%s_fragment_%s\n" %(name,fragmentNumber))
                for item in fragment:
                    output_file.writelines("%s" %item)
            lineNumber=contigEnd
            if len(fragment) < 24:
                break