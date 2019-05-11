# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 11:33:58 2016

@author: sophiali
"""


# Read files, add columns, join them
import os, glob, re

# Defined the path and change directory
currentPath = '/Volumes/Sophia_ResearchData/201800721_MitoRiboSeq_Serdpr_Chlorm_primerchange/data/codon_count_table/'
os.chdir(currentPath)

# Find all the files
ListofFiles = glob.glob('*codon_count_table.txt');
print(ListofFiles)
with open('AllCodonCount.txt',"w") as fout:
    # Loop through all the tabular files found in that folder
    i = 0
    for f in ListofFiles:
        with open(f,"r") as fin: 
            # Extract the condition and sequence type from the file name
            condition = re.split('\_codon_count_table.txt',f)[0]
            
            # For each line, add one extra value which is the condition
            first_line = fin.readline()
            if i ==0:
                first_line = first_line.rstrip('\n') + '\t' + 'sample' + '\n'
                fout.write(first_line)
            for line in fin:
                
                l = line.strip('\n').split('\t')
                new_line = '\t'.join(l) + '\t' + condition +  '\n' 
                                
#                new_line = line[:-3] + '\t' + condition + '\t' + seq_type + '\t.\n'
                fout.write(new_line)
        i = i + 1
