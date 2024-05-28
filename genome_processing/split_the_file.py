#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 16:40:11 2024

@author: m.wehrens
"""

# mydirectory='/Volumes/workdrive_m.wehrens_hubrecht/reference_genomes/GRCh38.p14/'
# myfilepath='GCF_000001405.40_GRCh38.p14_genomic.fna'

mydirectory='/Users/m.wehrens/Documents/git_repos/Side_Projects/Editor_Statistics_Thomas/genome_processing/'
myfilepath='mytest.fasta'

file1=None
c=0
with open(mydirectory+myfilepath) as f:
    for line in f:
        
        
        if (line[0]=='>'):
            
            if (file1 != None):
                file1.close()
                
            print(line)
            
            outputfile=line.split(' ')[0][1:]
            file1 = open(mydirectory+outputfile+".fna", "w") 
            
        
        file1.write(line)
        
        #c=c+1
        #if (c==100):
        #    break
        