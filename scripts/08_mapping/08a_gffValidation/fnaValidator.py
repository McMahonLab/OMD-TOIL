###############################################################################
# fnaValidator.py
# Copyright (c) 2015, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Read in fna files and update header lines.
################################################################################

#%%#############################################################################
### Import packages
################################################################################

from Bio import SeqIO
import os
import pandas as pd
import re

#%%#############################################################################
### Static folder structure
################################################################################
# Define fixed input and output files
genomeFolder = '../../../data/refGenomes/fna'
stdName = 'pFN18A_DNA_transcript'

#%%#############################################################################
### Step 0 - Read in directory lists and files needed for normalization
################################################################################

# Read in list of genomes. Ignore internal standard genome.
genomeList = []
for genome in os.listdir(genomeFolder):
    if stdName in genome or 'merged' in genome or genome.startswith('.'):
        next
    elif genome.endswith('.fna'):
       genomeList.append(genome)

genomeList = [genome.replace('.fna', '') for genome in genomeList]

#%%#############################################################################
### Step 1 - Loop over each fna file and update the header information.
################################################################################

for genome in genomeList:
    inFile = open(genomeFolder+'/'+genome+'.fna', 'r')
    outFile = open(genomeFolder+'/'+genome+'.new.fna', 'w')

    for record in SeqIO.parse(inFile, 'fasta'):
        record.description = record.id
        SeqIO.write(record, outFile, 'fasta')
    inFile.close()
    outFile.close()    

# Delete the old file and rename the new one
    os.remove(genomeFolder+'/'+genome+'.fna')
    os.rename(genomeFolder+'/'+genome+'.new.fna', genomeFolder+'/'+genome+'.fna')