###############################################################################
# processCompReadCounts
# Copyright (c) 2015, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Process read count data from mapping of OMD-TOIL MT reads to our reference
# genomes.
################################################################################

#%%#############################################################################
### Import packages
################################################################################

import os
import pandas as pd
import re

#%%#############################################################################
### Static folder structure
################################################################################
# Define fixed input and output files
genomeFolder = '../../../data/refGenomes/fna'
gffFolder = '../../../data/refGenomes/gff'
metadataFolder = '../../../metadata'
sampleFolder = '../../../data/rawData/'
countFolder = '../../../data/derivedData/mapping/competitive/readCounts'
normFolder = '../../../data/derivedData/mapping/competitive/RPKM'
stdName = 'pFN18A_DNA_transcript'

# Check that the new output directory exists and create if it doesn't
if not os.path.exists(normFolder):
        print "Creating output directory\n"
        os.makedirs(normFolder)
    
#%%#############################################################################
### Step 0 - Read in directory lists and files needed for normalization
################################################################################

# Read in list of MTs
mtList = []
for mt in os.listdir(sampleFolder):
    if mt.startswith('.'):
        next
    else:
       mtList.append(mt)

# Read in list of genomes. Ignore internal standard genome.
genomeList = []
for genome in os.listdir(genomeFolder):
    if stdName in genome or 'merged' in genome or genome.startswith('.'):
        next
    elif genome.endswith('.fna'):
       genomeList.append(genome)

genomeList = [genome.replace('.fna', '') for genome in genomeList]

# Create dataframe containing total read counts
mtReads = pd.read_csv(normFolder+'/countsPerFeature.csv', index_col=0)

#%%#############################################################################
### Step 1 - Read in each genome.counts.out file and compute the RPKM
################################################################################
for genome in genomeList:
    genomeRPKM = pd.read_csv(normFolder+'/'+genome+'.counts.out', index_col=0)  

    for MT in mtList:
        # Convert to RPKM
        # RPKM stands for 'Read per Kilobase of Transcript per Million Mapped Reads'
        # Kilobse of transcript is given by: K = genomeRPKM[Length] / 1000
        # Million mapped reads is given by: M = (mtReads[Total Reads] - mtReads[Int Std]) / 1000000
        # Therefore RPKM = (genomeRPKM[MT] / M) / K
        M = (mtReads['Reads'] - mtReads['Int Std']) / 1000000
        genomeRPKM[MT] = (genomeRPKM[MT] / M[MT]) / (genomeRPKM['Gene Length'] / 1000)

        # Drop the 'Gene Length' column and write to file
    genomeRPKM = genomeRPKM.drop('Gene Length',1)
    genomeRPKM.to_csv(normFolder+'/'+genome+'.RPKM.out', sep=',')
