###############################################################################
# processUncompReadCounts
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
countFolder = '../../../data/derivedData/mapping/uncompetitive/readCounts'
normFolder = '../../../data/derivedData/mapping/uncompetitive/RPKM'
stdName = 'pFN18A_DNA_transcript'

# Check that the new output directory exists and create if it doesn't
if not os.path.exists(normFolder):
        print "Creating output directory\n"
        os.makedirs(normFolder)
    
#%%#############################################################################
### Step 0 - Read in MT and genome lists. Create DF to be used for normalization.
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
mtReads = pd.read_csv(metadataFolder+'/totalReads.csv', index_col=0)

# Add additional empty colums
mtReads = pd.concat([mtReads, pd.DataFrame(0, index=mtList, columns=['Int Std', 'NormFact'])], axis=1)    
    
#%%#############################################################################
### Step 1 - Count reads that map to the internal standard
################################################################################

# Read in list of counts to the internal standard
for MT in mtList:
    genomeReadsStd = pd.read_csv(countFolder+'/'+MT+'-'+stdName+'.CDS.out', index_col=0, sep='\t', header=None)
    genomeReadsStd = genomeReadsStd.ix[:-5]
    totalReadsStd = genomeReadsStd.sum()[1]
    mtReads.loc[MT,'Int Std'] = totalReadsStd
    mtReads.loc[MT,'NormFact'] = mtReads.loc[MT,'Reads'] - mtReads.loc[MT,'Int Std']
    

# Create empty dataframe for genome read counts
alignedMatrix = pd.DataFrame(0, index=genomeList, columns=mtReads.index)

#%%#############################################################################
### Step 1 - Count reads which map to each genome. Perform counts for each
### feature type: CDS, rRNA, tRNA, RNA. For the CDS feature type, also
### calculate as the percent of total CDS reads from the metatranscriptome.
################################################################################
for MT in mtList:
    for genome in genomeList:
# Not all CDS file will exist and/or have content, so employ a check. Create an empty DF if the file doesn't exist.
    
        if os.path.isfile(countFolder+'/'+MT+'-'+genome+'.CDS.out') and os.path.getsize(countFolder+'/'+MT+'-'+genome+'.CDS.out'):
            genomeReadsCDS = pd.read_csv(countFolder+'/'+MT+'-'+genome+'.CDS.out', index_col=0, sep='\t', header=None)        
            genomeReadsCDS = genomeReadsCDS.ix[:-5]
            totalReadsCDS = genomeReadsCDS.sum()[1]
        else:
            genomeReadsCDS = []
            totalReadsCDS = 0

        if os.path.isfile(countFolder+'/'+MT+'-'+genome+'.rRNA.out') and os.path.getsize(countFolder+'/'+MT+'-'+genome+'.rRNA.out'):
            genomeReadsrRNA = pd.read_csv(countFolder+'/'+MT+'-'+genome+'.rRNA.out', index_col=0, sep='\t', header=None)
            genomeReadsrRNA = genomeReadsrRNA.ix[:-5]
            totalReadsrRNA = genomeReadsrRNA.sum()[1]
        else:
            genomeReadsrRNA = []
            totalReadsrRNA = 0

        if os.path.isfile(countFolder+'/'+MT+'-'+genome+'.tRNA.out') and os.path.getsize(countFolder+'/'+MT+'-'+genome+'.tRNA.out'):
            genomeReadstRNA = pd.read_csv(countFolder+'/'+MT+'-'+genome+'.tRNA.out', index_col=0, sep='\t', header=None)
            genomeReadstRNA = genomeReadstRNA.ix[:-5]
            totalReadstRNA = genomeReadstRNA.sum()[1]
        else:
            genomeReadstRNA = []
            totalReadstRNA = 0

        if os.path.isfile(countFolder+'/'+MT+'-'+genome+'.RNA.out') and os.path.getsize(countFolder+'/'+MT+'-'+genome+'.RNA.out'):
            genomeReadsRNA = pd.read_csv(countFolder+'/'+MT+'-'+genome+'.RNA.out', index_col=0, sep='\t', header=None)
            genomeReadsRNA = genomeReadsRNA.ix[:-5]
            totalReadsRNA = genomeReadsRNA.sum()[1]
        else:
            genomeReadsRNA = []
            totalReadsRNA = 0
            
# Add this info to the DF of alignment counts and update the count of total
        # CDS counts for the transcriptome
        alignedMatrix.loc[genome, MT] = totalReadsCDS + totalReadsRNA + totalReadsrRNA + totalReadstRNA

# Normalize and convert to a percent - coding sequences (CDS) only
for MT in mtList:
    alignedMatrix.loc[:, MT] = (alignedMatrix.loc[:, MT] / mtReads.loc[MT,'NormFact']) * 100

# Write to CSV file
alignedMatrix.to_csv(normFolder+'/percentReadsPerGenome.csv', sep=',')
