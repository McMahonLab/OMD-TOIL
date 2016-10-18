###############################################################################
# merge
# Copyright (c) 2016, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Merge reads using FLASH.
# https://ccb.jhu.edu/software/FLASH/
################################################################################

#%%#############################################################################
### Import packages
################################################################################

import os
import re
import subprocess

#%%#############################################################################
### Static folder structure
################################################################################

# Define fixed input and output files
inputFolder = '../../archivalData/trimmed'
outputFolder = '../../archivalData/merged'
moveFolder = '../../rawData'

# Check that the new output directory exists and create if it doesn't
if not os.path.exists(outputFolder):
        print "Creating output directory\n"
        os.makedirs(outputFolder)
        
if not os.path.exists(moveFolder):
        print "Creating output directory\n"
        os.makedirs(moveFolder)
        
#%%#############################################################################
### Step 1 - Read in samples to be trimmed
################################################################################
# Read in list of samples
sampleList = []
for sample in os.listdir(inputFolder):
    if sample.endswith('.fastq') and not '_singles' in sample:
       sampleList.append(sample)

# Because they're paired-ends, conver to a list of pairs
# Borrowed from stack-exchange
pairedEndList = zip(sampleList[::2], sampleList[1::2])

#%%#############################################################################
### Step 2 - Merge the reads using FLASH
################################################################################

# Create a dictionary for the -max--overlap parameter
overlapDict = {
    'ME150256': 150,
    'ME150270': 150,
    'ME150276': 150,
    'ME150286': 150,
    
    'ME150263': 100,
    'ME150266': 100,
    'ME150283': 100,
    'ME150290': 100,

    'ME-150289': 250,
    'TOIL-8-20-15-1': 250,
    'TOIL-8-20-15-2': 250,
}

mgList = ['ME-150289', 'TOIL-8-20-15-1','TOIL-8-20-15-2']

for pairList in pairedEndList:
    sampleName = re.search('([^_]+)', pairList[0]).group()
    subprocess.call(['flash',
                     inputFolder+'/'+pairList[0],
                     inputFolder+'/'+pairList[1],
                     '-d', outputFolder+'/'+sampleName,
                     '-M', str(overlapDict[sampleName])
                      ])
    
    # If the sample is a MG sample, move the output as no further processing is required
    if sampleName in mgList:
        os.rename(outputFolder+'/'+sampleName+'/out.extendedFrags.fastq', moveFolder+'/'+sampleName+'.fastq')
