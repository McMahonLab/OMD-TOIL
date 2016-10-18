###############################################################################
# trim
# Copyright (c) 2016, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Trim reads using Sickle.
# https://github.com/najoshi/sickle
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
inputFolder = '../../archivalData'
outputFolder = '../../archivalData/trimmed'

# Check that the new output directory exists and create if it doesn't
if not os.path.exists(outputFolder):
        print "Creating output directory\n"
        os.makedirs(outputFolder)
        
#%%#############################################################################
### Step 1 - Read in samples to be trimmed
################################################################################
# Read in list of samples from all three folders within 'inputFolder'
sampleList = []
for path, subdirs, files in os.walk(inputFolder):
    for name in files:
        if name.endswith('fastq'):
            sampleList.append(os.path.join(path, name))

# Because they're paired-ends, conver to a list of pairs
# Borrowed from stack-exchange
pairedEndList = zip(sampleList[::2], sampleList[1::2])

#%%#############################################################################
### Step 2 - Trim the metagenomes using sickle
################################################################################

for pairList in pairedEndList:
    leftSampleString = pairList[0].split('/')[-1]
    rightSampleString = pairList[1].split('/')[-1]

    sampleName = re.search('([^_]+)', leftSampleString).group()
    subprocess.call(['sickle', 'pe',
                      '-f', pairList[0],
                      '-r', pairList[1],
                      '-o', outputFolder+'/'+leftSampleString,
                      '-p', outputFolder+'/'+rightSampleString,
                      '-s', outputFolder+'/'+sampleName+'_singles.fastq',
                      '-t', 'sanger'
                      ])
