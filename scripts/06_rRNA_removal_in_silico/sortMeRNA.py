###############################################################################
# sortMeRNA
# Copyright (c) 2016, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Trim metagenomic reads using Sickle.
# https://github.com/najoshi/sickle
################################################################################

#%%#############################################################################
### Import packages
################################################################################

import os
import subprocess

#%%#############################################################################
### Static folder structure
################################################################################

# Define fixed input and output files
inputFolder = '../../archivalData/merged'
outputFolder = '../../archivalData/sortMeRNA'
moveFolder = '../../rawData'

refDBFolder = '/Applications/sortmerna-2.1-mac-64-multithread/rRNA_databases'

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

sampleList = ['ME150256', 'ME150270', 'ME150276', 'ME150286',
          'ME150263', 'ME150266', 'ME150283', 'ME150290']

#%%#############################################################################
### Step 2 - Index the SortMeRNA Databases
################################################################################

refDBList = []
for refDB in os.listdir(refDBFolder):
    if refDB.endswith('.fasta'):
       refDBList.append(refDB)

for refDB in refDBList:
    refDBName = refDB.replace('.fasta', '')
    subprocess.call(['indexdb_rna',
                     '--ref', refDBFolder+'/'+refDBName+'.fasta,'+refDBFolder+'/'+refDBName+'.idx'    
    ])

#%%#############################################################################
### Step 3 - Call SortMeRNA
################################################################################

dbArgList = []
for refDB in refDBList:
    refDBName = refDB.replace('.fasta', '')
    dbArgString = refDBFolder+'/'+refDBName+'.fasta,'+refDBFolder+'/'+refDBName+'.idx'
    dbArgList.append(dbArgString)

dbArgString = ':'.join(dbArgList)

for sample in sampleList:
    subprocess.call(['sortMeRNA',
                     '--ref', dbArgString,
                     '--reads', inputFolder+'/'+sample+'/out.extendedFrags.fastq',
                     '--aligned', outputFolder+'/'+sample+'_rRNA',
                     '--other', outputFolder+'/'+sample+'_not_rRNA',
                     '--fastx',
                     '-a', str(3),
                     '-v'
                      ])

    # Move the output    
    os.rename(outputFolder+'/'+sample+'_not_rRNA.fastq', moveFolder+'/'+sample+'.fastq')
