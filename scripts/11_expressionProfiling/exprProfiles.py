###############################################################################
# exprProfiles.py
# Copyright (c) 2016, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Code description.
################################################################################

#%%#############################################################################
### Import packages
################################################################################

import os
import pandas as pd

#%%#############################################################################
### Define folder structure
################################################################################

sampleFolder = '../../data/rawData/'
normFolder = '../../data/derivedData/mapping/uncompetitive/RPKM'
profileFolder = '../../data/derivedData/mapping/uncompetitive/profiles'
profileInfile = profileFolder+'/profiles.in'
profileOutfile = profileFolder+'/profiles.out'

#%%#############################################################################
### Create list of (genome, gene) pairs to process
################################################################################

# Import the list of (genome, gene) pairs to process
masterExprMatrix = pd.read_csv(profileInfile, index_col=[0,1], sep=',')

# Read in list of MTs
mtList = []
for mt in os.listdir(sampleFolder):
    if mt.startswith('.'):
        next
    else:
       mtList.append(mt)
       
# Add additonal columns for each MT
masterExprMatrix = pd.concat([masterExprMatrix, pd.DataFrame(0, index=masterExprMatrix.index, columns=mtList)], axis=1)

#%%#############################################################################
### Extract expression profiles for each (genome, gene) pair
################################################################################

for index in masterExprMatrix.index:
    genome = index[0]
    locus = index[1]
    # Read in the RPKM file for the genome
    genomeExprMatrix = pd.read_csv(normFolder+'/'+str(genome)+'.RPKM.out', index_col=[0], sep=',')
    masterExprMatrix.loc[(genome, locus),mtList] = genomeExprMatrix.loc[locus][2:]
    
masterExprMatrix.to_csv(profileOutfile)