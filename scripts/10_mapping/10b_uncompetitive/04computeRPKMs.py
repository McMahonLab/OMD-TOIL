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

from joblib import Parallel, delayed  
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

# Define number of processors for multi-processing
numCores = 30
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
### Step 2 - Rearrange mapped read counts.
################################################################################

def rearrangeReads(genome):
    
# Create an empty dataframe with the desired columns
    genomeRPKM = pd.DataFrame(columns=['Locus Tag', 'IMG Gene ID', 'Product', 'Gene Length'])

# Read in the GFF file. The file needs to be split along both tab and semi-
# colon characters. Because the number of fields will vary depending on the
# entries in the attributes column, the file cannot be directly read into a
# dataframe.
    myFile = open(gffFolder+'/'+genome+'.gff')
    for line in myFile:
        line = line.rstrip()
        if line == '##gff-version 3':
            next
        else:
# Split along the appropriate delimiters
            gffArray = re.split('\t|;', line)
# Assign elements to their proper location in the dataframe
            if len(gffArray) >= 11:
                genomeRPKM = genomeRPKM.append({'Locus Tag': gffArray[9].split('=')[1],
                                                'IMG Gene ID': gffArray[8].split('=')[1],
                                                'Product': gffArray[10].split('=')[1],
                                                'Gene Length': int(gffArray[4]) - int(gffArray[3]) + 1 },
                                                ignore_index = True)
            else:
                genomeRPKM = genomeRPKM.append({'Locus Tag': gffArray[9].split('=')[1],
                                                'IMG Gene ID': gffArray[8].split('=')[1],
                                                'Product': 'None Provided',
                                                'Gene Length': int(gffArray[4]) - int(gffArray[3]) + 1 },
                                                ignore_index = True)
    myFile.close()

# Reindex based on the locus tag for faster processing of read counts
    genomeRPKM = genomeRPKM.set_index('Locus Tag')

# Now read in the read counts from each genome-MT.feature.out file and add to the DF
    for MT in mtList:

# Read in the feature.out file and drop the unncessary rows
# Not all CDS file will exist and/or have content, so employ a check. Create an empty DF if the file doesn't exist.    
        if os.path.isfile(countFolder+'/'+MT+'-'+genome+'.CDS.out') and os.path.getsize(countFolder+'/'+MT+'-'+genome+'.CDS.out'):
            genomeReadsCDS = pd.read_csv(countFolder+'/'+MT+'-'+genome+'.CDS.out', index_col=0, sep='\t', header=None)        
            genomeReadsCDS = genomeReadsCDS.ix[:-5]
        else:
            genomeReadsCDS = pd.DataFrame(columns=['1'])

        if os.path.isfile(countFolder+'/'+MT+'-'+genome+'.rRNA.out') and os.path.getsize(countFolder+'/'+MT+'-'+genome+'.CDS.out'):
            genomeReadsrRNA = pd.read_csv(countFolder+'/'+MT+'-'+genome+'.rRNA.out', index_col=0, sep='\t', header=None)
            genomeReadsrRNA = genomeReadsrRNA.ix[:-5]
        else:
            genomeReadsrRNA = pd.DataFrame(columns=['1'])

        if os.path.isfile(countFolder+'/'+MT+'-'+genome+'.tRNA.out') and os.path.getsize(countFolder+'/'+MT+'-'+genome+'.CDS.out'):
            genomeReadstRNA = pd.read_csv(countFolder+'/'+MT+'-'+genome+'.tRNA.out', index_col=0, sep='\t', header=None)
            genomeReadstRNA = genomeReadstRNA.ix[:-5]
        else:
            genomeReadstRNA = pd.DataFrame(columns=['1'])

        if os.path.isfile(countFolder+'/'+MT+'-'+genome+'.RNA.out') and os.path.getsize(countFolder+'/'+MT+'-'+genome+'.CDS.out'):
            genomeReadsRNA = pd.read_csv(countFolder+'/'+MT+'-'+genome+'.RNA.out', index_col=0, sep='\t', header=None)
            genomeReadsRNA = genomeReadsRNA.ix[:-5]
        else:
            genomeReadsRNA = pd.DataFrame(columns=['1'])

        # Merge into a single genomeReads DF and rename the column with the MT name
        genomeReads = pd.concat([genomeReadsCDS, genomeReadsRNA, genomeReadsrRNA, genomeReadstRNA])
        genomeReads.columns = [MT]

        # Perform a left join with the RPKM matrix
        genomeRPKM = genomeRPKM.join(genomeReads, how='left')

    # Write to file
    genomeRPKM.to_csv(normFolder+'/'+genome+'.counts.out', sep=',')

    return

Parallel(n_jobs=numCores)(delayed(rearrangeReads)(genome) for genome in genomeList)

#%%#############################################################################
### Step 2 - Construct normalized read counts. Normalize to RPKM, reads per
### kilobase of sequence per million mapped reads. For mapped reads, consider
### only reads which don't map to the standard.
################################################################################

# Define a function to compute the RPKM for the (genome MT) pair
def computeRPKM(genome):
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

    return

Parallel(n_jobs=numCores)(delayed(computeRPKM)(genome) for genome in genomeList)