###############################################################################
# processReadCounts.py
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
genomeFolder = '../../../rawData/refGenomes/fna'
gffFolder = '../../../rawData/refGenomes/gff'
metadataFolder = '../../../metadata'
sampleFolder = '../../../rawData/'
countFolder = '../../../derivedData/mapping/uncompetitive/readCounts'
normFolder = '../../../derivedData/mapping/uncompetitive/RPKM'
stdName = 'pFN18A_DNA_transcript'

#%%#############################################################################
### Step 0 - Populate MT read frame. Create empty dataframes.
################################################################################

# Read in list of MTs
mtList = []
for mt in os.listdir(pathToMTs):
    if mt.startswith('.'):
        next
    else:
       mtList.append(mt)


# Read in list of genomes. Ignore internal standard genome.
genomeList = []
for genome in os.listdir(pathToGenomes):
    if stdName in genome or genome.startswith('.'):
        next
    else:
       genomeList.append(genome)

genomeList = [genome.replace('.fna', '') for genome in genomeList]

# Create dataframe containing total read counts
mtReads = pd.read_csv(pathToMetadata+'/totalReads.csv', index_col=0)

# Add additional empty colums
mtReads = pd.concat([mtReads, pd.DataFrame(0, index=mtList, columns=['CDS', 'RNA', 'rRNA', 'tRNA', 'Int Std'])], axis=1)

# Create empty dataframe for genome read counts
alignedMatrix = pd.DataFrame(0, index=genomeList, columns=mtReads.index)

#%%#############################################################################
### Step 1 - Count reads which map to each genome. Perform counts for each
### feature type: CDS, rRNA, tRNA, RNA. For the CDS feature type, also
### calculate as the percent of total CDS reads from the metatranscriptome.
################################################################################
for MT in mtList:
    for genome in genomeList:
    # Read in the .CDS.out file and find the total read count
    # Not all CDS file will exist, so employ a check. Create an empty DF if the file doesn't exist.
        if os.path.isfile(pathToReads+'/'+MT+'-'+genome+'.CDS.out'):
            genomeReadsCDS = pd.read_csv(pathToReads+'/'+MT+'-'+genome+'.CDS.out', index_col=0, sep='\t', header=None)
            totalReadsCDS = genomeReadsCDS.sum()[1]
        else:
            genomeReadsCDS = []
            totalReadsCDS = 0

        if os.path.isfile(pathToReads+'/'+MT+'-'+genome+'.rRNA.out'):
            genomeReadsrRNA = pd.read_csv(pathToReads+'/'+MT+'-'+genome+'.rRNA.out', index_col=0, sep='\t', header=None)
            totalReadsrRNA = genomeReadsrRNA.sum()[1]
        else:
            genomeReadsrRNA = []
            totalReadsrRNA = 0

        if os.path.isfile(pathToReads+'/'+MT+'-'+genome+'.tRNA.out'):
            genomeReadstRNA = pd.read_csv(pathToReads+'/'+MT+'-'+genome+'.tRNA.out', index_col=0, sep='\t', header=None)
            totalReadstRNA = genomeReadstRNA.sum()[1]
        else:
            genomeReadstRNA = []
            totalReadstRNA = 0

        if os.path.isfile(pathToReads+'/'+MT+'-'+genome+'.RNA.out'):
            genomeReadsRNA = pd.read_csv(pathToReads+'/'+MT+'-'+genome+'.RNA.out', index_col=0, sep='\t', header=None)
            totalReadsRNA = genomeReadsRNA.sum()[1]
        else:
            genomeReadsRNA = []
            totalReadsRNA = 0

# Add this info to the DF of alignment counts and update the count of total
        # CDS counts for the transcriptome
        alignedMatrix.loc[genome, MT] = totalReadsCDS
        mtReads.loc[MT]['CDS'] = mtReads.loc[MT]['CDS'] + totalReadsCDS

        mtReads.loc[MT]['rRNA'] = mtReads.loc[MT]['rRNA'] + totalReadsrRNA

        mtReads.loc[MT]['tRNA'] = mtReads.loc[MT]['tRNA'] + totalReadstRNA

        mtReads.loc[MT]['RNA'] = mtReads.loc[MT]['RNA'] + totalReadsRNA

# Read in list of counts to the internal standard
for MT in mtList:
    genomeReadsSTD = pd.read_csv(pathToReads+'/'+MT+'-'+stdName+'.CDS.out', index_col=0, sep='\t', header=None)
    totalReadsSTD = genomeReadsSTD.sum()[1]
    mtReads.loc[MT]['Int Std'] = mtReads.loc[MT]['Int Std'] + totalReadsSTD

# Normalize and convert to a percent - coding sequences (CDS) only
for MT in mtList:
    alignedMatrix.loc[:, MT] = (alignedMatrix.loc[:, MT] / mtReads.loc[MT]['CDS']) * 100

# Write to CSV file
alignedMatrix.to_csv(pathToOutput+'/percentReadsPerGenome.csv', sep=',')
mtReads.to_csv(pathToOutput+'/countsPerFeature.csv', sep=',')

#%%#############################################################################
### Step 2 - Construct normalized read counts. Normalize to RPKM, reads per
### kilobase of sequence per million mapped reads. For mapped reads, consider
### only reads which don't map to the standard.
################################################################################

for genome in genomeList:
# Create an empty dataframe with the desired columns
    genomeRPKM = pd.DataFrame(columns=['Locus Tag', 'IMG Gene ID', 'Product', 'Gene Length'])

# Read in the GFF file. The file needs to be split along both tab and semi-
# colon characters. Because the number of fields will vary depending on the
# entries in the attributes column, the file cannot be directly read into a
# dataframe.
    myFile = open(pathToGFF+'/'+genome+'.gff')
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
        # Create a new column corresponding to the MT
#        genomeRPKM[mt] = 0

# Read in the feature.out file
    # Not all CDS file will exist, so employ a check. Create an empty DF if the file doesn't exist.
        if os.path.isfile(pathToReads+'/'+MT+'-'+genome+'.CDS.out'):
            genomeReadsCDS = pd.read_csv(pathToReads+'/'+MT+'-'+genome+'.CDS.out', index_col=0, sep='\t', header=None)
        else:
            genomeReadsCDS = []

        if os.path.isfile(pathToReads+'/'+MT+'-'+genome+'.rRNA.out'):
            genomeReadsrRNA = pd.read_csv(pathToReads+'/'+MT+'-'+genome+'.rRNA.out', index_col=0, sep='\t', header=None)
        else:
            genomeReadsrRNA = []

        if os.path.isfile(pathToReads+'/'+MT+'-'+genome+'.tRNA.out'):
            genomeReadstRNA = pd.read_csv(pathToReads+'/'+MT+'-'+genome+'.tRNA.out', index_col=0, sep='\t', header=None)
        else:
            genomeReadstRNA = []

        if os.path.isfile(pathToReads+'/'+MT+'-'+genome+'.RNA.out'):
            genomeReadsRNA = pd.read_csv(pathToReads+'/'+MT+'-'+genome+'.RNA.out', index_col=0, sep='\t', header=None)
        else:
            genomeReadsRNA = []

        # Merge into a single genomeReads DF and rename the column with the MT name
        concatList = [genomeReadsCDS, genomeReadsRNA, genomeReadsrRNA, genomeReadstRNA]
        concatList = [x for x in concatList if type(x) == pd.core.frame.DataFrame]
        genomeReads = pd.concat(concatList)
        genomeReads.columns = [MT]

        # Perform a left join with the RPKM matrix
        genomeRPKM = genomeRPKM.join(genomeReads, how='left')

        # Convert to RPKM
        # RPKM stands for 'Read per Kilobase of Transcript per Million Mapped Reads'
        # Kilobse of transcript is given by: K = genomeRPKM[Length] / 1000
        # Million mapped reads is given by: M = (mtReads[Total Reads] - mtReads[Int Std]) / 1000000
        # Therefore RPKM = (genomeRPKM[MT] / M) / K
        M = (mtReads['Reads'] - mtReads['Int Std']) / 1000000
        genomeRPKM[MT] = (genomeRPKM[MT] / M[MT]) / (genomeRPKM['Gene Length'] / 1000)

    # Drop the 'Gene Length' column and write to file
    genomeRPKM = genomeRPKM.drop('Gene Length',1)
    genomeRPKM.to_csv(pathToOutput+'/'+genome+'.RPKM.out', sep=',')
