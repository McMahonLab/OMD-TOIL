OMD-TOIL: Quality Control, Trimming, and Filtering
===
Copyright (c) 2016, McMahon Lab  
URL: [https://mcmahonlab.wisc.edu/](https://mcmahonlab.wisc.edu/)  
URL: [https://github.com/McMahonLab/](https://github.com/McMahonLab/)  
All rights reserved.

Documentation and analysis by Joshua J. Hamilton  
URL: [https://github.com/joshamilton/](https://github.com/joshamilton/)  

Overview
--
This document outlines procedures for trimming and filtering metagenomic reads.

Prerequisites
--
Ensure the following software is installed:  

* [Sickle](https://github.com/najoshi/sickle) - Determines quality and length thresholds for trimming, and trims reads.
Ensure the following data are available:  

* `archivalData/MGs/sample_index_L00[0-9]_R[1-2]_001.fastq`, where `sample` is the sample ID, `index` is the barcode used for indexing, [0-9] is a digit referring to the lane number, and [1-2] are the unpaired ends.

* `archivalData/MTs_no_rRNA/sample_index_L00[0-9]_R[1-2]_001.fastq`, where `sample` is the sample ID, `index` is the barcode used for indexing, [0-9] is a digit referring to the lane number, and [1-2] are the unpaired ends.

* `archivalData/MTs_rRNA/sample_index_L00[0-9]_R[1-2]_001.fastq`, where `sample` is the sample ID, `index` is the barcode used for indexing, [0-9] is a digit referring to the lane number, and [1-2] are the unpaired ends.

Each of these folders contains the reads associated with

* MGs - metagenomic reads

* MTs_no_rRNA - metatranscriptomics reads for which rRNA depletion was performed prior to sequencing (e.g., these MTs contain no RNA)

* MTs_rRNA - metatranscriptomics reads for which rRNA depletion was NOT performed prior to sequencing (e.g., these MTs contain RNA)

Read Trimming
--

Reads are trimmed using Sickle using default parameters. For convenience, the script `trim.py` will run the following command on each set of unpaired ends:

    sickle pe
    -f archivalData/folder/sample_index_L00[0-9]_R1_001.fastq
    -r archivalData/folder/sample_index_L00[0-9]_R2_001.fastq
    -o archivalData/folder/trimmed/sample_index_L00[0-9]_R1_001.fastq
    -p archivalData/folder/trimmed/sample_index_L00[0-9]_R2_001.fastq
    -s archivalData/folder/trimmed/sample_singles.fastq
    -t sanger

where `folder` is one of the folders described above.
