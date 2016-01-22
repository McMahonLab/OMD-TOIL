### DNA reads: Quality-filtering & Trimming

Copyright (c) 2016, McMahon Lab
URL: https://mcmahonlab.wisc.edu/
URL: https://github.com/McMahonLab/
All rights reserved.

Documentation and analysis by Francisco Moya-Flores

[FLASH](http://ccb.jhu.edu/software/FLASH/) (Fast Length Adjustment of SHort reads) is a very accurate and fast tool to merge paired-end reads that were generated from DNA fragments whose lengths are shorter than twice the length of reads. Merged read pairs result in unpaired longer reads. Longer reads are generally more desired in genome assembly and genome analysis processes. They can also improve transcriptome assembly when FLASH is used to merge RNA-seq data.

##### Installation
Refer to [MANUAL](https://github.com/genome-vendor/FLASH/blob/master/MANUAL) file

##### Usage
flash <mates1.fastq> <mates2.fastq> [-m minOverlap] [-M maxOverlap] [-x mismatchRatio]
[-p phredOffset] [-o prefixOfOutputFiles] [-d pathToDirectoryForOutputFiles]
[-f averageFragment Length] [-s standardDeviationOfFragments] [-r averageReadLength]
[-h displayHelp]

mates1.fastq and mates2.fastq are fastq files of paired-end reads from the short fragment library (with insert size less than twice the length of reads). The corresponding mates should be in the same order in both files.

##### Example

In the following example, four samples were automatically quality filtered and merged with FLASH through a bash script. Samples are listed below.


  Sample name|File (Fw(1) and Rv(2))
--|--
  ME150263|ME150263_GAGTGG_L001_R[1/2]_001.fastq
  ME150266|ME150266_ATTCCT_L001_R[1/2]_001.fastq
  ME150283|ME150283_ACTGAT_L001_R[1/2]_001.fastq
  ME150290|ME150290_CGTACG_L001_R[1/2]_001.fastq

```
#!/bin/bash

for file in $(</MetaTranscriptomes_list.txt)

do

flash /${file}_R1_001.fastq /${file}_R2_001.fastq --interleaved-output -d /AHGJNTBCXX/${file}

done
```
