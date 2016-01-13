OMD-TOIL: Mapping Metatranscriptomic Reads
===
Copyright (c) 2016, McMahon Lab  
URL: [https://mcmahonlab.wisc.edu/](https://mcmahonlab.wisc.edu/)  
URL: [https://github.com/McMahonLab/](https://github.com/McMahonLab/)  
All rights reserved.

Documentation and analysis by Joshua J. Hamilton  
URL: [https://github.com/joshamilton/](https://github.com/joshamilton/)  

Overview
--
This document outlines procedures for mapping metatranscriptomic reads to the McMahon lab's collection of reference genomes. The analysis occurred in three parts, each with its own documentation:
  * [Validation of GFF files](10a_gffValidation/README.md)  
  * [Uncompetitive mapping of reads](10b_uncompetitive/README.md)    
  * [Competitive mapping of reads](10c_competitive/README.md)  

This document gives a high-level overview of the mapping process. Refer to individual READMEs for each step for specific details.

Prerequisites
--
Ensure the following software is installed:  

* [Burrows-Wheeler Aligner](http://bio-bwa.sourceforge.net/) (BWA) - As of 2015-11-09. BWA version 0.7.5a-r405 is on Zissou and accessible by typing `bwa` in the command line.

* [Genome Tools](http://genometools.org/pub/) - This will need to be installed in your `home` folder (assuming you cloned your fork of this repo to your `home` folder). Otherwise, Genome Tools will need to be accessible from wherever you cloned your fork.

* [samtools](http://www.htslib.org/download/) - This will need to be installed in your `home` folder  (assuming you cloned your fork of this repo to your `home` folder). Otherwise, samtools will need to be accessible from wherever you cloned your fork.

* [HTSeq](http://www-huber.embl.de/HTSeq/doc/overview.html) - As of 2015-11-17, this Python package is installed on Zissou.

* [Python](https://www.python.org/) with the following packages
    * [Pandas](http://pandas.pydata.org/) - As of 2015-12-16, this Python package is installed on Zissou.

Ensure the following data are available:  

* `../../rawData/sample/rRNA_removal/sample_non_rRNA.fastq`, where `sample` is the unique identifier for the sample (e.g., ME150256). One such file should exist for each set of metatranscriptomic reads you which to map. The path `../../rawData/` is referred throughout this documentation as `mtFolder`.

* `../../rawData/refGenomes/fna` -  a folder containing all reference genomes to be mapped. In FASTA nucleotide format (`file.fna`). __Note__: The FASTA nucleotide sequence for the internal standard should also be included here. This folder is referred throughout this documentation as `genomeFolder`.

* `../../rawData/refGenomes/gff` -  a folder containing all GFF files of all reference genomes to be mapped. In FASTA nucleotide format. __Note__: A dummy GFF file for the internal standard should also be included here. This folder is referred throughout this documentation as `gffFolder`.

Ensure the following folders are accessible, as data will be stored in these folders:

* `../../derivedData/mapping/uncompetitive/bamFiles` - results of the mapping will go here. Referred to as `mapFolder`.

* `../../derivedData/mapping/uncompetitive/readCounts` - total read counts will go here. Referred to as `countFolder`.

* `../../derivedData/mapping/uncompetitive/RPKM` - normalized read counts will go here. Referred to as `normFolder`.

* `../../derivedData/mapping/competitive/bamFiles` - results of the mapping will go here. Referred to as `mapFolder`.

* `../../derivedData/mapping/competitive/readCounts` - total read counts will go here. Referred to as `countFolder`.

* `../../derivedData/mapping/competitive/RPKM` - normalized read counts will go here. Referred to as `normFolder`.

Data Processing
--
In order for to properly count mapped reads, the GFF files must be in the proper format. The script `gffValidator` will analyze all GFF files in the `gffFolder` folder, and generate a file `gffFolder/all.out` describing errors in the GFF files. The GFF files must then be manually corrected.

Uncompetitive Mapping of Reads
--

### Mapping
This section describes the steps to map metatranscriptomic reads to the reference genomes using BWA. __This section describes uncompetitive mapping. In this approach, each genome is processed individually, and reads may may to multiple genomes.__ For convenience, the script `uncompReadMapping.pl` will execute the pipeline with a single command. The script maps each metatranscriptome to each reference genome.

### Counting
This section describes the steps to count the metatranscriptomic reads which mapped to each gene in a reference genome. For convenience, the script `uncompReadCountsFEATURE.pl` will execute the pipeline with a single command. The script counts reads mapped to each gene using [htseq-count](http://www-huber.embl.de/HTSeq/doc/count.html#count) for each (metatranscriptome, genome) pair.

### Normalization
This section describe the process for normalizing the read counts, as well as procedures for producing some other statistics. For each genome, the script `processUncompReadCounts.py` computes the total number of reads from each metatranscriptome which map to each gene locus of the genome, expressed on an RPKM basis.
