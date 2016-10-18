OMD-TOIL: rRNA Removal
===
Copyright (c) 2016, McMahon Lab  
URL: [https://mcmahonlab.wisc.edu/](https://mcmahonlab.wisc.edu/)  
URL: [https://github.com/McMahonLab/](https://github.com/McMahonLab/)  
All rights reserved.

Documentation and analysis by Francisco Moya-Flores and Joshua J. Hamilton  
URL: [https://github.com/fmoyaflores/](https://github.com/fmoyaflores/)  
URL: [https://github.com/joshamilton/](https://github.com/joshamilton/)

Prerequisites
--
Ensure the following software is installed:  

* [SortMeRNA](http://bioinfo.lifl.fr/RNA/sortmerna/) - removes rRNA reads from metatranscriptomic sequences

[SortMeRNA](http://bioinfo.lifl.fr/RNA/sortmerna/) takes as input a file of reads (fasta or fastq format) and one or multiple rRNA database file(s), and sorts apart rRNA and rejected reads into two files specified by the user. Optionally, it can provide high quality local alignments of rRNA reads against the rRNA database.

The following scripts were taken from the [SortMeRNA User Manual 2.0](http://bioinfo.lifl.fr/RNA/sortmerna/code/SortMeRNA-user-manual-v2.0.pdf). The previous four fasta files generated with FLASH were used as input for [SortMeRNA](http://bioinfo.lifl.fr/RNA/sortmerna/).

rRNA Removal
--
rRNA reads are removed using SortMeRNA. For convenience, the script `sortMeRNA.py` will run the following commands on each metatranscriptomic sample:

First, index the reference databases for each metatranscriptome sample:

    indexdb_rna --ref pathToSortMeRNA/rRNA_databases/db#.fasta,pathToSortMeRNA/rRNA_databases/db#.idx

where `db#` represents a database. The script identifies all databases present in `pathToSortMeRNA/rRNA_databases` and indexes them.


Second, run SortMeRNA, using the just-indexed databaes to identify and remove rRNA:

    sortMeRNA
    --ref pathToSortMeRNA/rRNA_databases/db#.fasta,pathToSortMeRNA/rRNA_databases/db#.idx:pathToSortMeRNA/rRNA_databases/db#.fasta,pathToSortMeRNA/rRNA_databases/db#.idx
    --reads archivalData/merged/sample/out.extendedFrags.fastq
    --aligned archivalData/sortMeRNA/sample_rRNA
    --other archivalData/sortMeRNA/sample_nor_rRNA
    --fastx
    -a 24
    -v

where the argument to `--ref` is a `:`-delimited list of `db#.fasta,db#.idx` pairs and `-a` is the number of processors to use.

The script also moves the file containing the non-rRNA reads to `rawData` and renames it to `sample.fastq`.
