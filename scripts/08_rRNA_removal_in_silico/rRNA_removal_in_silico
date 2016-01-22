### Ribosomal RNA (rRNA) removal

Copyright (c) 2016, McMahon Lab
URL: https://mcmahonlab.wisc.edu/
URL: https://github.com/McMahonLab/
All rights reserved.

Documentation and analysis by Francisco Moya-Flores


[SortMeRNA](http://bioinfo.lifl.fr/RNA/sortmerna/) is a biological sequence analysis tool for filtering, mapping and OTU-picking NGS reads. The core algorithm is based on approximate seeds and allows for fast and sensitive analyses of nucleotide sequences. The main application of SortMeRNA is filtering rRNA from metatranscriptomic data.

[SortMeRNA](http://bioinfo.lifl.fr/RNA/sortmerna/) takes as input a file of reads (fasta or fastq format) and one or multiple rRNA database file(s), and sorts apart rRNA and rejected reads into two files specified by the user. Optionally, it can provide high quality local alignments of rRNA reads against the rRNA database.

The following scripts were taken from the [SortMeRNA User Manual 2.0](http://bioinfo.lifl.fr/RNA/sortmerna/code/SortMeRNA-user-manual-v2.0.pdf). The previous four fasta files generated with FLASH were used as input for [SortMeRNA](http://bioinfo.lifl.fr/RNA/sortmerna/).

##### Installation

Refer to [README](https://github.com/biocore/sortmerna/blob/master/README.md) file.

##### Usage

###### Index multiple rRNA databases

The executable indexdb rna indexes an rRNA database.
```
./indexdb_rna --ref db.fasta,db.idx [OPTIONS]
```

There are eight rRNA representative databases provided in the __sortmerna-2.0/rRNA databases__ folder. All databases were derived from the SILVA SSU and LSU databases (release 119) and the RFAM databases using HMMER 3.1b1 and SumaClust v1.0.00. Additionally, the user can index their own database.

#### Example

In the following example, indexdb rna was used with multiple databases. Multiple databases can be indexed simultaneously by passing them as a ‘:’ separated list to --ref (no spaces allowed).
```
>./indexdb_rna --ref ./rRNA_databases/silva-bac-16s-id90.fasta,./index/silva-bac-16s-db:\
./rRNA_databases/silva-bac-23s-id98.fasta,./index/silva-bac-23s-db:\
./rRNA_databases/silva-arc-16s-id95.fasta,./index/silva-arc-16s-db:\
./rRNA_databases/silva-arc-23s-id98.fasta,./index/silva-arc-23s-db:\
./rRNA_databases/silva-euk-18s-id95.fasta,./index/silva-euk-18s-db:\
./rRNA_databases/silva-euk-28s-id98.fasta,./index/silva-euk-28s:\
./rRNA_databases/rfam-5s-database-id98.fasta,./index/rfam-5s-db:\
./rRNA_databases/rfam-5.8s-database-id98.fasta,./index/rfam-5.8s-db

```

###### Filter rRNA reads

The executable sortmerna can filter rRNA reads against an indexed rRNA database.

#### Usage

```
./sortmerna --ref db.fasta,db.idx --reads file.fa --aligned base_name_output [OPTIONS]
```

In the following example, four samples contained in folder AHGJNTBCXX were automatically rRNA filtered with [SortMeRNA](http://bioinfo.lifl.fr/RNA/sortmerna/) through a bash script. Samples are listed below.


MetaTranscriptomes_list_out_extendedFrags.txt

Sample name|File
--|--
ME150263|ME150263.fastq
ME150266|ME150266.fastq
ME150283|ME150283.fastq
ME150290|ME150290.fastq

```
#!/bin/bash

for file in $(</MetaTranscriptomes_list_out_extendedFrags.txt)

do

time ./sortmerna --ref ./rRNA_databases/silva-bac-16s-id90.fasta,./index/silva-bac-16s-db:\
./rRNA_databases/silva-bac-23s-id98.fasta,./index/silva-bac-23s-db:\
./rRNA_databases/silva-arc-16s-id95.fasta,./index/silva-arc-16s-db:\
./rRNA_databases/silva-arc-23s-id98.fasta,./index/silva-arc-23s-db:\
./rRNA_databases/silva-euk-18s-id95.fasta,./index/silva-euk-18s-db:\
./rRNA_databases/silva-euk-28s-id98.fasta,./index/silva-euk-28s:\
./rRNA_databases/rfam-5s-database-id98.fasta,./index/rfam-5s-db:\
./rRNA_databases/rfam-5.8s-database-id98.fasta,./index/rfam-5.8s-db\
 --reads /AHGJNTBCXX/out_extendedFrags/${file}.fastq --sam --num_alignments 1 --fastx --aligned /AHGJNTBCXX/out_extendedFrags/${file}_rRNA --other /AHGJNTBCXX/out_extendedFrags/${file}_non_rRNA --log -v -a 10

done
```
