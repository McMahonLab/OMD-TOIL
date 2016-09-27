OMD-TOIL: Extracting Expression Profiles
===
Copyright (c) 2016, McMahon Lab  
URL: [https://mcmahonlab.wisc.edu/](https://mcmahonlab.wisc.edu/)  
URL: [https://github.com/McMahonLab/](https://github.com/McMahonLab/)  
All rights reserved.

Documentation and analysis by Joshua J. Hamilton  
URL: [https://github.com/joshamilton/](https://github.com/joshamilton/)  

Overview
--
This document outlines procedures for extracting expression profiles for genes through a metatranscriptomic time-series.

Prerequisites
--
Ensure the following software is installed:  

* [Python](https://www.python.org/) with the following packages:
    * [Pandas](http://pandas.pydata.org/)

Ensure the following data are available:  

* `../../rawData/sample`, where `sample` is the unique identifier for the sample (e.g., ME150256). One such folder should exist for each set of metatranscriptomic reads. Referred to as `sampleFolder`.

* `../../derivedData/mapping/uncompetitive/RPKM/genome.RPKM.out` - normalized read counts for all genomes. Referred to as `normFolder`.

* `../../derivedData/mapping/uncompetitive/profiles` - expression profiles will go here. Referred to as `profileFolder`.

Data Input
--
In order for to properly extract expression profiles, gene information must be provided in the proper format. The input should be a file `normFolder/profiles.in`, formatted as a `csv` file with the following columns:

| Column | Description |
|---|---|
| GenomeOID | IMG ID for the genome |
| Locus Tag | locus tag of the gene |
| Optional Columns | optional info, such as a description of the protein product |

Columns beyond the first two are optional, and will be retained in the output.

Extracting Expression Profiles
--

The script `exprProfiles` reads in a list of genome OIDs and locus tags and extracts the RPKM value for that gene across the time-series. The output is a `csv`-formatted file with columns as follows:

| Column | Description |
|---|---|
| GenomeOID | IMG ID for the genome |
| Locus Tag | locus tag of the gene |
| Optional Columns | optional info, such as a description of the protein product |
| Sample | One per MT, RPKM value of the gene in that sample |

The inputs to the wrapper function are as follows:

| Argument | Description  |
|---|---|
| sampleFolder | location of the metatranscriptome reads (reads not required, just sample names) |
| normFolder | location of the RPKM values for each genome |
| profileFolder | folder to store expression profiles |
| profileInfile | file containing locus tags for expression profiling, should be within `profileFolder` |
| profileOutfile | file containing expression profiles for desired locus tags, should be within `profileFolder` |

__Note:__ The script is currently hard-coded to use the folder structure described in this repo, with input and output of `profile.in` and `profile.out`, respectively.
