OMD-TOIL: Competitive Mapping of Reads
===
Copyright (c) 2016, McMahon Lab  
URL: [https://mcmahonlab.wisc.edu/](https://mcmahonlab.wisc.edu/)  
URL: [https://github.com/McMahonLab/](https://github.com/McMahonLab/)  
All rights reserved.

Documentation and analysis by Joshua J. Hamilton  
URL: [https://github.com/joshamilton/](https://github.com/joshamilton/)  

Mapping of Reads: Competitive
--


Counting of Mapped Reads: Competitive
--


Processing of Mapped Reads: Competitive
--

  The results from the above script `compReadMapping.pl` contain read counts for the entire set of loci found in our reference genome collection. The script `splitCompReads.py` takes these files (`mt.FEATURE.out`) and creates a single file `genome.FEATURE.out` for each (genome, feature) combination in `countFolder`.

The python script `processCompReadCounts.py` aggregates these results into a number of files, as follows:

* `countFolder/percentReadsPerGenome.csv`

    For each (genome, metatranscriptome) pair, this file reports the ratio of:

    `total number of reads which map to a coding sequence in the genome`

    to

    `total number of reads which map to a coding sequence in any genome`

    expressed as a percentage. Because each genome is processed individually, it is possible that some reads in the denominator are being counted more than once. The section below describes a competitive mapping process, which counts only those reads which uniquely map to a single genome.

* `countFolder/countsPerFeature.csv`

    For each metatranscriptome, this file reports the total number of reads (Total Reads), as well as the number of reads which map to different genomic features:  
      * coding sequencs (CDS)  
      * rRNA  
      * tRNA  
      * other RNA (RNA)  
      * the internal standard (Int Std)

* `normFolder/genome.RPKM.out`

    For each genome, this file reports the total number of reads from each metatranscriptome which map to each gene locus of the genome, expressed on an RPKM basis.

    The file is a comma-separated file, with columns as follows:

| Columns | Description  |
|---|---|
| Locus Tag | locus tag of the gene |
| IMG Gene ID | IMG OID for the gene. Can be used in IMG's "Gene Search" to find more about the gene. |
| Product | A written description of the gene product, such as `actinorhodopsin` or `ammonia permease`. |
| ME150256 | Count of reads which map to the metatranscriptome, expressed on an RPKM basis. |

  __Calculating RPKM__. [RPKM](http://www.nature.com/nmeth/journal/v5/n7/abs/nmeth.1226.html) stands for 'Reads per Kilobase of Transcript per Million Mapped Reads' and was calculated as follows:

  * __Reads__. As described above, `htseq-count` was used to count the number of reads from the metatranscriptome which map to each gene locus. The mapping was competitive, so a read could map to only one genome. The counting was strict, so within a genome a read can only map to one gene.

  * __Kilobase of Transcript__. The start and stop sites of each gene are included in the gff files. I computed the transcript length from the difference between these two sites.

  * __Million mapped reads__. I chose to include _only_ those reads which do not map to either rRNA or the internal standard. Pancho used [SortMeRNA](http://bioinfo.lifl.fr/RNA/sortmerna/) to remove reads which map to rRNA, and I manually mapped reads to the internal standard. The file `countsPerFeature` contains the total reads remaining after the rRNA-removal step (column 'Total Reads') as well as the number which map to the internal standard ('Int Std'). I used the difference between these two values, realizing that some reads still map to rRNA (column 'rRNA').
