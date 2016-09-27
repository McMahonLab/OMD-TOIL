OMD-TOIL: Competitive Mapping of Reads
===
Copyright (c) 2016, McMahon Lab  
URL: [https://mcmahonlab.wisc.edu/](https://mcmahonlab.wisc.edu/)  
URL: [https://github.com/McMahonLab/](https://github.com/McMahonLab/)  
All rights reserved.

Documentation and analysis by Joshua J. Hamilton  
URL: [https://github.com/joshamilton/](https://github.com/joshamilton/)  

Mapping and Counting of Reads: Competitive
--
This section describes the steps to map and count metatranscriptomic reads to the reference genomes. __This section describes competitive mapping, in which reads are mapped to all reference genomes simultaneously. This enables identification of reads which map uniquely to a single genome.__

For convenience, the script `compReadMapping.pl` will execute the pipeline with a single command. The script takes as input the directories described below, and maps each metatranscriptome to each collection of reference genomes. The script also counts reads mapped to each gene using [htseq-count](http://www-huber.embl.de/HTSeq/doc/count.html#count) for each metatranscriptome via the following commands:

1. Concatenate the individual `fasta` and `gff` files. A single file is required for mapping and downstream counting:

    `cat gffFolder/*.gff > gffFolder/merged.out; mv gffFolder/merged.out gffFolder/merged.gff`

    `cat genomeFolder/*.gff > genomeFolder/merged.out; mv genomeFolder/merged.out genomeFolder/merged.fna`

2. Index the reference genomes.

    `bwa index -p genomeFolder/merged -a is genomeFolder/merged.fna`

3. Map the metatranscriptomes to the reference genomes. For each metatranscriptome `sample_non_rRNA.fastq` and reference genome `genome.fna`:

    `bwa mem -t numProcs genomeFolder/merged mtFolder/sample_non_rRNA.fastq > mapFolder/sample.sam`

    where `numProcs` specifies the number of processors. Zissou has 32, please don't use all of them.

5.  Count the reads which uniquely map to each `feature` type (CDS, rRNA, tRNA, other RNA) in our collection of reference genomes:

    `/usr/local/bin/htseq-count -f bam -r pos -s no -a 0 -t feature -i locus_tag -m intersection-strict -o countFolder/mt.feature.sam mapFolder/mt.sorted.bam gffFolder/merged.gff > countFolder/mt.feature.out`

The script is called as follows:

    `perl compReadCounts.pl`

__Note:__ The script is currently hard-coded to use the folder structure described in this repo. The script also specifies use of 30 processors.

The inputs to the wrapper function are as follows:

| Argument | Description  |
|---|---|
| mtFolder | location of the metatranscriptome reads to be mapped |
| genomeFolder | location of the reference genomes |
| gffFolder | location of the gff files for the reference genomes |
| mapFolder | desired output location of the mapped reads |
| countFolder | location to store read counts for each genome |
| numProcs | number of processors to use for mapping |

and flags to `htseq-count` are defined as follows:

| Flag | Description  |
|---|---|
| -f bam | mapped reads are stored in a `bam` file |
| -r pos | `bam` file has been sorted by position along the genome |
| -s no | reads are not strand-specific (e.g., they are from cDNA and can match to the same or opposite strand as the gene) |
| -t feature | type of feature to look for, options available in the the `gff` file |
| -i locus_tag | label to apply to each gene, options available in the the `gff` file |
| -m intersection-strict | mapped reads must fully lie within a single gene |
| -o | location of output `sam` file|
| > | location of file with read counts |

Additional info about these flags is available in the `htseq-count` [documentation](http://www-huber.embl.de/HTSeq/doc/count.html#count) and the `GFF` [file specification](http://gmod.org/wiki/GFF2).

The output of the script is a tab-delimited `.feature.out` file giving the locus tag and number of mapped reads for each metatranscriptome.


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

      * `normFolder/genome.counts.out`

      For each genome, this file reports the total number of reads from each metatranscriptome which map to each gene locus of the genome, expressed as the total number of counts. Users can then normalize these count data however they like.

      The file is a comma-separated file, with columns as follows:

| Columns | Description  |
|---|---|
| Locus Tag | locus tag of the gene |
| IMG Gene ID | IMG OID for the gene. Can be used in IMG's "Gene Search" to find more about the gene. |
| Product | A written description of the gene product, such as `actinorhodopsin` or `ammonia permease`. |
| ME150256 | Count of reads which map to the metatranscriptome, expressed on a per-gene basis. |

__Calculating RPKM__.

The python script `computeRPKMs` will then convert the total count data to normalized count data, on an RPKM basis:

 [RPKM](http://www.nature.com/nmeth/journal/v5/n7/abs/nmeth.1226.html) stands for 'Reads per Kilobase of Transcript per Million Mapped Reads' and was calculated as follows:

 * __Reads__. As described above, `htseq-count` was used to count the number of reads from the metatranscriptome which map to each gene locus. The mapping was uncompetitive, so a read could map to multiple genomes. The counting was strict, so within a genome a read can only map to one gene.

 * __Kilobase of Transcript__. The start and stop sites of each gene are included in the gff files. I computed the transcript length from the difference between these two sites.

 * __Million mapped reads__. I chose to include _only_ those reads which do not map to either rRNA or the internal standard. Pancho used [SortMeRNA](http://bioinfo.lifl.fr/RNA/sortmerna/) to remove reads which map to rRNA, and I manually mapped reads to the internal standard. The file `countsPerFeature` contains the total reads remaining after the rRNA-removal step (column 'Total Reads') as well as the number which map to the internal standard ('Int Std'). I used the difference between these two values, realizing that some reads still map to rRNA (column 'rRNA').

The output is a file `normFolder/genome.RPKM.out`.

Cleanup
  --

  By default, `bwa` generates SAM files for the aligned (mapped) reads. These files can be converted to BAM files to save space. The conversion process is very slow, so I recommend waiting until you are done with the pipeline to convert the files. For convenience, the script `samToBam.pl` will perform the conversion for you:

        `samtools view -b  -S -o sample-genome.bam sample-genome.sam`

        `samtools sort -o sample-genome.sorted.bam -O bam -T /temp/sample-genome sample-genome.bam`

        `samtools index sample-genome.sorted.bam`

        `rm sample.sam`

        `rm sample.bam`

  The script is called as follows:

  `perl samToBam.pl`
