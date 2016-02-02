OMD-TOIL: Uncompetitive Mapping of Reads
===
Copyright (c) 2016, McMahon Lab  
URL: [https://mcmahonlab.wisc.edu/](https://mcmahonlab.wisc.edu/)  
URL: [https://github.com/McMahonLab/](https://github.com/McMahonLab/)  
All rights reserved.

Documentation and analysis by Joshua J. Hamilton  
URL: [https://github.com/joshamilton/](https://github.com/joshamilton/)  

Mapping of Reads: Uncompetitive
--
This section describes the steps to map metatranscriptomic reads to the reference genomes using BBMap. __This section describes uncompetitive mapping. In this approach, each genome is processed individually, and reads may may to multiple genomes.__

For convenience, the script `uncompReadMapping.pl` will execute the pipeline with a single command. The script takes as input the directories described below, and maps each metatranscriptome to each reference genome with the following command:

    bbmap.sh t=numThreads in=sample_non_rRNA.fastq outm=sample-genome.sam ref=genome.fna nodisk sam=1.3

The script is called as follows:

    `perl uncompReadMapping.pl`

The script performs mapping in parallel, with `numThreads` allocated to each mapping run, parallelized across `numChildren` child processes.

__Note:__ The script is currently hard-coded to use the folder structure described in this repo. The script also specifies use of 1 thread per run, parallelized across 25 processors. For some reason, multithreading with BBMap causes the server to crash.

The inputs to the wrapper function are as follows:

| Argument | Description  |
|---|---|
| mtFolder | location of the metatranscriptome reads to be mapped |
| genomeFolder | location of the reference genomes |
| mapFolder | desired output location of the mapped reads |
| numThreads | number of processors to use for mapping |
| numChildren | number of parallel mapping runs |

Counting of Mapped Reads: Uncompetitive
--

This section describes the steps to count the metatranscriptomic reads which mapped to each gene in a reference genome using BBMap. For convenience, the script `uncompReadCountsFEATURE.pl` will execute the pipeline with a single command. The script takes as input the directories described below, and counts reads mapped to each gene using [htseq-count](http://www-huber.embl.de/HTSeq/doc/count.html#count) for each (metatranscriptome, genome) pair:

    `/usr/local/bin/htseq-count -f bam -r pos -s no -a 0 -t FEATURE -i locus_tag -m intersection-strict -o outputFolder/mt-genome.sam mapFolder/mt-genome.sorted.bam gffFolder/genome.gff > outputFolder/mt-genome.out`

where variables are defined as follows:

| Argument | Description  |
|---|---|
| outputFolder | location to store read counts for each genome |
| gffFolder | location of the gff files for the reference genomes |
| mapFolder | location of the mapped reads |
| mt | the current metatranscriptome|
| genome | the current genome |

and flags are defined as follows:

| Flag | Description  |
|---|---|
| -f bam | mapped reads are stored in a `bam` file |
| -r pos | `bam` file has been sorted by position along the genome |
| -s no | reads are not strand-specific (e.g., they are from cDNA and can match to the same or opposite strand as the gene) |
| -t FEATURE | type of feature to look for, options available in the the `gff` file |
| -i locus_tag | label to apply to each gene, options available in the the `gff` file |
| -m intersection-strict | mapped reads must fully lie within a single gene |
| -o | location of output `sam` file|
| > | location of file with read counts |

Additional info about these flags is available in the `htseq-count` [documentation](http://www-huber.embl.de/HTSeq/doc/count.html#count) and the `GFF` [file specification](http://gmod.org/wiki/GFF2).

The output of the script is a tab-delimited `.FEATURE.out` file giving the locus tag and number of mapped reads.

The script is called as follows:

    `perl uncompReadCountsFEATURE.pl`

The inputs to the wrapper function are as follows:

| Argument | Description  |
|---|---|
| genomeFolder | location of the reference genomes |
| gffFolder | location of the gff files for the reference genomes |
| mtFolder | location of the metatranscriptome reads |
| mapFolder | location of the mapped reads |
| outputFolder | location to store read counts for each genome |

__Note:__ The script is currently hard-coded to use the folder structure described in this repo.

Processing of Mapped Reads: Uncompetitive
--

The python script `percentReadsPerGenome` aggregates these results into the following file:

* `countFolder/percentReadsPerGenome.csv`

    For each (genome, metatranscriptome) pair, this file reports the ratio of:

    `total number of reads which map to a coding sequence in the genome`

    to

    `total number of reads which map to a coding sequence in any genome`

The python script `04computeRPKMs` aggregates these results into the following files:

* `normFolder/genome.counts.out`

    For each genome, this file reports the total number of reads from each metatranscriptome which map to each gene locus of the genome, expressed as the total number of counts. Users can then normalize these count data however they like.

    The file is a comma-separated file, with columns as follows:

| Columns | Description  |
|---|---|
| Locus Tag | locus tag of the gene |
| IMG Gene ID | IMG OID for the gene. Can be used in IMG's "Gene Search" to find more about the gene. |
| Product | A written description of the gene product, such as `actinorhodopsin` or `ammonia permease`. |
| ME150256 | Count of reads which map to the metatranscriptome, expressed on a per-gene basis. |

  The script then converts the total count data to normalized count data, on an RPKM basis:

 [RPKM](http://www.nature.com/nmeth/journal/v5/n7/abs/nmeth.1226.html) stands for 'Reads per Kilobase of Transcript per Million Mapped Reads' and was calculated as follows:

  * __Reads__. As described above, `htseq-count` was used to count the number of reads from the metatranscriptome which map to each gene locus. The mapping was uncompetitive, so a read could map to multiple genomes. The counting was strict, so within a genome a read can only map to one gene.

  * __Kilobase of Transcript__. The start and stop sites of each gene are included in the gff files. I computed the transcript length from the difference between these two sites.

  * __Million mapped reads__. I chose to include _only_ those reads which do not map to either rRNA or the internal standard. Pancho used [SortMeRNA](http://bioinfo.lifl.fr/RNA/sortmerna/) to remove reads which map to rRNA, and I manually mapped reads to the internal standard. The file `countsPerFeature` contains the total reads remaining after the rRNA-removal step (column 'Total Reads') as well as the number which map to the internal standard ('Int Std'). I used the difference between these two values, realizing that some reads still map to rRNA (column 'rRNA').

The output is a file `normFolder/genome.RPKM.out` with columns the same as above.

Cleanup
--

By default, `BBMap` generates SAM files for the aligned (mapped) reads. These files can be converted to BAM files to save space. The conversion process is very slow, so I recommend waiting until you are done with the pipeline to convert the files. For convenience, the script `samToBam.pl` will perform the conversion for you:

      `samtools view -b  -S -o sample-genome.bam sample-genome.sam`

      `samtools sort -o sample-genome.sorted.bam -O bam -T /temp/sample-genome sample-genome.bam`

      `samtools index sample-genome.sorted.bam`

      `rm sample.sam`

      `rm sample.bam`

The script is called as follows:

`perl samToBam.pl`
