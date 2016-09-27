#########################
rRNA Analysis of OMD-TOIL data
Copyright (c) Katherine McMahon + lab members
Analysis by Alex Linz <amlinz@wisc.edu> on 11-25-15
README last updated by Alex Linz on 1-12-16
#########################

SUMMARY

On Aug 20, the McMahon Lab performed 24 hour sampling for metatranscriptomic sequencing (please see the main repo README document for more information on this project). Four samples were sequenced without rRNA depletion. Greater than 95% of reads in these samples were identified as rRNA.
In order to assess the expression levels of rRNA, we analyzed this rRNA data and classified sequences to get an estimate of the activity levels of freshwater taxa over the 24 hour sampling effort.

Prior to this analysis, sequencing reads had been merged and rRNA had been filtered out using SortmeRNA, using the Silva database as a reference. (Please see the README docs for each of these steps for more information). 
We wanted to classify these rRNA reads and count the number attributed to each taxon in order to estimate expression levels of taxa over diel scales. 
The general approach was to align rRNA reads so that only the V4 region is kept, and then treat the remaining reads as 16S sequencing data, classifying and counting them in mothur.

More specifically, here are the steps in the workflow:
- run SortmeRNA on the non-rRNA files using our customized freshwater database (see McMahon Lab Github page) to make sure we didn't miss any rRNA
- concatenate the results of the SortmeRNA runs (from Silva run and freshwater database run)
- use 16S sequencing data generated as part of the Mendota metagenomes project as a reference for the V4 region (available on JGI's IMG)
- run SortmeRNA on concatenated rRNA files using Mendota iTags as reference
- Convert resulting fastq output to fasta using the fastx-toolkit
- Begin mothur analysis by making group files for each fasta file
- Merge group and fasta files
- Align rRNA V4 regions and trim to only this region (run the alignment step on the Mendota iTAGs reference to determine where to trim the rRNA reads)
- Continue standard McMahon Lab 16S analysis workflow (includes chimera removal, chloroplast removal, OTU clustering at 98% sequence identity, and dual classification with both Greengenes and customized freshwater database as explained in McMahon Lab/16Tax-Ass repo)
- A major deviation from this protocol was the absence of a rarefaction step; instead, proportional data was used to normalize samples
- Analyze resulting OTU table and taxonomic assignments in R using package OTUtable (see McMahonLab/NorthTemperateLakes-MicrobialObservatory repo)

In order to accurately classify sequences, we needed to look at only a single region. Classifying all of the rRNA reads would have created issues in the constant regions. We chose to only look at reads covering the V4 region because of our existing 16S data from Mendota on this region.
SortmeRNA was used to identify reads containing the V4 region because it uses a user-specified database (in this case, our V4 region iTags from the same lake) and because of its speed.
Once fasta files of V4 rRNA reads were created, they were analyzed in mothur because of the similarity between 16S workflows and what we wanted to run on these samples. This approach allowed sequence quality control and clustering of similar sequences. 
One change from our standard 16S workflow was the use of proportional data instead of rarefaction. This was chosen because of the small amount of data available and the similar sizes of all four samples.

REQUIRED SOFTWARE

SortmeRNA (v2.0)		Installation instructions: https://github.com/biocore/sortmerna/blob/master/README.md
fastx-toolkit (v0.0.12)		Installation instructions: http://hannonlab.cshl.edu/fastx_toolkit/download.html
mothur (v1.36.0)		Installation instructions: http://www.mothur.org/wiki/Installation

WORKFLOW

First step: use SortmeRNA (version 2.0). Done in Bio-Linux8 virtual machine. FW trainingset is used to make sure we didn't miss any rRNA reads with the Silva database.
The non_rRNA and _rRNA files are the output of SortmeRNA run on the quality filtered reads using the default silva databases

Input files: non_RNA fastq files from timpoints 6AM (ME15025616A), 2PM (ME15027052P), 6PM (ME15027676P), and 12AM (ME15020861012A). Freshwater database (FW_trainingset_MMBR_strict_12July12.fasta) as references (see McMahonLab repo for copy)
	indexdb_rna --ref FW_trainingset_MMBR_strict_12July12.fasta,FW.idx
	sortmerna --ref FW_trainingset_MMBR_strict_12July12.fasta,FW.idx --reads ME15025616A/ME15025616A_non_rRNA.fastq --fastx --log -v --aligned add_6AM
	sortmerna --ref FW_trainingset_MMBR_strict_12July12.fasta,FW.idx --reads ME15027052P/ME15027052P_non_rRNA.fastq --fastx --log -v --aligned add_2PM
	sortmerna --ref FW_trainingset_MMBR_strict_12July12.fasta,FW.idx --reads ME15027676P/ME15027676P_non_rRNA.fastq --fastx --log -v --aligned add_6PM
	sortmerna --ref FW_trainingset_MMBR_strict_12July12.fasta,FW.idx --reads ME1502861012A/ME1502861012A_non_rRNA.fastq --fastx --log -v --aligned add_12AM

What each flag means:
	--ref		Reference sequences to align experimental reads to
	--fastx		Write output in .fastq format
	--log		Write output to log file
	-v		Run in verbose mode
	--aligned	Name of output files


Now that I've found a few more 16S sequences, concatenate the results into a single file

	cat add_6AM.fastq ME15025616A/ME15025616A_rRNA.fastq > 6AM_rRNA.fastq
	cat add_2PM.fastq ME15027052P/ME15027052P_rRNA.fastq > 2PM_rRNA.fastq
	cat add_6PM.fastq ME15027676P/ME15027676P_rRNA.fastq > 6PM_rRNA.fastq
	cat add_12AM.fastq ME1502861012A/ME1502861012A_rRNA.fastq > 12AM_rRNA.fastq

Use some V4 tags as the reference to pull out 16S seqs aligning with the V4 region. "otus.fasta" is from the ME metagenomic sample iTAGs (available on IMG).
It contains approx. 5000 rep seqs of OTUs from Mendota and was produced by JGI's iTagger pipeline

	indexdb_rna --ref otu.fasta,otus.idx	#Creates index of the Mendota V4 reads to for SortmeRNA to align to

	sortmerna --ref otu.fasta,otus.idx --reads 6AM_rRNA.fastq --fastx --log -v --aligned 6AM_V4
	sortmerna --ref otu.fasta,otus.idx --reads 2PM_rRNA.fastq --fastx --log -v --aligned 2PM_V4 
	sortmerna  --ref otu.fasta,otus.idx --reads 6PM_rRNA.fastq --fastx --log -v --aligned 6PM_V4
	sortmerna  --ref otu.fasta,otus.idx --reads 12AM_rRNA.fastq --fastx --log -v --aligned 12AM_V4

Next, I want to classify and count the V4 regions from metatranscriptomes using mothur, but mothur will not accept .fastq files. Use the fastx-toolkit (v0.0.12) to convert the .fastq files into .fasta format. 
	
	fastq_to_fasta -v -i 6AM_V4.fastq -o 6AM_V4.fasta
	fastq_to_fasta -v -i 2PM_V4.fastq -o 2PM_V4.fasta
	fastq_to_fasta -v -i 6PM_V4.fastq -o 6PM_V4.fasta
	fastq_to_fasta -v -i 12AM_V4.fastq -o 12AM_V4.fasta

What each flag means:
	-v 	Run in verbose mode
	-i 	Input file
	-o 	Output file

I now have 4 fasta files of reads that aligned to the V4 16S region: 6AM_V4.fasta, 2PM_V4.fasta, 12AM_V4.fasta, and 6PM_V4.fasta. The Mendota V4 tags (otus.fasta) were also used as a reference for where to trim alignments. I ran the remaning commands on Zissou (McMahon Lab server, Linux environment).

	mothur
	
	make.group(fasta=6AM_V4.fasta, groups=6AM)		#Make a group file for each fasta file
	make.group(fasta=2PM_V4.fasta, groups=2PM)
	make.group(fasta=6PM_V4.fasta, groups=6PM)
	make.group(fasta=12AM_V4.fasta, groups=12AM)


	rename.seqs(fasta=6AM_V4.fasta, group=6AM_V4.groups)		#Rename so that seqs do not overlap between samples
	rename.seqs(fasta=2PM_V4.fasta, group=2PM_V4.groups)
	rename.seqs(fasta=6PM_V4.fasta, group=6PM_V4.groups)
	rename.seqs(fasta=12AM_V4.fasta, group=12AM_V4.groups)

	merge.files(input=6AM_V4.renamed.groups-2PM_V4.renamed.groups-6PM_V4.renamed.groups-12AM_V4.renamed.groups, output=dielrRNA.groups)	#Combine fasta and groups files
	merge.files(input=6AM_V4.renamed.fasta-2PM_V4.renamed.fasta-6PM_V4.renamed.fasta-12AM_V4.renamed.fasta, output=dielrRNA.fasta)


	unique.seqs(fasta=dielrRNA.fasta)		#Make unique sequences
	count.seqs(name=dielrRNA.names, group=dielrRNA.groups, processors=8)
	summary.seqs()		#803,515 unique sequences

	#Checked using system(head dielrRNA.count_table) - looks like all the samples are there!

	align.seqs(fasta=dielrRNA.unique.fasta, reference=../Classification_Databases/silva.bacteria/silva.bacteria.fasta)	#Tried using silva seed alignment, but got error. Using silva bacterial only template instead

	summary.seqs(fasta=dielrRNA.unique.align, count=dielrRNA.count_table)
	
	#Side step: run alignment on otu.fasta to find the starts and ends on this alignment database	
	align.seqs(fasta=otu.fasta, reference=../Classification_databases/silva.bacteria/silva.bacteria.fasta)
	summary.seqs(fasta=otu.align)		#Trim at positions 13862 to 23444 for V4 region when using silva.bacteria.fasta as reference for alignment

	#Filter alignment of V4 rRNA region based on results from tag data V4 alignment	
	screen.seqs(fasta=dielrRNA.unique.align, count=dielrRNA.count_table, start=13862, end=23444, maxhomop=8)
	summary.seqs(fasta=dielrRNA.unique.good.align, count=dielrRNA.good.count_table)		#Down to 867 unique seqs out of 888 seqs. Less than I'd like, but the program will run quickly
	filter.seqs(fasta=dielrRNA.unique.good.align, count=dielrRNA.good.count_table, vertical=T, trump=.)
	unique.seqs(fasta=dielrRNA.unique.good.filter.fasta, count=dielrRNA.good.count_table)		#442 unique seqs

	#Now onto the clustering and chimera removal steps. Switching to the standard lab 16S workflow for this part (see McMahonLab/16STax-Ass repo).
	pre.cluster(diffs=2, fasta=dielrRNA.unique.good.filter.unique.fasta, count=dielrRNA.unique.good.filter.count_table)
	chimera.uchime(count=dielrRNA.unique.good.filter.unique.precluster.count_table, fasta=dielrRNA.unique.good.filter.unique.precluster.fasta)	#22 chimeras found
	remove.seqs(fasta=dielrRNA.unique.good.filter.unique.precluster.fasta, accnos=dielrRNA.unique.good.filter.unique.precluster.denovo.uchime.accnos, count=dielrRNA.unique.good.filter.unique.precluster.count_table)

	#Do a classification to identify and remove chloroplasts
	classify.seqs(fasta=dielrRNA.unique.good.filter.unique.precluster.pick.fasta, template=../Classification_Databases/gg_13_5_97_otus_noFW.fasta, taxonomy=../Classification_Databases/gg_13_5_97_otus_noFW.taxonomy, cutoff=60)
	remove.lineage(taxon=k__Bacteria;p__Cyanobacteria;c__Chloroplast;, count=dielrRNA.unique.good.filter.unique.precluster.pick.count_table, fasta=dielrRNA.unique.good.filter.unique.precluster.pick.fasta)

	system(mv dielrRNA.unique.good.filter.unique.precluster.pick.pick.fasta qc.dielrRNA.fasta)		#Shorten up those file names
	system(mv dielrRNA.unique.good.filter.unique.precluster.pick.pick.count_table qc.dielrRNA.count_table)

	#On to OTU clustering. This is a pretty small dataset, so I'll use dist.seqs() + cluster() (de novo clustering into OTUs by distance)
	dist.seqs(fasta=qc.dielrRNA.fasta, cutoff=0.020)
	cluster(column=qc.dielrRNA.dist, count=qc.dielrRNA.count_table)
	make.shared(list=qc.dielrRNA.an.unique_list.list, count=qc.dielrRNA.count_table, label=0.02)
	get.oturep(column=qc.dielrRNA.dist, count=qc.dielrRNA.count_table, list=qc.dielrRNA.an.unique_list.list, fasta=qc.dielrRNA.fasta, label=0.02)

	#Note: I'll use the proportional method instead of the rarefaction method this time. See SUMMARY for why. This will be implemented in downstream R analysis

	#mothur output is looking great. Now to classify using Robin's method (McMahonLab/16Tax-Ass repo)
	#First step is to remove white space in fasta file headers. replacing them with a colon. White space confuses BLAST.
	tr '\t' ':' <qc.dielrRNA.an.unique_list.0.02.rep.fasta> qc.dielrRNA.repseqs.nowhitespace.fasta

	#Make blast database of freshwater database	
	makeblastdb -dbtype nucl -in ../Classification_Databases/FW_trainingset_MMBR_strict_12July12.fasta -input_type fasta -parse_seqids -out FW.db

	#blast the rRNA representative sequences against the freshwater database
	blastn -query qc.dielrRNA.repseqs.nowhitespace.fasta -task megablast -db FW.db -out otus.custom.blast -outfmt 11 -max_target_seqs 1
	blast_formatter -archive dielrRNA.blast -outfmt "6 qseqid pident length qlen qstart qend" -out dielrRNA.blast.table

	#Sort sequences into good hits and not good hits, then format for classification in mothur
	Rscript ../16TaxAss-master/tax-scripts/find_seqIDs_with_pident.R dielrRNA.blast.table ids.above.98 stats1.txt 98 TRUE
	Rscript ../16TaxAss-master/tax-scripts/find_seqIDs_with_pident.R dielrRNA.blast.table ids.below.98 stats2.txt 98 FALSE
	python ../16TaxAss-master/tax-scripts/fetch_seqIDs_blast_removed.py qc.dielrRNA.repseqs.nowhitespace.fasta dielrRNA.blast.table ids.missing
	cat ids.below.98 ids.missing > ids.below.98.all
	python ../16TaxAss-master/tax-scripts/fetch_fastas_with_seqIDs.py ids.above.98 qc.dielrRNA.repseqs.nowhitespace.fasta otus.above.98.fasta
	python ../16TaxAss-master/tax-scripts/fetch_fastas_with_seqIDs.py ids.below.98.all qc.dielrRNA.repseqs.nowhitespace.fasta otus.below.98.fasta

	#Classify sequences that hit the FW database with the FW database (custom), and sequences that did not hit with the Greengenes database (general)
	mothur
	classify.seqs(fasta=otus.above.98.fasta, template=custom.fasta,  taxonomy=custom.taxonomy, method=wang, probs=T, processors=2)
	classify.seqs(fasta=otus.below.98.fasta, template=general.fasta, taxonomy=general.taxonomy, method=wang, probs=T, processors=2)
	quit()
	
	#Combine into a single taxonomy file
	cat otus.above.98.custom.wang.taxonomy otus.below.98.general.wang.taxonomy > otus.taxonomy

Finished! The output files are:

qc.dielrRNA.an.unique_list.0.02.rep.fasta	Representative sequences for each OTU
qc.dielrRNA.an.unique_list.shared		OTU table
qc.dielrRNA.taxonomy				Classification of each OTU.

Further processing took place in R - see script dielrRNA_intro.Rmd or .html



DIRECTORY STRUCTURE
OMD-TOIL
|- rRNA_analysis
||- README.txt						This document
||- dielrRNA_workflow.txt				Commands used to perform this analysis
||- rRNA_V4region_files					fasta files of reads covering the V4 region
|||- 6AM_V4.fasta
|||- 6PM_V4.fasta
|||- 2PM_V4.fasta
|||- 12AM_V4.fasta
||- rRNA_mothur_results					Output of mothur analysis and classification
|||- qc.dielrRNA.taxonomy				Taxonomic classifications
|||- qc.dielrRNA.an.unique_list.0.02.rep.fasta		Representative sequences
|||- qc.dielrRNA.an.unique_list.shared			OTU table
||- rRNA_downstream_analysis
|||- dielrRNA_intro.html				HTML R Markdown document; contains overview of taxa present in the dataset
|||- dielrRNA_intro.Rmd					Code to produce R Markdown document
|||- dielrRNA_expanded_taxonomy.csv			Taxonomy file with assignments split into multiple columns, each representing a phylogenetic level
|||- dielrRNA_parsed_taxonomy.csv			Taxonomy file with shortened OTU names
|||- dielrRNA_otutable.csv				OTU table derived from .shared file that is prepped for analysis with R package OTUtable

Not shown: original data files, intermediate files from SortmeRNA.
These files are too large to put on Github. Please contact trina.mcmahon@wisc.edu or amlinz@wisc.edu if you would like these files.

