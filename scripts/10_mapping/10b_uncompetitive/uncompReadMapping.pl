###############################################################################
# uncompReadMapping.pl
# Copyright (c) 2015, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# This is a wrapper funtion to map metagenomics reads to reference genomes
# using BWA. It was developed to map reads from OMD-TOIL to our SAGs and MAGs.
################################################################################
#!/usr/bin/perl
use strict;
use warnings;

################################################################################
### User-defined files and folder structure
################################################################################

my $mtFolder = '../../../rawData/';
my $genomeFolder = '../../../rawData/refGenomes/fna';
my $mapFolder = '../../../derivedData/mapping/uncompetitive/bamFiles';
my $numProcs = 30;

my @genomeList = glob($genomeFolder.'/*.fna');
my @sampleList = glob($mtFolder.'/*');

################################################################################
### Step 1: Index Genomes
################################################################################

my $int = 1;
foreach my $genomePath (@genomeList) {
  print "Indexing genome ".$int." of ".@genomeList."\n";
  if ($genomePath =~ /(.+).fna/) {
    my $genome = $1;
    system("bwa index -p $genome -a is $genomePath");
    $int++;
    }
}

################################################################################
### Step 2: Map Metatranscriptomes to Reference Genomes
################################################################################

$int = 1;
my $total = @genomeList*@sampleList;

foreach my $samplePath (@sampleList) {
  foreach my $genomePath (@genomeList) {
    if ($genomePath =~ /.+\/(.+).fna/) {
      my $genome = $1;
      if ($samplePath =~ /.+\/(.+)/) {
	my $sample = $1;
	print "Mapping sample ".$sample." against genome ".$genome." (".$int." of ".$total."). \n";
	$int ++;
	system("bwa mem -t ".$numProcs." ".$genomeFolder."/".$genome." ".$samplePath."/rRNA_removal/".$sample."_non_rRNA.fastq > ".$mapFolder."/".$sample."-".$genome.".sam");
      }
    }
  }
}

################################################################################
### Step 3: Manipulate Output SAM files. Convert to BAM, sort and index. Delete
### unsorted SAM and BAM files to save space.
################################################################################

my @samList = glob($mapFolder.'/*.sam');
$int = 1;

foreach my $samPath (@samList) {
  print "Processing SAM file ".$int." of ".@samList."\n";
  if ($samPath =~ /.+\/(.+).sam/) {
    my $sam = $1;
    system("samtools view -b -S -o ".$mapFolder."/".$sam.".bam ".$mapFolder."/".$sam.".sam");
    system("samtools sort -o ".$mapFolder."/".$sam.".sorted.bam -O bam -T ".$mapFolder."/temp/".$sam." ".$mapFolder."/".$sam.".bam");
    system("samtools index ".$mapFolder."/".$sam.".sorted.bam");
    system("rm ".$mapFolder."/".$sam.".bam");
    system("rm ".$mapFolder."/".$sam.".sam");
    $int ++;
  }
}
