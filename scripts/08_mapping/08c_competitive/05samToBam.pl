###############################################################################
# compReadMapping.pl
# Copyright (c) 2015, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# This is a wrapper function for competitive mapping of metagenomic reads to
# reference genomes using BWA. The function also counts mapped reads. It was
# developed for the OMD-TOIL project.
################################################################################
#!/usr/bin/perl
use strict;
use warnings;

# Import packages
use File::Path qw(make_path);

################################################################################
### User-defined files and folder structure
################################################################################
my $mtFolder = '../../../data/rawData/';
my $genomeFolder = '../../../data/refGenomes/fna';
my $gffFolder = '../../../data/refGenomes/gff';
my $mapFolder = '../../../data/derivedData/mapping/competitive/bamFiles';
my $countFolder = '../../../data/derivedData/mapping/competitive/readCounts';
my $numProcs = 30;

# Check if the output directories exists and create them if necessary
if (-d $mapFolder) {
  print "Output directory exists \n";
}
# If not, create it
else {
  print "Creating output directory \n";
  make_path($mapFolder);
}

if (-d $countFolder) {
  print "Output directory exists \n";
}
# If not, create it
else {
  print "Creating output directory \n";
  make_path($countFolder);
}

my @genomeList = glob($genomeFolder."/*.fna");
my @sampleList = glob($mtFolder."/*");
my @featureList = ("CDS", "rRNA", "tRNA", "RNA");

################################################################################
### Step 0: Merge FASTA and GFF files
################################################################################
# Remove these files if already present
#system("rm ".$gffFolder."/merged.gff");
#system("rm ".$genomeFolder."/merged.fna");
#system("cat ".$gffFolder."/*.gff > ".$gffFolder."/merged.out");
#system("mv ".$gffFolder."/merged.out ".$gffFolder."/merged.gff");
#system("cat ".$genomeFolder."/*.fna > ".$genomeFolder."/merged.out");
#system("mv ".$genomeFolder."/merged.out ".$genomeFolder."/merged.fna");

################################################################################
### Step 1: Index Genomes
################################################################################

#system("bwa index -p ".$genomeFolder."/merged -a is ".$genomeFolder."/merged.fna");

################################################################################
### Step 2: Map Metatranscriptomes to Reference Genomes
################################################################################

my $int = 1;
#my $total = @sampleList;

#foreach my $samplePath (@sampleList) {
#  if ($samplePath =~ /.+\/(.+)/) {
#    my $sample = $1;
#    if ($sample eq "ME15025616A") {
#      print "Mapping sample ".$sample." (".$int." of ".$total."). \n";
#      $int ++;
#      system("bwa mem -t ".$numProcs." ".$genomeFolder."/merged  ".$samplePath."/".$sample."_non_rRNA.fastq > ".$mapFolder."/".$sample.".sam");
#    }
#  }
#}

################################################################################
### Step 3: Manipulate Output SAM files. Convert to BAM, sort and index. Delete
### unsorted SAM and BAM files to save space.
################################################################################

my @samList = glob($mapFolder."/*.sam");
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

################################################################################
### Step 4: Obtain count data.
################################################################################

#$int = 1;
#foreach my $feature (@featureList) {
#  foreach my $samPath (@samList) {
#    print "Counting reads for MT ".$int." of ".@samList."\n";
#    if ($samPath =~ /.+\/(.+).sam/) {
#      my $sam = $1;
#      system("/usr/local/bin/htseq-count -f bam -r pos -s no -a 0 -t ".$feature." -i locus_tag -m intersection-strict -o ". $countFolder."/".$sam.".".$feature.".sam ".$mapFolder."/".$sam.".sorted.bam ".$gffFolder."/merged.gff >  ".$countFolder."/".$sam.".".$feature.".out");
#    }
#  }
#  $int ++;
#}
