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

# Import packages
use File::Path qw(make_path);

################################################################################
### User-defined files and folder structure
################################################################################

my $mtFolder = '../../../data/rawData/';
my $genomeFolder = '../../../data/refGenomes/fna';
my $mapFolder = '../../../data/derivedData/mapping/uncompetitive/bamFiles';
my $numProcs = 30;

# Check if the output directory exists and create it if necessary
if (-d $mapFolder) {
  print "Output directory exists \n";
}
# If not, create it
else {
  print "Creating output directory \n";
  make_path($mapFolder);
}

my @genomeList = glob($genomeFolder.'/*.fna');
@genomeList = grep !/merged.fna/, @genomeList;
my @sampleList = glob($mtFolder.'/*');

################################################################################
### Step 2: Simultaneously Index and Map Metatranscriptomes to Reference Genomes
################################################################################

my $int = 1;
my $total = @genomeList*@sampleList;

foreach my $samplePath (@sampleList) {
  foreach my $genomePath (@genomeList) {
    if ($genomePath =~ /.+\/(.+).fna/) {
      my $genome = $1;
      if ($samplePath =~ /.+\/(.+)/) {
	my $sample = $1;
	print "Mapping sample ".$sample." against genome ".$genome." (".$int." of ".$total."). \n";
	$int ++;
	system("bbmap.sh t=".$numProcs." in=".$samplePath."/rRNA_removal/".$sample."_non_rRNA.fastq outm=".$mapFolder."/".$sample."-".$genome.".sam ref=".$genomeFolder."/".$genome.".fna nodisk sam=1.3");
      }
    }
  }
}
