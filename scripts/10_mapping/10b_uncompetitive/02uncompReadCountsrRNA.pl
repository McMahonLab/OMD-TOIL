###############################################################################
# readCounts.py
# Copyright (c) 2015, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# This is a wrapper funtion to count the MT reads which map to each gene in a
# genome. Run this after you have run MTwrapperFunction.pl
################################################################################
#!/usr/bin/perl
use strict;
use warnings;

# Import packages
use File::Path qw(make_path);

################################################################################
### User-defined files and folder structure
################################################################################
my $genomeFolder = '../../../data/refGenomes/fna';
my $gffFolder = '../../../data/refGenomes/gff';
my $mtFolder = '../../../data/rawData/';
my $mapFolder = '../../../data/derivedData/mapping/uncompetitive/bamFiles';
my $countFolder = '../../../data/derivedData/mapping/uncompetitive/readCounts';

# Check if the output directory exists and create it if necessary
if (-d $countFolder) {
  print "Output directory exists \n";
}
# If not, create it
else {
  print "Creating output directory \n";
  make_path($countFolder);
}

################################################################################
### Initialize lists of refGenomes and MTs
################################################################################

my @genomeList = glob($genomeFolder.'/*.fna');
@genomeList = grep !/merged.fna/, @genomeList;
for (@genomeList) {
  s/.+\/(.+).fna/$1/
  }

my @mtList;
opendir (DIR, $mtFolder);
while (my $dir = readdir(DIR)) {
  next if ($dir =~ m/^\./);
  push @mtList, $dir;
}
closedir(DIR);

################################################################################
### Loop to obtain count data
################################################################################

# Loop over each genome and MT
foreach my $genome (@genomeList) {
    print "Processing genome ".$genome."\n";
    foreach my $mt (@mtList) {
      print "Processing metranscriptome ".$mt."\n";
      # Call HTseq-count
      system('/usr/local/bin/htseq-count -f sam -r pos -s no -a 0 -t rRNA -i locus_tag -m intersection-strict -o '. $countFolder.'/'.$mt.'-'.$genome.'.rRNA.sam '.$mapFolder.'/'.$mt.'-'.$genome.'.sam '.$gffFolder.'/'.$genome.'.gff >  '.$countFolder.'/'.$mt.'-'.$genome.'.rRNA.out');
      }
  }
