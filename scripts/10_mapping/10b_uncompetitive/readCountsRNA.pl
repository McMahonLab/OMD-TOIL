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

################################################################################
### User-defined files and folder structure
################################################################################
my $genomeFolder = '../../../rawData/refGenomes/fna';
my $gffFolder = '../../../rawData/refGenomes/gff';
my $mtFolder = '../../../rawData/';
my $mapFolder = '../../../derivedData/mapping/uncompetitive/bamFiles';
my $countFolder = '../../../derivedData/mapping/uncompetitive/readCounts';

mkdir $countFolder unless -d $countFolder;

################################################################################
### Initialize lists of refGenomes and MTs
################################################################################

my @genomeList = glob($genomeFolder.'/*.fna');
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
      system('/usr/local/bin/htseq-count -f bam -r pos -s no -a 0 -t RNA -i locus_tag -m intersection-strict -o '. $countFolder.'/'.$mt.'-'.$genome.'.RNA.sam '.$mapFolder.'/'.$mt.'-'.$genome.'.sorted.bam '.$gffFolder.'/'.$genome.'.gff >  '.$countFolder.'/'.$mt.'-'.$genome.'.RNA.out');
      }
  }
