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
use Parallel::ForkManager;

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
### Generate a list of commands to run
################################################################################

my @commandArray;

foreach my $samplePath (@sampleList) {
  foreach my $genomePath (@genomeList) {
    push @commandArray, '/usr/local/bin/htseq-count -f sam -r pos -s no -a 0 -t tRNA -i locus_tag -m intersection-strict -o '. $countFolder.'/'.$mt.'-'.$genome.'.tRNA.sam '.$mapFolder.'/'.$mt.'-'.$genome.'.sam '.$gffFolder.'/'.$genome.'.gff >  '.$countFolder.'/'.$mt.'-'.$genome.'.tRNA.out';
  }
}
################################################################################
### Loop to obtain count data
################################################################################

my $numChildren = 25;
my $pm = new Parallel::ForkManager($numChildren);

foreach my $command (@commandArray) {
  # Forks and returns the pid for the child:
  my $pid = $pm->start and next;

  system($command);
  
  $pm->finish; # Terminates the child process
}

$pm->wait_all_children;
print "Counting complete!\n";
