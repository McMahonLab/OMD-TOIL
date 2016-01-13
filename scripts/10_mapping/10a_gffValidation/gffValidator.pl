###############################################################################
# gffValidator.pl
# Copyright (c) 2015, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# This function identifies errors in the GFF files which will prevent them from
# being properly parsed by ht-seq. It relies on the Genome Tools software
# package, available at http://genometools.org/pub/.
################################################################################
#!/usr/bin/perl
use strict;
use warnings;

################################################################################
### User-defined files and folder structure
################################################################################
my $gffFolder = '../../../rawData/refGenomes/gff';

################################################################################
### Initialize lists of refGenomes
################################################################################
my @genomeList = glob($gffFolder.'/*.gff');
for (@genomeList) {
  s/.+\/(.+).gff/$1/
}

################################################################################
### Validate the gff Files
################################################################################
foreach my $genome (@genomeList) {
  system('gt gff3validator '.$gffFolder.'/'.$genome.'.gff &> '.$genome.'.out');
}

################################################################################
### Loop over the output .out files and write errors to all.out
################################################################################

# Create output file
open(my $file, '>', $gffFolder.'/all.out');

# Stream individual gff files to the all.gff file
foreach my $genome (@genomeList) {
  my $gff = $genome.".out";
  open(my $gffFile, '<', $gff);
  foreach my $line (<$gffFile>) {
    if ($line =~ m/.+error.+/) { print $file $line; };
  }
  close $gffFile;
}
close $file;

system('rm *.out');
