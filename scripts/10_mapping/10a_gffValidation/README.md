OMD-TOIL: Mapping Metatranscriptomic Reads: Validation of GFF files
===
Copyright (c) 2016, McMahon Lab  
URL: [https://mcmahonlab.wisc.edu/](https://mcmahonlab.wisc.edu/)  
URL: [https://github.com/McMahonLab/](https://github.com/McMahonLab/)  
All rights reserved.

Documentation and analysis by Joshua J. Hamilton  
URL: [https://github.com/joshamilton/](https://github.com/joshamilton/)  

Validation of GFF files
--
In order for to properly count mapped reads, the GFF files must be in the proper format. The script `gffValidator` will analyze all GFF files in the `gffFolder` folder, and generate a file `gffFolder/all.out` describing errors in the GFF files. The GFF files must then be manually corrected. In our experience, errors will be one of the following messages:

* `token "string" on line XYZ in file “genome.gff" does not contain exactly one ‘=‘`

  This error occurs when a gene annotation contains a semi-colon. According to the GFFv3 [standard](http://www.sequenceontology.org/gff3.shtml), semi-colons have a reserved meaning and must be escaped. The standard recommends using the URL escape `%3B`.

* `error: could not parse score '' on line XYZ in file ‘genome.gff’`

  This error occurs because a line of the GFF file is incomplete (e.g., missing the score information). This error only seems to occur for CRISPR arrays defined in the final line of the GFF file. Unless you are concerned with mapping to the arrays, we suggest deleting those line.

* `error: 'unmatched quote'`

  This error occurs because an annotation contains a (single) double quote. For example, `ADP-ribose 1"-phosphate phophatase related protein`. This error can be resolved by replacing the double quote with a single quote: `ADP-ribose 1'-phosphate phophatase related protein`.

Validation of Fasta header files
--
In order for to properly count mapped reads, the `seqID` field in the GFF files must match the description line in the fasta files. The script `fnaValidator` updates the fasta description to use the short form used in the `seqID` field.

I don't know why the gff and fasta files from IMG have different values for this field, but they do. Briefly, the fasta files contain a description of the form:

`ME10739DRAFT_MEint_metabat_10739_757002236.1 Composite genome from Lake Mendota Epilimnion pan-assembly MEint.metabat.10739 : ME10739DRAFT_MEint_metabat_10739_757002236.1`

while the GFF files have a seqID of the form:
`ME10739DRAFT_MEint_metabat_10739_757002236.1`

The desired fasta description line is limited to the material after the semi-colon. The script `fnaValidator` rewrites the fasta files with the correct description line.
