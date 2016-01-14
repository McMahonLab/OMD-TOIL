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
In order for to properly count mapped reads, the GFF files must be in the proper format. The script `gffValidator` will analyze all GFF files in the `gffFolder` folder, and generate a file `gffFolder/all.out` describing errors in the GFF files. The GFF files must then be manually corrected. In my experience, errors will be one of the following messages:

* `token "string" on line XYZ in file “genome.gff" does not contain exactly one ‘=‘`

  This error occurs when a gene annotation contains a semi-colon. According to the GFFv3 [standard](http://www.sequenceontology.org/gff3.shtml), semi-colons have a reserved meaning and must be escaped. The standard recommends using the URL escape `%3B`.

* `error: could not parse score '' on line XYZ in file ‘genome.gff’`

  This error occurs because a line of the GFF file is incomplete (e.g., missing the score information). I found this to occur only for CRISPR arrays defined in the final line of the GFF file. Since I am not concerned with mapping to CRISPR arrays, I deleted those lines.
