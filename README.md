OMD-TOIL
===
Copyright (c) 2016, McMahon Lab  
URL: [https://mcmahonlab.wisc.edu/](https://mcmahonlab.wisc.edu/)  
URL: [https://github.com/McMahonLab/](https://github.com/McMahonLab/)  
All rights reserved.

Overview
--
During August 20-21, 2015, the McMahon Lab carried out a 24 hour sampling of Lake Mendota. This repo describes the data and analysis associated with this expedition, known as OMD-TOIL (Operation Mendota Drain: Transcriptomics of Inland Lakes).

Samples were collected approximately every two hours from 6am on August 20 to 6am on August 21. A subset of samples were subject to RNA extraction, rRNA depletion, spiking with an internal standard, and sent to the [UW Biotech Center](https://www.biotech.wisc.edu/) for sequencing. Resulting sequences were subject to computational rRNA removal and downstream analysis.

This repository ('repo') describes the experimental and analytical procedures associated with this project.

Data Organization
--
We have established separate locations for storage of data and analysis associated with this project. __Data__ are stored on the McMahon Lab's private server (Zissou) in the `data_lakes/OMD-TOIL` folder. __Protocols and scripts__ are stored in the [OMD-TOIL](https://github.com/McMahonLab/OMD-TOILv2) repo on the McMahon Lab's GitHub account. A clone of this repo exists on Zissou, in `/home/shared/OMD-TOIL`. This folder also contains symbolic links to the data. That way, your scripts can access the data folder using the appropriate relative path within `/home/shared/OMD-TOIL`.

Repo Table of Contents
--
This repo describes the following analyses. Clicking on each link will take you to the appropriate location.

1. Sample collection
2. Internal standard addition
3. Sample RNA extraction
4. cDNA synthesis & rRNA removal
5. Library preparation
6. DNA sequencing
7. Quality-filtering & trimming
8. rRNA removal in silico
9. rRNA analysis
10. [Mapping](scripts/10_mapping/README.md)  
  10a. [Validation of GFF files](scripts/10_mapping/10a_gffValidation/README.md)  
  10b. [Uncompetitive Mapping](scripts/10_mapping/10b_uncompetitive/README.md)  
  10c. [Competitive Mapping](scripts/10_mapping/10c_competitive/README.md)  

Protocols associated with the __experimental__ steps (1 to 6) are in the __protocols__ folder. Scripts and workflows associated with the __computational__ steps (7 forward) are in the __scripts__ folder. Metadata about the samples are in the __metadata__ folder. Each step is contained within a numbered and named folder, such as `protocols/01_sample_collection`.

How To Contribute
--
Lab members who wish to contribute to the project should contact [Josh Hamilton](https://github.com/joshamilton).

To contribute to this repo, create a fork in your personal GitHub account and clone the repo to your personal computer (or your `home` folder on Zissou.) Changes to the repo should be made in a new branch. When you have made your changes, submit a pull request, requesting your new branch be merged into `McMahonLab:master`.

When submitting a pull request, please adhere to the following guidelines:

* __Follow the proper file structure__, as described in the Table of Contents section.
* __All files should be text-based__. Git cannot track changes made to binary file formats such as doc, xls, pdf, etc.
* __Each set of scripts should be accompanied by a README file__, which contains the information described below. We recommend using [MarkDown](https://help.github.com/articles/markdown-basics/), as it renders nicely on GitHub.
* __Data are not to be stored in this repo__. Instead, upload your data to `/home/shared/OMD-TOIL/data` on Zissou, and indicate these data in your pull request. The data will be moved to `data_lakes` and symbolic links created. Additionally, you may wish to retain your data files in your local branch of the repo; in that case, you may indicate these data files in a [.gitignore](https://help.github.com/articles/ignoring-files/) file committed in your local repository.
* __Update the master README.md file__ with a new Table of Contents entry and link.

Contents of a Good README File
--
If you are contributing scripts, workflows, or analysis to this repo, please provide a README which contains the following information:

* __A header__, indicating who wrote the README and did the analysis, along with your contact information
* __A summary of the prodecure__, explaining the purpose of this step and a justification of the approach
* __A list of required software__, with instructions on how to install it
* __Scripts or workflows__ that run the analysis (a script) or walk the user through the necessary commands (a workflow).
* __Documentation__ for each script, including the expected inputs and outputs

The goal is to make your analysis fully reproducible. The documentation should enable another scientist to fully recreate your analysis and results! Feel free to check out other McMahon lab repos for README ideas.

Repo Structure
--
    project
    |- README                                             # Overview of the repo and its contents
    ||- scripts						  
    ||- metadata					  # Metadata collected at each sampling timepoint
