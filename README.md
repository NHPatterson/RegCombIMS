# RegCombIMS
Code accompanying ["Co-registration and analysis of multiple imaging mass spectrometry datasets targeting different analytes"](https://doi.org/10.1093/bioinformatics/bty780) Please cite this work if you use this tool in your own analysis.

This repository contains code to demo the analysis pipeline described in "Co-registration and analysis of multiple imaging mass spectrometry datasets targeting different analytes"

RegCombIMS extends Cardinal's excellent data structures for IMS for computational registration tasks using the RNiftyReg library.

## Installation
This package can be installed from github in R using the devtools package and following code:

```R
#install Cardinal following instructions at https://cardinalmsi.org/
#install RNiftyReg
install.packages("RNiftyReg")
install.packages("devtools") # if you have not installed "devtools" package
devtools::install_github("nhpatterson/RegCombIMS")
```

## Demo
[See demo example combining two datasets](https://github.com/NHPatterson/RegCombIMS/blob/master/scripts/DAN_agLDI_combination_script.R)

### Contact:
Heath Patterson nathan.h.patterson@vanderbilt.edu

Please use the GitHub issue tracker for any errors encountered using this tool. I will try to solve them ASAP.
