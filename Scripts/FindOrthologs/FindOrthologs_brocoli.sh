#!/bin/bash

module add Phylogeny/ete3/3.1.1
module add intel/intelpython3

# Scripts
broccoli="Scripts/FindOrthologs/Broccoli-master/broccoli.py"
fasttree="/software/Phylogeny/FastTree/2.1.10/bin/FastTree"
diamond="Scripts/FindOrthologs/diamond"

# Files & parameters
ResultsFolder="Results/FindOrthologs"

################################################
### Run BROCCOLI
gunzip ${ResultsFolder}/CheckedProteomes/*.fa.gz
mkdir -p ${ResultsFolder}/broccoli
cd ${ResultsFolder}/broccoli
python ../../../${broccoli} -dir ../CheckedProteomes -ext fa -path_diamond ../../../${diamond} -path_fasttree ${fasttree}
cd ../../..
gzip ${ResultsFolder}/CheckedProteomes/*.fa






