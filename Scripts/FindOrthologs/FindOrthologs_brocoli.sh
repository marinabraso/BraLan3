#!/bin/bash

# Scripts
broccoli="Scripts/FindOrthologs/Broccoli-master/broccoli.py"
fasttree="/software/Phylogeny/FastTree/2.1.10/bin/FastTree"
diamond="Scripts/FindOrthologs/diamond"

# Files & parameters
ResultsFolder="Results/FindOrthologs"

################################################
### Create and activate conda environment
#conda create -n env-broccoli python=3.6 ete3
#conda init bash
#conda activate env-broccoli

################################################
### Run BROCCOLI
gunzip ${ResultsFolder}/CheckedProteomes/*.fa.gz
mkdir -p ${ResultsFolder}/broccoli
cd ${ResultsFolder}/broccoli
python ../../../${broccoli} -dir ../CheckedProteomes -ext fa -path_diamond ../../../${diamond} -path_fasttree ${fasttree}
cd ../../..
gzip ${ResultsFolder}/CheckedProteomes/*.fa





