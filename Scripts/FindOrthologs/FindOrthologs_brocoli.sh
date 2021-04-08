#!/bin/bash

module add Phylogeny/ete3/3.1.1
module add intel/intelpython3

# Scripts
broccoli="Scripts/FindOrthologs/Broccoli-master/broccoli.py"
fasttree="/software/Phylogeny/FastTree/2.1.10/bin/FastTree"
diamond="Scripts/FindOrthologs/diamond"

# Files & parameters
ResultsFolder="Results/FindOrthologs"
ProteomesFolder="Results/FilteringGeneSets/Proteomes"

################################################
### Run BROCCOLI
gunzip ${ProteomesFolder}/*.fa.gz 2> ~/null
mkdir -p ${ResultsFolder}/broccoli
cd ${ResultsFolder}/broccoli
python ../../../${broccoli} -dir ../../../${ProteomesFolder} -ext fa -path_diamond ../../../${diamond} -path_fasttree ${fasttree}
cd ../../..
gzip ${ProteomesFolder}/*.fa 2> ~/null






