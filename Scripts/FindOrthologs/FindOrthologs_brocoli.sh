#!/bin/bash

pwd=$(pwd)
conda activate broccoli

# Scripts
broccoli="../Software/Broccoli-master/broccoli.py"
diamond="../Software/diamond"
FastTree="../Software/FastTree"

# Files & parameters
RunName=$1
ResultsFolder="Results/FindOrthologs"
ProteomesFolder="Results/FilteringGeneSets/Proteomes${RunName}"

################################################
mkdir -p ${ResultsFolder}/${RunName}_broccoli
cd ${ResultsFolder}/${RunName}_broccoli
python ${pwd}/${broccoli} -dir ${pwd}/${ProteomesFolder} -ext fa -path_diamond ${pwd}/${diamond} -path_fasttree ${pwd}/${FastTree}
cd ${pwd}





