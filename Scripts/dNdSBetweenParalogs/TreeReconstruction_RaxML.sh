#!/bin/bash

# Scripts
module add Phylogeny/raxml/8.2.12

# Files & parameters
OrthologsFolder="Results/FindOrthologs"
ResultsFolder="Results/dNdSBetweenParalogs"
mkdir -p ${ResultsFolder}/TreeReconstruction
pwd=$(pwd)

cd ${ResultsFolder}/TreeReconstruction
for group in $(cat ../Groups_wChodates.txt)
do
	echo ${group}
	if [[ ! -s RAxML_result.${group}_GTRGAMMA ]]
	then
		raxmlHPC -s ../MSA_mafft/${group}_DNA.fa -n ${group}_GTRGAMMA -m GTRGAMMA -p 12345
	fi
done
cd $pwd














