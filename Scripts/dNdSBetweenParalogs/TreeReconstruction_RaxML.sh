#!/bin/bash

# Scripts
module add Phylogeny/raxml/8.2.12

# Files & parameters
OrthologsFolder="Results/FindOrthologs"
ResultsFolder="Results/dNdSBetweenParalogs"
mkdir -p ${ResultsFolder}/TreeReconstruction
pwd=$(pwd)

for group in $(cut -f1 ${OrthologsFolder}/broccoli/dir_step3/table_OGs_protein_counts.txt | tail -n +2 | head -2)
do
	echo ${group}
	#if [[ ! -s ${ResultsFolder}/TreeReconstruction/RAxML_result.${group}_GTRGAMMA ]]
	#then
		cd ${ResultsFolder}/TreeReconstruction
		raxmlHPC -s ../MSA_mafft/${group}_DNA.fa -n ${group}_GTRGAMMA -m GTRGAMMA -p 12345
		cd $pwd
	#fi
done














