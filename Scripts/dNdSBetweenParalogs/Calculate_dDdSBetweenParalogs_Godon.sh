#!/bin/bash

# Scripts
pwd=$(pwd)
godon=${pwd}/../Software/godon
RAxML=${pwd}/../Software/standard-RAxML-master/raxmlHPC

# Files & parameters
RunName=$1
OrthologsFolder="Results/FindOrthologs/${RunName}_broccoli"
ResultsFolder="Results/dNdSBetweenParalogs"
mkdir -p ${ResultsFolder}/Godon_M8 ${ResultsFolder}/TreeReconstruction


################################################
### RAxML & godon M8
for group in $(cut -f1 ${OrthologsFolder}/dir_step3/table_OGs_protein_counts.txt | tail -n +2 | head -500)
do
	echo ${group}
	missingOG=$(grep -w ${group} ${ResultsFolder}/MissingOG_TooFewGenes.txt | wc -l)
	if [[ ! -s ${ResultsFolder}/Godon_M8/${group}_M8_likelihood_ratio.txt ]] && [[ ${missingOG} -eq 0 ]]
	then
		cd ${ResultsFolder}/TreeReconstruction
		${RAxML} -s ${pwd}/${ResultsFolder}/MSA_cleaning_TCoffee/${group}_DNA_clean.fa -n ${group}_GTRGAMMA -m GTRGAMMA -p 12345
		cd ${pwd}
		if [[ -s ${ResultsFolder}/TreeReconstruction/RAxML_result.${group}_GTRGAMMA ]]
		then
			${godon} test M8 --m0-tree ${ResultsFolder}/MSA_cleaning_TCoffee/${group}_DNA_clean.fa ${ResultsFolder}/TreeReconstruction/RAxML_result.${group}_GTRGAMMA > ${ResultsFolder}/Godon_M8/${group}_M8_likelihood_ratio.txt
		else
			echo ${group} >> ${ResultsFolder}/MissingOG_TooFewGenes.txt
		fi

	fi
done





