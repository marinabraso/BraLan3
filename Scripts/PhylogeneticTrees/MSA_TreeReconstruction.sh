#!/bin/bash

# Scripts
module add mafft/7.475
pwd=$(pwd)
RAxML=${pwd}/../Software/standard-RAxML-master/raxmlHPC

# Files & parameters
RunName=$1
Start=$2
Step=$3
OrthologsFolder="Results/FindOrthologs/${RunName}_broccoli"
ResultsFolder="Results/PhylogeneticTrees"
MSAFolder=${ResultsFolder}/MSA_mafft
TreeFolder=${ResultsFolder}/RAxML_Trees
SeqFolder=${ResultsFolder}/OGSequences_${RunName}
mkdir -p ${MSAFolder} ${TreeFolder}

for group in $(cut -f1 ${OrthologsFolder}/dir_step3/table_OGs_protein_counts.txt | tail -n +2 | tail -n +${Start} | head -${Step})
do
	echo ${group}
	## Alignment of AA sequences with MAFFT
	if [[ ! -s ${MSAFolder}/${group}_AA.fa ]]
	then
		echo "mafft for ${group}"
		mafft --globalpair --quiet --anysymbol --allowshift ${SeqFolder}/${group}_AA.fa | awk '{if($1 ~ />/){print seq; print $0; seq=""}else{seq=seq$0}}END{print seq}' | tail -n +2 > ${MSAFolder}/${group}_AA.fa
	fi

	## RAxML tree with PROTGAMMAAUTO
	if [[ -s ${MSAFolder}/${group}_AA.fa ]] && [[ ! -s ${TreeFolder}/RAxML_result.${group}_PROTGAMMAAUTO ]]
	then
		echo "RAxML for ${group}" 
		cd ${TreeFolder}
		${RAxML} -s ${pwd}/${MSAFolder}/${group}_AA.fa -n ${group}_PROTGAMMAAUTO -m PROTGAMMAAUTO -p 12345
		cd ${pwd}
	fi
done














