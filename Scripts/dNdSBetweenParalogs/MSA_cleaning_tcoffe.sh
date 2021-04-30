#!/bin/bash

pwd=$(pwd)

# Scripts
module add SequenceAnalysis/MultipleSequenceAlignment/T-Coffee/11.00.8cbe486
module add Blast/ncbi-blast/2.10.1+
IndividualOGScript=${pwd}"/Scripts/dNdSBetweenParalogs/MSA_cleaning_tcoffe_PerOG.sh"

# Files & parameters
OrthologsFolder="Results/FindOrthologs"
ResultsFolder="Results/dNdSBetweenParalogs"
mkdir -p ${ResultsFolder}/MSA_cleaning_TCoffee

cd ${ResultsFolder}/MSA_cleaning_TCoffee
for group in $(cat ../Groups_wChodates.txt | grep -w 'OG_1')
do
	echo ${group}
	if [[ ! -s ${group}_AA.score_ascii ]]
	then
		sbatch -t 00:10:00 --mem=4000 -J TC${group} -o tmp/TC${group}.out -e tmp/TC${group}.err ${IndividualOGScript} ../MSA_mafft/${group}_AA.fa
		#seqScore=$(cat ${group}_AA.score_ascii | tail -n +8 | grep -v ':' | grep -v '^$' | sed 's/ \+/\t/g' | grep 'cons' | awk '{seq=seq$2}END{print seq}')
		#echo ${seqScore}
	fi
done
cd $pwd














