#!/bin/bash

pwd=$(pwd)
export MAX_N_PID_4_TCOFFEE=10000000

# Scripts
module add tcoffee/13.45.0.4846264
module add blast-plus/2.11.0
IndividualOGScript=${pwd}"/Scripts/dNdSBetweenParalogs/MSA_cleaning_tcoffe_PerOG.sh"

# Files & parameters
RunName=$1
OrthologsFolder="Results/FindOrthologs/${RunName}_broccoli"
ResultsFolder="Results/dNdSBetweenParalogs"
mkdir -p ${ResultsFolder}/MSA_cleaning_TCoffee

cd ${ResultsFolder}/MSA_cleaning_TCoffee
for group in $(cut -f1 ${OrthologsFolder}/dir_step3/table_OGs_protein_counts.txt | tail -n +2 | head -500)
do
	echo ${group}
	if [[ ! -s ${group}_DNA_clean.fa ]]
	then
		t_coffee -infile=${pwd}/${ResultsFolder}/MSA_mafft/${group}_AA.fa  -mode evaluate  -output=score_ascii
		# Extract alignment quality score
		seqScore=$(cat ${group}_AA.score_ascii | grep 'cons' | grep -v ':' | sed 's/ \+/\t/g' | awk '{seq=seq$2}END{print seq}')
		echo ${seqScore} > ${group}_AA.score
		# each individual AA score x3 to fit DNA length
		echo ${seqScore} | awk '{for(i=1;i<=length($0);i++){l=substr($0,i,1); seq=seq""l""l""l}}END{print seq}' > ${group}_DNA.score
		# Clean DNA alignment with the alignment score only scores >5 are considered
		awk '{if(NR==FNR){for(i=1;i<=length($0);i++){a[i]=substr($0,i,1)};next} if($1 ~ />/){print $0; seq=""}else{for(i=1;i<=length($0);i++){if(a[i]>5 || substr($0,i,1)=="-"){seq=seq""substr($0,i,1)}else{seq=seq"N"}}print seq}}' ${group}_DNA.score ${pwd}/${ResultsFolder}/MSA_mafft/${group}_DNA.fa > ${group}_DNA_clean.fa
	fi
done
cd $pwd














