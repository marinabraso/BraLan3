#!/bin/bash


# Scripts
Script=Scripts/dNdS/MSA_AA_backtranslation_DNA.sh

# Parameters
RunName=$1
OrthologsFolder="Results/FindOrthologs/${RunName}_broccoli"
Step=100
Total=$(cut -f1 ${OrthologsFolder}/dir_step3/table_OGs_protein_counts.txt | wc -l)

for j in `seq 0 ${Step} ${Total}`
do
	echo $j
	sbatch -t 10:00:00 --mem=8000 -J MSAMaff${j} -o tmp/MSAMaff${j}.out -e tmp/MSAMaff${j}.err ${Script} ${RunName} ${j} ${Step}
done






