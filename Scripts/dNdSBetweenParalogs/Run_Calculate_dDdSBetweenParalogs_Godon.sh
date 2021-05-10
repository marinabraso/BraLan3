#!/bin/bash


# Scripts
Script=Scripts/dNdSBetweenParalogs/Calculate_dDdSBetweenParalogs_Godon.sh

# Parameters
RunName=$1
OrthologsFolder="Results/FindOrthologs/${RunName}_broccoli"
Step=50
Total=$(cut -f1 ${OrthologsFolder}/dir_step3/table_OGs_protein_counts.txt | wc -l)

for j in `seq 250 ${Step} ${Total}`
do
	echo $j
	sbatch -t 10:00:00 --mem=8000 -J GodM8${j} -o tmp/GodM8${j}.out -e tmp/GodM8${j}.err ${Script} ${RunName} ${j} ${Step}
done




