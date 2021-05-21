#!/bin/bash


# Scripts
Script=Scripts/dNdSBetweenParalogs/MSA_cleaning_tcoffe.sh

# Parameters
RunName=$1
OrthologsFolder="Results/FindOrthologs/${RunName}_broccoli"
Step=10000000
Total=$(cut -f1 ${OrthologsFolder}/dir_step3/table_OGs_protein_counts.txt | wc -l)

for j in `seq 0 ${Step} ${Total}`
do
	echo $j
	sbatch -t 10:00:00 --mem=35000 -J Tcoff${j} -o tmp/Tcoff${j}.out -e tmp/Tcoff${j}.err ${Script} ${RunName} ${j} ${Step}
done






