#!/bin/bash


# Scripts
Script=Scripts/PhylogeneticTrees/MSA_TreeReconstruction.sh

# Parameters
RunName=$1
OrthologsFolder="Results/FindOrthologs/${RunName}_broccoli"
Step=100
Total=$(cut -f1 ${OrthologsFolder}/dir_step3/table_OGs_protein_counts.txt | wc -l)

for j in `seq 0 ${Step} ${Total}`
#for j in `seq 0 ${Step} 500`
do
	echo $j
	sbatch -t 10:00:00 --mem=8000 -J MSAMaff${j} -o tmp/MSAMaff${j}.out -e tmp/MSAMaff${j}.err ${Script} ${RunName} ${j} ${Step}
	#./${Script} ${RunName} ${j} ${Step}
done






