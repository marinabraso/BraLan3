#!/bin/bash

# Scripts

# Files & parameters
RunName=$1
SequencesFolder="Results/FilteringGeneSets"
OrthologsFolder="Results/FindOrthologs/${RunName}_broccoli"
BaseResultsFolder="Results/PhylogeneticTrees"
ResultsFolder="${BaseResultsFolder}/OGSequences_${RunName}"
mkdir -p ${ResultsFolder}
SelectedSpecies=("Branchiostoma_floridae" "Branchiostoma_lanceolatum" "Branchiostoma_belcheri" "Danio_rerio" "Mus_musculus" "Homo_sapiens" "Gallus_gallus")

################################################
### For each orthologous grup extract AA and DNA sequences
SpeciesOrder=(`head -1 ${OrthologsFolder}/dir_step3/table_OGs_protein_counts.txt | awk '{for(i=2;i<=NF;i++){print $i}}' | cut -f1,2 -d'_' | cut -f1 -d'.'`)
for group in $(cut -f1 ${OrthologsFolder}/dir_step3/table_OGs_protein_counts.txt | tail -n +2)
do
	echo ${group}
	if [[ ! -s ${ResultsFolder}/${group}_AA.fa ]] || [[ ! -s ${ResultsFolder}/${group}_DNA.fa ]]
	then
		rm ${ResultsFolder}/${group}_*.fa 2> ~/null # in case there was previous runnings
		withchordates=0
		for species in "${!SpeciesOrder[@]}"
		do
			if [[ " ${SelectedSpecies[@]} " =~ " ${SpeciesOrder[$species]} " ]]
			then
				fileDNA=$(ls ${SequencesFolder}/DNASequences/${SpeciesOrder[$species]}*_DNA* | head -1)
				fileAA=$(ls ${SequencesFolder}/Proteomes/${SpeciesOrder[$species]}* | head -1)
				colnum=$(( ${species}+2 )) # cut indexation starts with 1 and the first column is the group name
				for gene in $(grep -w ${group} ${OrthologsFolder}/dir_step3/table_OGs_protein_names.txt | cut -f ${colnum} | awk -F ' ' '{for (i=1;i<=NF;i++) print $i}')
				do
					$withchordates=1
					grep -A 1 ${gene} $fileDNA >> ${ResultsFolder}/${group}_DNA.fa
					grep -A 1 ${gene} $fileAA >> ${ResultsFolder}/${group}_AA.fa
				done
			fi
		done
		if [[ $withchordates -eq 1 ]]
		then
			echo ${group} >> ${BaseResultsFolder}/OG_wChodates_${RunName}.txt
		fi
	fi
done
















