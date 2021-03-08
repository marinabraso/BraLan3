#!/bin/bash

# Scripts
bedtools="/scratch/wally/FAC/FBM/DEE/mrobinso/default/mbrasovi/Software/bedtools.static.binary"

# Files & parameters
OrthologsFolder="Results/FindOrthologs"
ResultsFolder="Results/dNdSBetweenParalogs"
mkdir -p ${ResultsFolder}/Sequences
GenomesFolder="Data/Genomes"
GTFFolder="Data/Transcriptomes"

################################################
###
gunzip ${GenomesFolder}/*gz 2> ~/null
SpeciesGenomes=()
SpeciesGTF=()
SpeciesOrder=(`head -1 Results/FindOrthologs/broccoli/dir_step3/table_OGs_protein_counts.txt | awk '{for(i=2;i<=NF;i++){print $i}}' | cut -f1,2 -d'_' | cut -f1 -d'.'`)
for species in "${!SpeciesOrder[@]}"
	do
	Genome=$(ls ${GenomesFolder}/${SpeciesOrder[${species}]}*)
	SpeciesGenomes+=(${Genome})
	GTF=$(ls ${GTFFolder}/${SpeciesOrder[${species}]}*)
	SpeciesGTF+=(${GTF})
done

for group in $(cut -f1 ${OrthologsFolder}/broccoli/dir_step3/table_OGs_protein_counts.txt | tail -n +2 | head -2)
do
	echo ${group}
	#if [[ ! -s ${ResultsFolder}/Sequences/${group}_dna.fa ]]
	#then
		rm ${ResultsFolder}/Sequences/${group}_dna.fa 2> ~/null
		for species in "${!SpeciesOrder[@]}"
		do
			#echo ${species}" "${SpeciesOrder[${species}]}" "${SpeciesGenomes[${species}]}" "${SpeciesGTF[${species}]}
			colnum=$(( ${species}+2 )) # cut indexation starts with 1 and the first column is the group name
			genes=$(grep -w ${group} ${OrthologsFolder}/broccoli/dir_step3/table_OGs_protein_names.txt | cut -f ${colnum})
			IFS=' ' read -r -a GenesA <<< "$genes"
			for gene in "${GenesA[@]}"
			do
				echo '>'${gene} >> ${ResultsFolder}/Sequences/${group}_dna.fa
				${bedtools} getfasta -s -fi ${SpeciesGenomes[${species}]} -bed <(grep ${gene} <(zcat ${SpeciesGTF[${species}]}) | awk '$3 ~ /CDS/' | sort -k2,2V) | grep -v '>' | awk '{seq=seq$0}END{print seq}'  >> ${ResultsFolder}/Sequences/${group}_dna.fa
			done
		done
	#fi
done


















