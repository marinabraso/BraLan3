#!/bin/bash

# Scripts
bedtools="/scratch/wally/FAC/FBM/DEE/mrobinso/default/mbrasovi/Software/bedtools.static.binary"

# Files & parameters
OrthologsFolder="Results/FindOrthologs"
ResultsFolder="Results/dNdSBetweenParalogs"
mkdir -p ${ResultsFolder}/GroupSequences

################################################
###
SpeciesOrder=(`head -1 ${OrthologsFolder}/broccoli/dir_step3/table_OGs_protein_counts.txt | awk '{for(i=2;i<=NF;i++){print $i}}' | cut -f1,2 -d'_' | cut -f1 -d'.'`)

for group in $(cut -f1 ${OrthologsFolder}/broccoli/dir_step3/table_OGs_protein_counts.txt | tail -n +2 | head -2)
do
	echo ${group}
	#if [[ ! -s ${ResultsFolder}/GroupSequences/${group}_AA.fa ]] || [[ ! -s ${ResultsFolder}/GroupSequences/${group}_DNA.fa ]]
	#then
		rm ${ResultsFolder}/GroupSequences/${group}_*.fa 2> ~/null
		for species in "${!SpeciesOrder[@]}"
		do
			echo ${SpeciesOrder[$species]}
			fileDNA=$(ls ${OrthologsFolder}/CheckedDNASequences/${SpeciesOrder[$species]}*_DNA* | head -1)
			echo $fileDNA
			if [[ $(echo ${fileDNA} | sed 's/.*\.//g') =~ "gz" ]]; then cmdReadDNA="zcat"; else  cmdReadDNA="cat"; fi
			fileAA=$(ls ${OrthologsFolder}/CheckedProteomes/${SpeciesOrder[$species]}*_AA* | head -1)
			echo $fileAA
			if [[ $(echo ${fileAA} | sed 's/.*\.//g') =~ "gz" ]]; then cmdReadAA="zcat"; else  cmdReadAA="cat"; fi

			colnum=$(( ${species}+2 )) # cut indexation starts with 1 and the first column is the group name
			for gene in $(grep -w ${group} ${OrthologsFolder}/broccoli/dir_step3/table_OGs_protein_names.txt | cut -f ${colnum} | awk -F ' ' '{for (i=1;i<=NF;i++) print $i}')
			do
				echo ${gene}
				${cmdReadDNA} $fileDNA | awk -v gene=${gene} '{if($1 ~ />/){if($1 == ">"gene){valid=1}else{valid=0}}if(valid==1){print $0}}' >> ${ResultsFolder}/GroupSequences/${group}_DNA.fa
				${cmdReadAA} $fileAA | awk -v gene=${gene} '{if($1 ~ />/){if($1 == ">"gene){valid=1}else{valid=0}}if(valid==1){print $0}}' >> ${ResultsFolder}/GroupSequences/${group}_AA.fa
				#zcat ${OrthologsFolder}/CheckedDNASequences/${SpeciesOrder[$species]}*_DNA.fa.gz | awk -v gene=${gene} '{if($1 ~ />/){if($1 == ">"gene){valid=1}else{valid=0}}if(valid==1){print $0}}' >> ${ResultsFolder}/GroupSequences/${group}_DNA.fa
				#zcat ${OrthologsFolder}/CheckedProteomes/${SpeciesOrder[$species]}*_AA.fa.gz | awk -v gene=${gene} '{if($1 ~ />/){if($1 == ">"gene){valid=1}else{valid=0}}if(valid==1){print $0}}' >> ${ResultsFolder}/GroupSequences/${group}_AA.fa
			done
		done
	#fi
done


















