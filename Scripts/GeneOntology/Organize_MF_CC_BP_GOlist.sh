#!/bin/bash

# Scripts

# Files & parameters
OBO="Data/GeneOntology/go.obo"
ENSAnnot="Data/GeneOntology/HumanGeneToGoTerm_Ensembl.txt"
ResultsFolder="Results/GeneOntology"
GOlistFile=${ResultsFolder}"/GOlist_MF.list"
rm ${GOlistFile} 2> ~/null

################################################
### Molecular function
BigGOterms=( GO:0005198 GO:0016887 GO:0016209 GO:0003824 GO:0140110 GO:0045182 GO:0005215 GO:0038024 GO:0004857 GO:0008047 GO:0019207 GO:0099106 GO:0098772 GO:0030546 GO:0019888 GO:0044092 GO:0044093 GO:0065009 GO:0005488 GO:0060090 GO:0140104 GO:0060089 )
BigGOclass=( "Structural proteins" "Enzymes" "Enzymes" "Enzymes" "Transcription factors" "Translation factors" "Transporters" "Transporters" "Regulators" "Regulators" "Regulators" "Regulators" "Regulators" "Regulators" "Regulators" "Regulators" "Regulators" "Regulators" "Binding" "Binding" "Carriers" "Signaling" )

for bGO in ${!BigGOterms[@]}; do
	echo ${BigGOterms[${bGO}]} ${BigGOclass[${bGO}]}

	# For each big GO term retrieve its childs and create its gene file
	ChildGOterms=()
	while IFS= read -r line; do
		ChildGOterms+=( "$line" )
	done < <( ./Scripts/GeneOntology/GetAllGenesFromTerm_from_OBO_ENS.sh ${BigGOterms[${bGO}]} ${ResultsFolder}"/"${BigGOterms[${bGO}]}".txt" ${OBO} ${ENSAnnot} )

	# For each child GO term create its gene file
	for cGO in ${!ChildGOterms[@]}; do
		if [[ ! -e ${ResultsFolder}"/"${ChildGOterms[${cGO}]}".txt" ]]; then
			echo  ${ResultsFolder}"/"${ChildGOterms[${cGO}]}".txt"
			./Scripts/GeneOntology/GetAllGenesFromTerm_from_OBO_ENS.sh ${ChildGOterms[${cGO}]} ${ResultsFolder}"/"${ChildGOterms[${cGO}]}".txt" ${OBO} ${ENSAnnot} > ~/null
		fi
	done

	# For each child GO term retrieve its name and print it
	ChildGOnames=()
	while IFS= read -r line; do
		ChildGOnames+=( "$line" )
	done < <( awk -F '\t' '{if(NR==FNR){a[$1]=1;next} if($1 ~ /^id: /){id=substr($1, 5, length($1));next} if(a[id] && $1 ~ /^name: /){name=substr($1, 7, length($1)); print name}}' <(for cGO in ${!ChildGOterms[@]}; do echo ${ChildGOterms[${cGO}]}; done) ${OBO} )

	for cGO in ${!ChildGOterms[@]}; do
		echo ${ChildGOterms[${cGO}]}"	"${BigGOclass[${bGO}]}"	"${ChildGOnames[${cGO}]} >> ${GOlistFile}
	done

done









