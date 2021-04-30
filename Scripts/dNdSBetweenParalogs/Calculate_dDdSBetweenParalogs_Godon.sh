#!/bin/bash

module add Phylogeny/godon/20200613
module add SequenceAnalysis/MultipleSequenceAlignment/mafft/7.471

# Scripts
guidance=../Software/guidance.v2.02/www/Guidance/guidance.pl

# Files & parameters
OrthologsFolder="Results/FindOrthologs"
ResultsFolder="Results/dNdSBetweenParalogs"
mkdir -p ${ResultsFolder}/Godon_M0 ${ResultsFolder}/Godon_M8 ${ResultsFolder}/Guidance



#perl ${guidance} --seqFile ${ResultsFolder}/GroupSequences_Chordates/OG_1_DNA.fa --msaProgram MAFFT --seqType codon --outDir ${ResultsFolder}/Guidance --msaFile ${ResultsFolder}/MSA_mafft/OG_1_DNA.fa

################################################
### Recalculate branch lengths with M0
for group in $(cat ${ResultsFolder}/Groups_wChodates.txt)
do
	echo ${group}
	if [[ -s ${ResultsFolder}/TreeReconstruction/RAxML_result.${group}_GTRGAMMA ]]
	then
		if [[ ! -s ${ResultsFolder}/Godon_M8/${group}_M8_likelihood_ratio.txt ]]
		then
			godon test M8 --m0-tree ${ResultsFolder}/MSA_mafft/${group}_DNA.fa ${ResultsFolder}/TreeReconstruction/RAxML_result.${group}_GTRGAMMA > ${ResultsFolder}/Godon_M8/${group}_M8_likelihood_ratio.txt
		fi
	else
		echo ${group} > ${ResultsFolder}/MissingTrees.txt
	fi
done







