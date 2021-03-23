#!/bin/bash

# Scripts
module add SequenceAnalysis/MultipleSequenceAlignment/mafft/7.471

# Files & parameters
OrthologsFolder="Results/FindOrthologs"
ResultsFolder="Results/dNdSBetweenParalogs"
mkdir -p ${ResultsFolder}/MSA_mafft

for group in $(cut -f1 ${OrthologsFolder}/broccoli/dir_step3/table_OGs_protein_counts.txt | tail -n +2 | head -2)
do
	echo ${group}
	#if [[ ! -s ${ResultsFolder}/MSA_mafft/${group}_AA.fa ]] || [[ ! -s ${ResultsFolder}/MSA_mafft/${group}_DNA.fa ]]
	#then
		rm ${ResultsFolder}/MSA_mafft/${group}_*.fa 2> ~/null
		mafft --globalpair --quiet --anysymbol --allowshift ${ResultsFolder}/GroupSequences/${group}_AA.fa > ${ResultsFolder}/MSA_mafft/${group}_AA.fa

		for gene in $(grep '>' ${ResultsFolder}/GroupSequences/${group}_AA.fa)
		do
			echo ${gene}
			echo ${gene} >> ${ResultsFolder}/MSA_mafft/${group}_DNA.fa
			AAalnseq=$(awk -v g=gene '{if($1 ~ />/){if($1 == ">"g){valif=1}else{valid=0}next;} if(valid==1){print $0}}' ${ResultsFolder}/MSA_mafft/${group}_AA.fa)
			echo ${AAalnseq}
			DNAseq=$(awk -v g=gene '{if($1 ~ />/){if($1 == ">"g){valif=1}else{valid=0}next;} if(valid==1){print $0}}' ${ResultsFolder}/GroupSequences/${group}_DNA.fa)
			echo ${DNAseq}

			perl -le '
				my $DNAalnseq=""; # back-translated aligned sequence
				my $j=0;
				for(my$i=0;$i<=length($ARGV[0]);$i++){
					print substr($ARGV[0], $i, 1)."  ";
					if(substr($ARGV[0], $i, 1) == "-"){
						$DNAalnseq.="---";
						print "---";
					}else{
						$DNAalnseq.=substr($ARGV[1], $j, 3)};
						print substr($ARGV[1], $j, 3);
						$j+=3;
					}
				}
				print $DNAalnseq;
			' -- "$AAalnseq" "$DNAseq" >> ${ResultsFolder}/MSA_mafft/${group}_DNA.fa

		done
	#fi
done














