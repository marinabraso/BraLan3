#!/bin/bash

# Scripts
module add SequenceAnalysis/MultipleSequenceAlignment/mafft/7.471

# Files & parameters
OrthologsFolder="Results/FindOrthologs"
ResultsFolder="Results/dNdSBetweenParalogs"
MSAFolder=${ResultsFolder}/MSA_mafft
SeqFolder=${ResultsFolder}/GroupSequences_Chordates
mkdir -p ${MSAFolder}

for group in $(cut -f1 ${OrthologsFolder}/broccoli/dir_step3/table_OGs_protein_counts.txt | tail -n +2 | head -1)
do
	echo ${group}
	rm ${MSAFolder}/${group}_DNA.fa 

	## Alignment of AA sequences with MAFFT
	if [[ ! -s ${MSAFolder}/${group}_AA.fa ]]
	then
		mafft --globalpair --quiet --anysymbol --allowshift ${SeqFolder}/${group}_AA.fa > ${MSAFolder}/${group}_AA.fa
	fi

	## Backtranslation of aligned AA sequences to DNA sequences
	if [[ ! -s ${MSAFolder}/${group}_DNA.fa ]]
	then
		for gene in $(grep '>' ${SeqFolder}/${group}_AA.fa)
		do
			echo ${gene} >> ${MSAFolder}/${group}_DNA.fa
			echo ${gene}
			AAalnseq=$(awk -v g=${gene} '{if($1 ~ />/){if($1 == g){valid=1}else{valid=0}next;} if(valid==1){seq=seq$0}}END{print seq}' ${MSAFolder}/${group}_AA.fa)
			DNAseq=$(awk -v g=${gene} '{if($1 ~ />/){if($1 == g){valid=1}else{valid=0}next;} if(valid==1){seq=seq$0}}END{print seq}' ${SeqFolder}/${group}_DNA.fa)
			perl -le '
				my $DNAalnseq=""; # back-translated aligned sequence
				my $j=0;
				my $tmp=0;
				for(my$i=0;$i<length($ARGV[0]);$i++){
					if(substr($ARGV[0], $i, 1) =~ /-/){
						$DNAalnseq.="---";
					}else{
						$DNAalnseq.=substr($ARGV[1], $j, 3);
						$j+=3;
						$tmp++;
					}
				}
				$tmp*=3;
				print STDERR length($ARGV[0])."\t".length($DNAalnseq)."\t".$tmp."\t".$j."\t".length($ARGV[1]);
				if(length($DNAalnseq)==2772){
					print STDERR $ARGV[0]."\n".$DNAalnseq;
				}
				print $DNAalnseq;
			' -- "$AAalnseq" "$DNAseq" >> ${MSAFolder}/${group}_DNA.fa
		done
	fi
done














