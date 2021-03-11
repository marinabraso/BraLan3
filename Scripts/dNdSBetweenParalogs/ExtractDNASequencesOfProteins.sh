#!/bin/bash

# Scripts
bedtools="/scratch/wally/FAC/FBM/DEE/mrobinso/default/mbrasovi/Software/bedtools.static.binary"

# Files & parameters
OrthologsFolder="Results/FindOrthologs"
ResultsFolder="Results/dNdSBetweenParalogs"
mkdir -p ${ResultsFolder}/DNASequences
mkdir -p ${ResultsFolder}/OrthologousGroupsSequences
mkdir -p ${ResultsFolder}/Checking_DNA_AA_sequences
GenomesFolder="Data/Genomes"

################################################
###
for Species in Branchiostoma_lanceolatum.BraLan3 Homo_sapiens.GRCh38 Mus_musculus.GRCm39 Danio_rerio.GRCz11 Gallus_gallus.GRCg6a Branchiostoma_belcheri.Haploidv18h27 Branchiostoma_floridae.Bfl_VNyyK Strongylocentrotus_purpuratus.Spur5.0 Asterias_rubens.eAstRub1.3 Saccoglossus_kowalevskii.Skow1.1
#for Species in Homo_sapiens.GRCh38
do
	echo ${Species}
	Genome=$(ls ${GenomesFolder}/${Species}* | grep -v '.fai')
	GTF=$(ls ${OrthologsFolder}/ProteomesGTFs/${Species}*.gtf.gz)
	ProteomeBroccoli=$(ls ${OrthologsFolder}/Proteomes/${Species}*)
	format=$(echo ${ProteomeBroccoli} | sed 's/.*\.//g')
	if [[ ${format} =~ "gz" ]];	then gunzip ${ProteomeBroccoli}; fi
	ProteomeBroccoli=$(ls ${OrthologsFolder}/Proteomes/${Species}*.fa)

	# Extract DNA sequences of proteins from the genome reading GTF
	if [[ ! -s ${ResultsFolder}/DNASequences/${Species}.fa ]] && [[ ! -s ${ResultsFolder}/DNASequences/${Species}.fa.gz ]] 
	then
		echo "Extracting sequences"
		format=$(echo ${Genome} | sed 's/.*\.//g')
		if [[ ${format} =~ "gz" ]];	then gunzip ${Genome}; fi
		Genome=$(ls ${GenomesFolder}/${Species}* | grep -v '.fai')
		rm ${ResultsFolder}/DNASequences/${Species}.fa 2> ~/null
		for id in $(cat ${ProteomeBroccoli} | grep '>' | sed 's/>//g')
		do
			echo '>'${id} >> ${ResultsFolder}/DNASequences/${Species}.fa
			# Use bedtoos getfasta to extract sequences in a bed format (greped and formated from the GTF)
			PrimarySeq=$(${bedtools} getfasta -fi ${Genome} -bed <(grep -w ${id} <(zcat ${GTF}) | awk '{if($3 ~ /CDS/){print $1"\t"$4-1"\t"$5}}' | sort -k2,2V) | grep -v '>' | awk '{seq=seq$0}END{print seq}')
			Strand=$(zcat ${GTF} | grep -w ${id} | head -1 | cut -f7)
			if [[ ${Strand} == "+" ]]
			then
				echo ${PrimarySeq}	>> ${ResultsFolder}/DNASequences/${Species}.fa
			else
				# Reverse complement
				echo ${PrimarySeq} | rev | sed 's/A/Z/g' | sed 's/T/A/g' | sed 's/Z/T/g' | sed 's/C/Z/g' | sed 's/G/C/g' | sed 's/Z/G/g' >> ${ResultsFolder}/DNASequences/${Species}.fa
			fi
		done
	fi

	# Check that extracted DNA sequences correspond to protein sequences
	echo "Comparing DNA and AA sequences"
	rm ${ResultsFolder}/Checking_DNA_AA_sequences/${Species}_DNA_AA_seq_check.txt 2> ~/null
	for id in $(cat ${ProteomeBroccoli} | grep '>' | sed 's/>//g')
	do
		#echo ${id}
		ProtSeq=$(awk -v id=${id} '{if($1 ~ />/){if($1 == ">"id){valid=1;next;}else{valid=0}} if(valid == 1){print $0}}' ${ProteomeBroccoli} | awk '{seq=seq""$1}END{print seq}')
		cDNASeq=$(awk -v id=${id} '{if($1 ~ />/){if($1 == ">"id){valid=1;next;}else{valid=0}} if(valid == 1){print $0}}' ${ResultsFolder}/DNASequences/${Species}.fa | awk '{seq=seq""$1}END{print seq}')
		chr=$(zcat ${GTF} | grep -w ${id} | head -1 | cut -f1)
		perl -le '
			my %GeneticCode;
	    
			my   @AAs=qw/F F L L S S S S Y Y 1 1 C C 1 W L L L L P P P P H H Q Q R R R R I I I M T T T T N N K K S S R R V V V V A A A A D D E E G G G G/;
			my @Base1=qw/T T T T T T T T T T T T T T T T C C C C C C C C C C C C C C C C A A A A A A A A A A A A A A A A G G G G G G G G G G G G G G G G/;
			my @Base2=qw/T T T T C C C C A A A A G G G G T T T T C C C C A A A A G G G G T T T T C C C C A A A A G G G G T T T T C C C C A A A A G G G G/;
			my @Base3=qw/T C A G T C A G T C A G T C A G T C A G T C A G T C A G T C A G T C A G T C A G T C A G T C A G T C A G T C A G T C A G T C A G/;

			for(my$c=0;$c<=scalar(@AAs);$c++){
				$GeneticCode{$Base1[$c].$Base2[$c].$Base3[$c]}=$AAs[$c];
			}

			if($ARGV[2]=="MT"){
				$GeneticCode{"AGA"}="1";
				$GeneticCode{"AGG"}="1";
				$GeneticCode{"ATA"}="M";
				$GeneticCode{"TGA"}="W";
			}

			my $BTseq=""; # back-translated sequence
			for(my$c=0;$c<=length($ARGV[1]);$c+=3){
				$BTseq.=$GeneticCode{substr($ARGV[1], $c, 3)};
			}

			if($ARGV[0]==$BTseq){
				print $ARGV[3]."\tIdentical";
			}else{
				print $ARGV[3]."\tError\t".$ARGV[0]."\t".$BTseq."\t".$ARGV[1];
			}
		' -- "$ProtSeq" "$cDNASeq" "$chr" "$id" >> ${ResultsFolder}/Checking_DNA_AA_sequences/${Species}_DNA_AA_seq_check.txt
	done
done




























