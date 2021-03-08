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
#for Species in Branchiostoma_lanceolatum.BraLan3 Homo_sapiens.GRCh38 Mus_musculus.GRCm39 Danio_rerio.GRCz11 Gallus_gallus.GRCg6a Branchiostoma_belcheri.Haploidv18h27 Branchiostoma_floridae.Bfl_VNyyK Strongylocentrotus_purpuratus.Spur5.0 Asterias_rubens.eAstRub1.3 Saccoglossus_kowalevskii.Skow1.1
for Species in Branchiostoma_lanceolatum.BraLan3
do
	echo ${Species}
	if [[ ! -s ${ResultsFolder}/DNASequences/${Species}.fa ]]
	then
		echo "Extracting sequences"
		Genome=$(ls ${GenomesFolder}/${Species}* | grep -v '.fai')
		format=$(echo ${Genome} | sed 's/.*\.//g')
		if [[ ${format} =~ "gz" ]]
		then
			echo ${format}
			gunzip ${Genome}
		fi
		Genome=$(ls ${GenomesFolder}/${Species}* | grep -v '.fai')
		GTF=$(ls ${OrthologsFolder}/ProteomesGTFs/${Species}*.gtf.gz)
		ProteomeBroccoli=$(ls ${OrthologsFolder}/Proteomes/${Species}*.fa)
		#rm ${ResultsFolder}/DNASequences/${Species}.fa 2> ~/null
		#for id in $(cat ${ProteomeBroccoli} | grep '>' | sed 's/>//g')
		#do
		#	echo '>'${id} >> ${ResultsFolder}/DNASequences/${Species}.fa
		#	${bedtools} getfasta -s -fi ${Genome} -bed <(grep -w ${id} <(zcat ${GTF}) | awk '{if($3 ~ /CDS/){print $1"\t"$4-1"\t"$5}}' | sort -k2,2V) | grep -v '>' | awk '{seq=seq$0}END{print seq}'  >> ${ResultsFolder}/DNASequences/${Species}.fa
		#done
		#gzip ${ResultsFolder}/DNASequences/${Species}.fa
	fi
	echo "Comparing DNA and AA sequences"
	./Scripts/dNdSBetweenParalogs/Check_DNA_AA_seq.pl --DNAfasta=${ResultsFolder}/DNASequences/${Species}.fa --AAfasta=${ProteomeBroccoli} --output=${ResultsFolder}/Checking_DNA_AA_sequences/${Species}_DNA_AA_seq_check.txt
done














