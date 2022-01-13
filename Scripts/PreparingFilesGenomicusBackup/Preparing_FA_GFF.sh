#!/bin/bash

# Scripts
bedtools="../Software/bedtools.static.binary"

# Files & parameters
Species="Branchiostoma_lanceolatum.BraLan3"
ResultsFolder="Results/PreparingFilesGenomicus"
mkdir -p ${ResultsFolder}
GTF=Data/Transcriptomes/Branchiostoma_lanceolatum.BraLan3.gtf.gz
Genome=Data/Genomes/Branchiostoma_lanceolatum.BraLan3_genome.fa


# GTF to GFF
if [[ ! -s ${ResultsFolder}/${Species}.gff ]]
then
	echo "Extracting gff"
	zcat ${GTF} | grep 'strong' | sed 's/gene_id "\([A-Za-z0-9]\+\)".*transcript_id "\([A-Za-z0-9_]\+\)".*gene_name "\([A-Za-z0-9\.\-]\+\)".*protein_id "\([A-Za-z0-9_]\+\)".*/gene_id=\1;transcript_id=\2;gene_name=\3;protein_id=\4/g' > ${ResultsFolder}/${Species}.gff
fi

# GTF to CDS GFF
if [[ ! -s ${ResultsFolder}/${Species}_CDS.gff ]]
then
	echo "Extracting gff"
	zcat ${GTF} | awk -F'\t' '{if($3 ~ /CDS/){print $0}}' | grep 'strong' | sed 's/gene_id "\([A-Za-z0-9]\+\)".*transcript_id "\([A-Za-z0-9_]\+\)".*gene_name "\([A-Za-z0-9\.\-]\+\)".*protein_id "\([A-Za-z0-9_]\+\)".*/gene_id=\1;transcript_id=\2;gene_name=\3;protein_id=\4/g' > ${ResultsFolder}/${Species}_CDS.gff
fi

# Extract DNA sequences
echo "Extracting sequences"
rm ${ResultsFolder}/${Species}_CDS_3.fa
if [[ ! -s ${ResultsFolder}/${Species}_CDS_3.fa ]]
then
	for id in $(cut -f9 ${ResultsFolder}/${Species}_CDS.gff | sort -u)
	do
		echo '>'${id} >> ${ResultsFolder}/${Species}_CDS_3.fa
		# Use bedtoos getfasta to extract sequences
		PrimarySeq=$(${bedtools} getfasta -fi ${Genome} -bed <(grep -w ${id} ${ResultsFolder}/${Species}_CDS.gff | sort -k4,4V) | grep -v '>' | awk '{seq=seq$0}END{print seq}')
		Strand=$(cat ${ResultsFolder}/${Species}_CDS.gff | grep -w ${id} | head -1 | cut -f7)
		if [[ ${Strand} == "+" ]]
		then
			echo ${PrimarySeq}	>> ${ResultsFolder}/${Species}_CDS_3.fa
		else
			# Reverse complement
			echo ${PrimarySeq} | rev | sed 's/A/Z/g' | sed 's/T/A/g' | sed 's/Z/T/g' | sed 's/C/Z/g' | sed 's/G/C/g' | sed 's/Z/G/g' >> ${ResultsFolder}/${Species}_CDS_3.fa
		fi
	done
fi




#gzip ${ResultsFolder}/Proteomes/*.fa ${ResultsFolder}/ProteomesGTFs/*.gtf



























