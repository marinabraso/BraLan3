#!/bin/bash

# Scripts

# Files & parameters
ProteomesFolder="Data/Proteomes"
TranscriptomesFolder="Data/Transcriptomes"
ResultsFolder="Results/FindOrthologs"
mkdir -p ${ResultsFolder}/Proteomes
mkdir -p ${ResultsFolder}/ProteomesGTFs


################################################
### Branchiostoma lanceolatum
echo "Branchiostoma_lanceolatum.BraLan3"
if [[ ! -s ${ResultsFolder}/Proteomes/Branchiostoma_lanceolatum.BraLan3_ProteinCoding_StrongEvidence.fa.gz ]] || [[ ! -s ${ResultsFolder}/ProteomesGTFs/Branchiostoma_lanceolatum.BraLan3_ProteinCoding_StrongEvidence.gtf.gz ]]
then
	# List of protein coding genes with "strong" evidence in any of the testedevidencies
	zcat ${TranscriptomesFolder}/Branchiostoma_lanceolatum.BraLan3.gtf.gz | awk '{if($3 ~ /CDS/){print $0}}' | grep 'protein_coding' | cut -f9 | sed 's/gene_id \"\([A-Z0-9]\+\)\"; transcript_id.*bgee_evidence \"\([a-z]\+\)\"; uniprot_evidence \"\([a-z]\+\)\"; bflo_evidence \"\([a-z]\+\)\"; bbel_evidence \"\([a-z]\+\)\";/\1\t\2\t\3\t\4\t\5/g' |  sort | uniq | grep 'strong' | cut -f1 > ${ResultsFolder}/BraLan3_ProteinCodingGenes_StrongEvidence.txt
	# Select the fasta sequences of those in the list & substitute . by * in fasta sequences (?)
	awk '{if(NR==FNR){a[">"$1]=1;next;} if(a[$1]==1){valid=1;}else{if($1 ~ />/){valid=0;}} if(valid==1){print $0}}' ${ResultsFolder}/BraLan3_ProteinCodingGenes_StrongEvidence.txt <(zcat ${ProteomesFolder}/Branchiostoma_lanceolatum.BraLan3.fa.gz) | sed 's/\./*/g' > ${ResultsFolder}/Proteomes/Branchiostoma_lanceolatum.BraLan3_ProteinCoding_StrongEvidence.fa
	# Extract GTF info of the selected genes
	awk '{if(NR==FNR){a[$1]=1;next;} if(a[$9]==1){print $0}}' <(cat ${ResultsFolder}/Proteomes/Branchiostoma_lanceolatum.BraLan3_ProteinCoding_StrongEvidence.fa | grep '>' | sed 's/>//g') <(zcat ${TranscriptomesFolder}/Branchiostoma_lanceolatum.BraLan3.gtf.gz | awk '{if($3 == "CDS"){print $0}}' | sed 's/gene_id "\([A-Z0-9]\+\)".*/\1/g') > ${ResultsFolder}/ProteomesGTFs/Branchiostoma_lanceolatum.BraLan3_ProteinCoding_StrongEvidence.gtf
fi

################################################
### Species from ENSEMBL 103: 
#		- Homo sapiens
# 		- Mus musculus
# 		- Gallus gallus
# 		- Danio rerio
for Species in Homo_sapiens.GRCh38 Mus_musculus.GRCm39 Danio_rerio.GRCz11 Gallus_gallus.GRCg6a
do
	echo ${Species}
	if [[ ! -s ${ResultsFolder}/Proteomes/${Species}_ProteinCoding_LongestTx_inGTF.fa.gz ]] || [[ ! -s ${ResultsFolder}/ProteomesGTFs/${Species}_ProteinCoding_LongestTx_inGTF.gtf.gz ]]
	then
		# List all the protein coding genes with CDS annotation in the GTF
		zcat ${TranscriptomesFolder}/${Species}.103.gtf.gz | grep 'protein_coding' | awk '{if($3 == "CDS"){print $0}}' | cut -f9 | sed 's/gene_id "\([A-Z0-9]\+\)".*transcript_id "\([A-Z0-9]\+\)".*/\1\t\2/g' | sort -u > ${ResultsFolder}/Proteomes/${Species}_ProteinCoding_AllTx_inGTF.txt
		# Extract the protein sequences properly annotated as protein coding and adapt the sequence names
		zcat ${ProteomesFolder}/${Species}.pep.all.fa.gz | awk -F ' ' '{if($1 ~ />/){if($6 == "gene_biotype:protein_coding"){print ">"$4"|"$5; valid=1;}else{valid=0}}else{if(valid==1){print $0}}}' | sed 's/gene://g' | sed 's/transcript://g' | sed 's/\.[0-9]//g' > ${ResultsFolder}/Proteomes/${Species}_ProteinCoding_AllTx.fa
		# Extract the longest tx from the ones present in the GTF file
		awk '{if(NR==FNR){a[$1]=$2;next;} if(a[$1]==$2){print $0}}' ${ResultsFolder}/Proteomes/${Species}_ProteinCoding_AllTx_inGTF.txt <(awk '{if($1 ~ />/){print name"\t"len; name=$1; len=0; next;} len=len+length($1);}END{print name"\t"len;}' ${ResultsFolder}/Proteomes/${Species}_ProteinCoding_AllTx.fa | tail -n +2 | sed 's/>//g' | sed 's/|/\t/g') | sort -k1,1 -k3,3V | awk '{if(g==$1){t=$2}else{print g"\t"t;g=$1;t=$2}}END{print g"\t"t;}' | tail -n +2 > ${ResultsFolder}/${Species}_ProteinCoding_LongestTx_inGTF.txt
		# Extract sequences from the selected longest tx
		awk '{if(NR==FNR){a[">"$1]=$2;next;} if($1 ~ />/){if(a[$1]==$2){valid=1}else{valid=0}} if(valid==1){print $0}}' ${ResultsFolder}/${Species}_ProteinCoding_LongestTx_inGTF.txt <(sed 's/|/\t/g' ${ResultsFolder}/Proteomes/${Species}_ProteinCoding_AllTx.fa) | sed 's/\t/|/g' > ${ResultsFolder}/Proteomes/${Species}_ProteinCoding_LongestTx_inGTF.fa
		# Remove intermediate files
		rm ${ResultsFolder}/Proteomes/${Species}_ProteinCoding_AllTx.fa ${ResultsFolder}/Proteomes/${Species}_ProteinCoding_AllTx_inGTF.txt
		# Extract GTF info of the selected genes
		awk '{if(NR==FNR){a[$1]=1;next;} if(a[$9]==1){print $0}}' <(cat ${ResultsFolder}/Proteomes/${Species}_ProteinCoding_LongestTx_inGTF.fa | grep '>' | sed 's/>//g') <(zcat Data/Transcriptomes/${Species}.103.gtf.gz | awk '{if($3 == "CDS"){print $0}}' | sed 's/gene_id "\([A-Z0-9]\+\)".*transcript_id "\([A-Z0-9]\+\)".*/\1|\2\t\1\t\2/g') > ${ResultsFolder}/ProteomesGTFs/${Species}_ProteinCoding_LongestTx_inGTF.gtf
	fi
done

################################################
### Species from NCBI Genome Database: 
# 		- Branchiostoma floridae
# 		- Branchiostoma belcheri
# 		- Strongylocentrotus purpuratus
#		- Asterias rubens
#		- Saccoglossus kowalevskii

for Species in Branchiostoma_belcheri.Haploidv18h27 Branchiostoma_floridae.Bfl_VNyyK Strongylocentrotus_purpuratus.Spur5.0 Asterias_rubens.eAstRub1.3 Saccoglossus_kowalevskii.Skow1.1
do
	echo ${Species}
	if [[ ! -s ${ResultsFolder}/Proteomes/${Species}_ProteinCoding_inGTF_CleanFormat.fa.gz ]] || [[ ! -s ${ResultsFolder}/ProteomesGTFs/${Species}_ProteinCoding_inGTF_CleanFormat_LongestProt.gtf.gz ]]
	then
		SpeciesShortName=$(echo ${Species^^} | awk -F '_' '{print substr($1, 1, 1)substr($2, 1, 3)}')
		# Extract info from GFF (filter for CDS & clean format of gene ID and protein ID)
		zcat ${TranscriptomesFolder}/${Species}.gff.gz | grep -v '^#' | awk '{if($3=="CDS"){print $0}}' | sed 's/ID.*;gene=\([A-Za-z0-9\.\/\-]\+\);.*;protein_id=\([A-Za-z0-9_\.]\+\)$/\1\t\2/g' | sed 's/ID.*;gene=\([A-Za-z0-9\.\/\-]\+\);.*;protein_id=\([A-Za-z0-9_\.]\+\);.*/\1\t\2/g' | awk -v sn=${SpeciesShortName} '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"sn"_"$9"\t"sn"_"$10}' > ${ResultsFolder}/${Species}_ProteinCoding_inGTF_CleanFormat.txt
		# Extract protein sequences present in filtered gff and add geneID to protID
		awk '{if(NR==FNR){a[">"$10]=1;b[">"$10]=$9;next} if($1 ~ /^>/){if(a[$1]==1){valid=1; print $1"|"b[$1]}else{valid=0}}else{if(valid==1){print $0}}}' ${ResultsFolder}/${Species}_ProteinCoding_inGTF_CleanFormat.txt <(zcat ${ProteomesFolder}/${Species}.faa.gz | awk -v sn=${SpeciesShortName} '{if($1 ~ /^>/){print ">"sn"_"substr($1, 2, length($1))}else{print $0}}') > ${ResultsFolder}/Proteomes/${Species}_ProteinCoding_inGTF_CleanFormat.fa
		# Extract longest protein for each gene
		awk '{if($1 ~ />/){print name"\t"len; name=$1; len=0; next;} len=len+length($1);}END{print name"\t"len;}' ${ResultsFolder}/Proteomes/${Species}_ProteinCoding_inGTF_CleanFormat.fa | tail -n +2 | sed 's/>//g' | sed 's/|/\t/g' | sort -k2,2 -k3,3V | awk '{if(g==$2){p=$1}else{print g"\t"p;g=$2;p=$1}}END{print g"\t"p;}' | tail -n +2 > ${ResultsFolder}/${Species}_ProteinCoding_inGTF_CleanFormat_LongestProt.txt
		# Extract sequences of the longest protein per gene
		awk '{if(NR==FNR){a[">"$2"|"$1]=1;next;} if($1 ~ /^>/){if(a[$1]==1){valid=1}else{valid=0}}if(valid==1){print $0}}' ${ResultsFolder}/${Species}_ProteinCoding_inGTF_CleanFormat_LongestProt.txt ${ResultsFolder}/Proteomes/${Species}_ProteinCoding_inGTF_CleanFormat.fa > ${ResultsFolder}/Proteomes/${Species}_ProteinCoding_inGTF_CleanFormat_LongestProt.fa
		# Remove unnecessary intermediate files
		rm ${ResultsFolder}/Proteomes/${Species}_ProteinCoding_inGTF_CleanFormat.fa
		# Extract GTF info of the selected genes
		awk '{if(NR==FNR){a[$1]=1;next;} if(a[$9]==1){print $0}}' <(cat ${ResultsFolder}/Proteomes/${Species}_ProteinCoding_inGTF_CleanFormat_LongestProt.fa | grep '>' | sed 's/>//g') <(cat ${ResultsFolder}/${Species}_ProteinCoding_inGTF_CleanFormat.txt | awk '{if($3 == "CDS"){print $0}}' | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$10"|"$9"\t"$9"\t"$10}') > ${ResultsFolder}/ProteomesGTFs/${Species}_ProteinCoding_inGTF_CleanFormat_LongestProt.gtf
	fi
done

gzip ${ResultsFolder}/Proteomes/*.fa ${ResultsFolder}/ProteomesGTFs/*.gtf




