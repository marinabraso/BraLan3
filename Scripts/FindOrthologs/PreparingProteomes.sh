#!/bin/bash

# Scripts

# Files & parameters
ProteomesFolder="Data/Proteomes"
TranscriptomesFolder="Data/Transcriptomes"
ResultsFolder="Results/FindOrthologs"
mkdir -p ${ResultsFolder}/Proteomes

################################################
### Branchiostoma lanceolatum

# List of protein coding genes with "strong" evidence in any of the testedevidencies
zcat ${TranscriptomesFolder}/BranchiostomaLanceolatum_BraLan3.gtf.gz | awk '{if($3 ~ /CDS/){print $0}}' | grep 'protein_coding' | cut -f9 | sed 's/gene_id \"\([A-Z0-9]\+\)\"; transcript_id.*bgee_evidence \"\([a-z]\+\)\"; uniprot_evidence \"\([a-z]\+\)\"; bflo_evidence \"\([a-z]\+\)\"; bbel_evidence \"\([a-z]\+\)\";/\1\t\2\t\3\t\4\t\5/g' |  sort | uniq | grep 'strong' | cut -f1 > ${ResultsFolder}/BraLan3_ProteinCodingGenes_StrongEvidence.txt
# Select the fasta sequences of those in the list & sibstitute . by * in fasta sequences (?)
awk '{if(NR==FNR){a[">"$1]=1;next;} if(a[$1]==1){valid=1;}else{if($1 ~ />/){valid=0;}} if(valid==1){print $0}}' ${ResultsFolder}/BraLan3_ProteinCodingGenes_StrongEvidence.txt <(zcat ${ProteomesFolder}/Branchiostoma_lanceolatum.BraLan3.fa.gz) | sed 's/\./*/g' > ${ResultsFolder}/Proteomes/Branchiostoma_lanceolatum.BraLan3_ProteinCoding_StrongEvidence.fa

################################################
### Species from ENSEMBL 103: 
#		- Homo sapiens
# 		- Mus musculus
# 		- Gallus gallus
# 		- Danio rerio
for Species in Homo_sapiens.GRCh38 Mus_musculus.GRCm39 Danio_rerio.GRCz11 Gallus_gallus.GRCg6a
do
	# List all the protein coding genes with CDS annotation in the GTF
	zcat ${TranscriptomesFolder}/${Species}.103.gtf.gz | grep 'protein_coding' | awk '{if($3 == "CDS"){print $0}}' | cut -f9 | sed 's/gene_id "\([A-Z0-9]\+\)".*transcript_id "\([A-Z0-9]\+\)".*/\1\t\2/g' | sed 's/gene_id "\([A-Z0-9]\+\)".*/\1/g' | uniq | sort -u > ${ResultsFolder}/Proteomes/${Species}_ProteinCoding_AllTx_inGTF.txt
	# Extract the protein sequences properly annotated as protein coding and adapt the sequence names
	zcat ${ProteomesFolder}/${Species}.pep.all.fa.gz | awk -F ' ' '{if($1 ~ />/){if($6 == "gene_biotype:protein_coding"){print ">"$4"|"$5; valid=1;}else{valid=0}}else{if(valid==1){print $0}}}' | sed 's/gene://g' | sed 's/transcript://g' | sed 's/\.[0-9]//g' > ${ResultsFolder}/Proteomes/${Species}_ProteinCoding_AllTx.fa
	# Extract the longest tx from the ones present in the GTF file
	awk '{if(NR==FNR){a[$1]=$2;next;} if(a[$1]==$2){print $0}}' ${ResultsFolder}/Proteomes/${Species}_ProteinCoding_AllTx_inGTF.txt <(awk '{if($1 ~ />/){print name"\t"len; name=$1; len=0; next;} len=len+length($1);}END{print name"\t"len;}' ${ResultsFolder}/Proteomes/${Species}_ProteinCoding_AllTx.fa | tail -n +2 | sed 's/>//g' | sed 's/|/\t/g') | sort -k1,1 -k3,3V | awk '{if(g==$1){t=$2}else{print g"\t"t;g=$1;t=$2}}END{print g"\t"t;}' | tail -n +2 > ${ResultsFolder}/${Species}_ProteinCoding_LongestTx_inGTF.txt
	# Extract sequences from the selected longest tx
	awk '{if(NR==FNR){a[">"$1]=$2;next;} if($1 ~ />/){if(a[$1]==$2){valid=1}else{valid=0}} if(valid==1){print $1}}' ${ResultsFolder}/${Species}_ProteinCoding_LongestTx_inGTF.txt <(sed 's/|/\t/g' ${ResultsFolder}/Proteomes/${Species}_ProteinCoding_AllTx.fa) | awk -F ' ' '{if($1 ~ /^>/){print substr($1, 1, length($1)-2)}else{print $0}}' > ${ResultsFolder}/Proteomes/${Species}_ProteinCoding_LongestTx_inGTF.fa
	# Remove intermediate files
	rm ${ResultsFolder}/Proteomes/${Species}_ProteinCoding_AllTx.fa ${ResultsFolder}/Proteomes/${Species}_ProteinCoding_AllTx_inGTF.txt
done

################################################
### Species from NCBI Genome Database: 
# 		- Branchiostoma floridae
# 		- Branchiostoma belcheri
# 		- Strongylocentrotus purpuratus

for Species in Branchiostoma_belcheri.Haploidv18h27 Branchiostoma_floridae.Bfl_VNyyK
do
	SpeciesShortName=$(echo ${Species^^} | awk -F '_' '{print substr($1, 1, 1)substr($2, 1, 3)}')
	zcat ${ProteomesFolder}/${Species}.faa.gz | awk -F ' ' -v sn=${SpeciesShortName} '{if($1 ~ /^>/){print ">"sn"_"substr($1, 2, length($1))}else{print $0}}' | sed 's/\.[0-9]//g'


	
	zcat ${ProteomesFolder}/${Species}.faa.gz | awk -F ' ' '{if($1 ~ /^>/){print ">BBEL_"substr($1, 2, length($1)-3)}else{print $0}}' > ${ResultsFolder}/Proteomes/${Species}.fa
done

gzip ${ResultsFolder}/Proteomes/*.fa






