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
### Homo sapiens, Mus musculus & Danio rerio
for Species in Homo_sapiens.GRCh38 Mus_musculus.GRCm39 Danio_rerio.GRCz11
do
	zcat ${ProteomesFolder}/${Species}.pep.all.fa.gz | awk -F ' ' '{if($1 ~ />/){if($6 == "gene_biotype:protein_coding"){print ">"$4"|"$5; valid=1;}else{valid=0}}else{if(valid==1){print $0}}}' | sed 's/gene://g' | sed 's/transcript://g' > ${ResultsFolder}/Proteomes/${Species}_ProteinCoding_AllTx.fa
	awk '{if($1 ~ />/){print name"\t"len; name=$1; len=0; next;} len=len+length($1);}END{print name"\t"len;}' ${ResultsFolder}/Proteomes/${Species}_ProteinCoding_AllTx.fa | tail -n +2 | sed 's/>//g' | sed 's/|/\t/g' | sort -k1,1 -k3,3V | awk '{if(g==$1){t=$2}else{print g"\t"t;g=$1;t=$2}}END{print g"\t"t;}' | tail -n +2 > ${ResultsFolder}/${Species}_ProteinCoding_LongestTx.txt
	awk '{if(NR==FNR){a[">"$1]=$2;next;} if($1 ~ />/){if(a[$1]==$2){valid=1}else{valid=0}} if(valid==1){print $1}}' ${ResultsFolder}/${Species}_ProteinCoding_LongestTx.txt <(sed 's/|/\t/g' ${ResultsFolder}/Proteomes/${Species}_ProteinCoding_AllTx.fa) | awk -F ' ' '{if($1 ~ /^>/){print substr($1, 1, length($1)-2)}else{print $0}}' > ${ResultsFolder}/Proteomes/${Species}_ProteinCoding_LongestTx.fa
	rm ${ResultsFolder}/Proteomes/${Species}_ProteinCoding_AllTx.fa
done

################################################
### Branchiostoma floridae & belcheri
zcat ${ProteomesFolder}/Branchiostoma_belcheri.Haploidv18h27.faa.gz | awk -F ' ' '{if($1 ~ /^>/){print ">BBEL_"substr($1, 2, length($1)-3)}else{print $0}}' > ${ResultsFolder}/Proteomes/Branchiostoma_belcheri.Haploidv18h27.fa
zcat ${ProteomesFolder}/Branchiostoma_floridae.Bfl_VNyyK.faa.gz | awk -F ' ' '{if($1 ~ /^>/){print ">BFLO_"substr($1, 2, length($1)-3)}else{print $0}}' > ${ResultsFolder}/Proteomes/Branchiostoma_floridae.Bfl_VNyyK.fa


gzip ${ResultsFolder}/Proteomes/*.fa






