#!/bin/bash



module add blast-plus/2.11.0
module add bedtools2/2.29.2


makeblastdb -in Branchiostoma_lanceolatum.BraLan3_genome.fa -dbtype nucl -parse_seqids

blastn -task blastn -db Branchiostoma_lanceolatum.BraLan3_genome.fa -query BLAG12000123_HOXA1_dna.fa -out BLAG12000123_HOXA1_2Genome.out -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qframe sframe evalue bitscore"



makeblastdb -in Branchiostoma_lanceolatum.BraLan3_proteome.fa -dbtype prot -parse_seqids

blastp -db Branchiostoma_lanceolatum.BraLan3_proteome.fa -query BLAG12000123_HOXA1_aa.fa -out BLAG12000123_HOXA1_2Proteome.out -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qframe sframe evalue bitscore"

grep -f <(cat BLAG12000123_HOXA1_2Proteome.out | head -20 | cut -f2) <(zcat ../../Data/Transcriptomes/Branchiostoma_lanceolatum.BraLan3.gtf.gz | grep 'start_codon') | cut -f1,4,9 | sed 's/gene_id "\([A-Z0-9]\+\)".*gene_name "\([A-Z0-9a-z_]\+\)".*/\1\t\2/g'



zcat ../../Data/Transcriptomes/Branchiostoma_lanceolatum.BraLan3.gtf.gz | grep '_codon' | grep '"HOX' | cut -f1,3,4,7,9 | sed 's/gene_id "\([A-Z0-9]\+\)".*gene_name "\([A-Z0-9a-z_]\+\)".*/\1\t\2/g' | awk '{if($2 == "stop_codon"){end=$3; if(start > end){end=start;start=$3}print $1"\t"start"\t"end"\t"$5":"$6"\t\t"$4}else{start=$3;}}' > HOX_Genes.bed 
bedtools getfasta -name -fi Branchiostoma_lanceolatum.BraLan3_genome.fa -bed HOX_Genes.bed > HOX_Genes_DNA.fa

awk -F '\t' '{if(NR==FNR){a[$1]=1;next;}if(a[FNR] || a[FNR-1]){print $0}}' <(grep -nf <(cat HOX_Genes.bed | cut -f4 | cut -f1 -d':') ../FilteringGeneSets/DNASequences/Branchiostoma_lanceolatum.BraLan3_DNA.fa | cut -f1 -d":") ../FilteringGeneSets/DNASequences/Branchiostoma_lanceolatum.BraLan3_DNA.fa > HOX_Genes_cDNA.fa
awk -F '\t' '{if(NR==FNR){a[$1]=1;next;}if(a[FNR] || a[FNR-1]){print $0}}' <(grep -nf <(cat HOX_Genes.bed | cut -f4 | cut -f1 -d':') ../FilteringGeneSets/Proteomes/Branchiostoma_lanceolatum.BraLan3.fa | cut -f1 -d":") ../FilteringGeneSets/Proteomes/Branchiostoma_lanceolatum.BraLan3.fa > HOX_Genes_AA.fa 
awk '{if(NR==FNR){a[">"$4]=">"$4":"$5"::"$1":"$2"-"$3; next} if(a[$1]){print a[$1]}else{print $0}}' <(cat HOX_Genes.bed | sed 's/:/\t/g') HOX_Genes_AA.fa > HOX_Genes_AA_names.fa
awk '{if(NR==FNR){a[">"$4]=">"$4":"$5"::"$1":"$2"-"$3; next} if(a[$1]){print a[$1]}else{print $0}}' <(cat HOX_Genes.bed | sed 's/:/\t/g') HOX_Genes_cDNA.fa > HOX_Genes_cDNA_names.fa



blastp -db Branchiostoma_lanceolatum.BraLan3_proteome.fa -query Human_HOXA1.fa -out Human_HOXA1_2Blan3Proteome.out -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qframe sframe evalue bitscore"


zcat ../../Data/Transcriptomes/Branchiostoma_lanceolatum.BraLan3.gtf.gz | grep 'start_codon' | grep '"HOX' | cut -f1,4,9 | sed 's/gene_id "\([A-Z0-9]\+\)".*gene_name "\([A-Z0-9a-z_]\+\)".*/\1\t\2/g'


makeblastdb -in Branchiostoma_floridae.Bfl_VNyyK.fa -dbtype prot -parse_seqids
blastp -db Branchiostoma_floridae.Bfl_VNyyK.fa -query BLAG12000123_HOXA1_aa.fa -out BLAG12000123_HOXA1_2BfloProteome.out -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qframe sframe evalue bitscore"

blastp -db Branchiostoma_floridae.Bfl_VNyyK.fa -query BFLOLOC118403096_HOX3.fa -out BFLOLOC118403096_HOX3_2BfloProteome.out -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qframe sframe evalue bitscore"
blastp -db Branchiostoma_floridae.Bfl_VNyyK.fa -query HOX_Genes_AA.fa -out HOX_Genes_AA_2BfloProteome.out -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qframe sframe evalue bitscore"




zcat ../FilteringGeneSets/GTFs/Branchiostoma_floridae.Bfl_VNyyK.gtf.gz | awk '{if($3 == "CDS"){print $0}}' | cut -f1,4,9 | sed 's/ID .*gene_id "\(.*\)"; product "\(.*\)"; protein_id .*/\1\t\2/g' | sort -k3,3 -k2,2V | uniq | awk '{if(g != $3){g=$3;print $0}}' | awk '{if($1==16 && $2 > 10040000 && $2 < 14693000){print $0}}' | grep 'homeobox protein ' | sort -k2,2V > Bflo_HOX_genes.txt
awk '{if(NR==FNR){a[$1]=1;next;} if(a[FNR] || a[FNR-1]){print $0}}' <(grep -n -f <(cut -f3 Bflo_HOX_genes.txt) ../FilteringGeneSets/Proteomes/Branchiostoma_floridae.Bfl_VNyyK.fa | cut -f1 -d':') ../FilteringGeneSets/Proteomes/Branchiostoma_floridae.Bfl_VNyyK.fa > Bflo_HOX_genes_AA.fa
awk '{if(NR==FNR){a[$1]=1;next;} if(a[FNR] || a[FNR-1]){print $0}}' <(grep -n -f <(cut -f3 Bflo_HOX_genes.txt) ../FilteringGeneSets/DNASequences/Branchiostoma_floridae.Bfl_VNyyK_DNA.fa | cut -f1 -d':') ../FilteringGeneSets/DNASequences/Branchiostoma_floridae.Bfl_VNyyK_DNA.fa > Bflo_HOX_genes_cDNA.fa


makeblastdb -in Braflo.fa -dbtype prot -parse_seqids
blastp -db Braflo.fa -query HOX_Genes_AA.fa -out HOX_Genes_AA_2BfloFerdiProteome.out -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qframe sframe evalue bitscore"

blastp -db Branchiostoma_floridae.Bfl_VNyyK.fa -query HOX_Genes_AA_2BfloFerdiProteome.fa -out BfloFerdi_HOX_Genes_AA_2BfloProteome.out -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qframe sframe evalue bitscore"







