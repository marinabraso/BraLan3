#!/bin/bash



module add blast-plus/2.11.0
module add bedtools2/2.29.2
module add mafft/7.475




# Filtered GTF to tab: Blan > Results/HoxGenes/Blan_FilteredGenes.bed
zcat Results/FilteringGeneSets/GTFs/Branchiostoma_lanceolatum.BraLan3.gtf.gz | awk '{if($3 == "CDS"){print $0}}' | cut -f1,4,7,9 | sed 's/gene_id "\([A-Z0-9]\+\)".*gene_name "\([A-Za-z0-9\._-]\+\)".*/\1\t\2/g' | sort -k3 -k2V | awk '{if(g != $4){print chr"\t"start"\t"end"\t"info;chr=$1;start=$2;end=$2;info=$3"\t"$4"\t"$5; g=$4}else{end=$2}}END{print chr"\t"start"\t"end"\t"info;}' | tail -n +2 | awk '{if($2>$3){print $1"\t"$3"\t"$2"\t"$5" "$6"\t\t"$4}else{print $1"\t"$2"\t"$3"\t"$5" "$6"\t\t"$4}}' > Results/HoxGenes/Blan_FilteredGenes.bed
# Filtered GTF to tab: Bflo > Results/HoxGenes/Bflo_FilteredGenes.bed
zcat Results/FilteringGeneSets/GTFs/Branchiostoma_floridae.Bfl_VNyyK.gtf.gz | awk '{if($3 == "CDS"){print $0}}' | cut -f1,4,7,9 | tr -d \' | sed 's/\(.*\)\t\(.*\)\t\(.*\)\t.*gene_id "\([A-Z0-9]\+\)".*product "\([][A-Za-z0-9 %\/()+:>\._-]\+\)".*/\1\t\2\t\3\t\4\t\5/g' | sort -k3 -k2V | awk -F '\t' '{if(g != $4){print chr"\t"start"\t"end"\t"info;chr=$1;start=$2;end=$2;info=$3"\t"$4"\t"$5; g=$4}else{end=$2}}END{print chr"\t"start"\t"end"\t"info;}' | tail -n +2 | awk -F '\t' '{if($2>$3){print $1"\t"$3"\t"$2"\t"$5" "$6"\t\t"$4}else{print $1"\t"$2"\t"$3"\t"$5" "$6"\t\t"$4}}' > Results/HoxGenes/Bflo_FilteredGenes.txt



# Get Blan Hox genes info:
#	- Results/HoxGenes/Blan_HoxGenes.bed 
#	- Results/HoxGenes/Blan_HoxGenes_AA.fa
#	- Results/HoxGenes/Blan_HoxGenes_cDNA.fa
#	- Results/HoxGenes/Blan_HoxGenes_gDNA.fa
#	- Results/HoxGenes/Blan_HoxGenes_WholeRegion.fa
grep 'HOX' Results/HoxGenes/Blan_FilteredGenes.txt | sort -k2V > Results/HoxGenes/Blan_HoxGenes.bed
awk '{if(NR==FNR){a[$1]=1;next}if(a[FNR] || a[FNR-1]){print $0}}' <(grep -n -f <(cut -f4 Results/HoxGenes/Blan_HoxGenes.bed | cut -f1 -d' ') Results/FilteringGeneSets/Proteomes/Branchiostoma_lanceolatum.BraLan3.fa | cut -f1 -d':') Results/FilteringGeneSets/Proteomes/Branchiostoma_lanceolatum.BraLan3.fa > Results/HoxGenes/Blan_HoxGenes_AA.fa
awk '{if(NR==FNR){a[$1]=1;next}if(a[FNR] || a[FNR-1]){print $0}}' <(grep -n -f <(cut -f4 Results/HoxGenes/Blan_HoxGenes.bed | cut -f1 -d' ') Results/FilteringGeneSets/DNASequences/Branchiostoma_lanceolatum.BraLan3_DNA.fa | cut -f1 -d':') Results/FilteringGeneSets/DNASequences/Branchiostoma_lanceolatum.BraLan3_DNA.fa > Results/HoxGenes/Blan_HoxGenes_cDNA.fa
bedtools getfasta -name -fi Data/Genomes/Branchiostoma_lanceolatum.BraLan3_genome.fa -bed Results/HoxGenes/Blan_HoxGenes.bed > Results/HoxGenes/Blan_HoxGenes_gDNA.fa
bedtools getfasta -name -fi Data/Genomes/Branchiostoma_lanceolatum.BraLan3_genome.fa -bed <(cat Results/HoxGenes/Blan_HoxGenes.bed | awk '{if(NR==1){chr=$1;start=$2}else{end=$3}}END{print chr"\t"start"\t"end}') > Results/HoxGenes/Blan_HoxGenes_WholeRegion.fa



# Get Bflo Hox genes info:
#	- Results/HoxGenes/Bflo_HoxGenes.bed 
#	- Results/HoxGenes/Bflo_HoxGenes_AA.fa
#	- Results/HoxGenes/Bflo_HoxGenes_cDNA.fa
#	- Results/HoxGenes/Bflo_HoxGenes_gDNA.fa
#	- Results/HoxGenes/Bflo_HoxGenes_WholeRegion.fa
cat Results/HoxGenes/Bflo_FilteredGenes.txt | sort -k1V -k2V | awk '{if($1 == "16" && $2 > 10160000){print $0}}' | grep 'homeobox protein' > Results/HoxGenes/Bflo_HoxGenes.bed
awk '{if(NR==FNR){a[$1]=1;next}if(a[FNR] || a[FNR-1]){print $0}}' <(grep -n -f <(cut -f4 Results/HoxGenes/Bflo_HoxGenes.bed | cut -f1 -d' ') Results/FilteringGeneSets/Proteomes/Branchiostoma_floridae.Bfl_VNyyK.fa | cut -f1 -d':') Results/FilteringGeneSets/Proteomes/Branchiostoma_floridae.Bfl_VNyyK.fa > Results/HoxGenes/Bflo_HoxGenes_AA.fa
awk '{if(NR==FNR){a[$1]=1;next}if(a[FNR] || a[FNR-1]){print $0}}' <(grep -n -f <(cut -f4 Results/HoxGenes/Bflo_HoxGenes.bed | cut -f1 -d' ') Results/FilteringGeneSets/DNASequences/Branchiostoma_floridae.Bfl_VNyyK_DNA.fa | cut -f1 -d':') Results/FilteringGeneSets/DNASequences/Branchiostoma_floridae.Bfl_VNyyK_DNA.fa > Results/HoxGenes/Bflo_HoxGenes_cDNA.fa
cat Data/Genomes/Branchiostoma_floridae.Bfl_VNyyK_genome.fa | sed 's/>\([0-9]\+\) .*/>\1/g' > Results/HoxGenes/Branchiostoma_floridae.Bfl_VNyyK_genome.fa
bedtools getfasta -name -fi Results/HoxGenes/Branchiostoma_floridae.Bfl_VNyyK_genome.fa -bed Results/HoxGenes/Bflo_HoxGenes.bed > Results/HoxGenes/Bflo_HoxGenes_gDNA.fa
bedtools getfasta -name -fi Results/HoxGenes/Branchiostoma_floridae.Bfl_VNyyK_genome.fa -bed <(cat Results/HoxGenes/Bflo_HoxGenes.bed | awk '{if(NR==1){chr=$1;start=$2}else{end=$3}}END{print chr"\t"start"\t"end}') > Results/HoxGenes/Bflo_HoxGenes_WholeRegion.fa




makeblastdb -in Results/HoxGenes/Bflo_HoxGenes_cDNA.fa -dbtype nucl -parse_seqids
blastn -task blastn -db Results/HoxGenes/Bflo_HoxGenes_cDNA.fa -query Results/HoxGenes/Blan_HoxGenes_cDNA.fa -out Results/HoxGenes/Blan_HoxGenes_2Bflo_HoxGenes_blast.out -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qframe sframe evalue bitscore"









