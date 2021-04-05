#!/bin/bash

# Scripts
bedtools="../Software/bedtools.static.binary"
ExtractcDNAPerlScript="Scripts/FilteringGeneSets/ExtractDNAseqFromGenome.pl"

# Files & parameters
ProteomesFolder="Data/Proteomes"
GTFFolder="Data/Transcriptomes"
GenomesFolder="Data/Genomes"
ResultsFolder="Results/FilteringGeneSets"
mkdir -p ${ResultsFolder}/FilteredProteomes
mkdir -p ${ResultsFolder}/FilteredGTFs
mkdir -p ${ResultsFolder}/FilteredDNASequences
mkdir -p ${ResultsFolder}/BacktranslatedSequences
mkdir -p ${ResultsFolder}/Checking_DNA_AA_sequences
mkdir -p ${ResultsFolder}/Proteomes
mkdir -p ${ResultsFolder}/GTFs
mkdir -p ${ResultsFolder}/DNASequences

# This script:
# - Homogenizes format of the downloaded proteomes
# - Filters for protein coding genes and selects one transcript (longest) per gene
# - Extracts the cDNA sequence of the filtered proteins
# - Checks correspondande of DNA - AA sequences
# - Filters again for the sequences that do correspond 


################################################
### Homogenize format
echo "### HOMOGENIZING FORMAT, EXTRACTING GENE SET..."
################################################
### Branchiostoma lanceolatum
echo "	Branchiostoma_lanceolatum.BraLan3"
if [[ ! -s ${ResultsFolder}/Branchiostoma_lanceolatum.BraLan3_gene2tx_tmp.txt ]]; then
	#	Filter for protein coding, "strong" evidence in any of the fields & print the longest transcript per gene
	zcat ${GTFFolder}/Branchiostoma_lanceolatum.BraLan3.gtf.gz | grep 'strong' | grep 'protein_coding' | awk '{if($3 == "CDS"){print $0}}' | sed 's/"; /"\t/g' | cut -f4,5,9,10 | sed 's/\t[a-z_]\+ "/\t/g' | sed 's/"//g' | sort -k3,3 -k4,4 -k1,1V | awk '{if(n != $3"\t"$4){print n"\t"len; n=$3"\t"$4; len=0}else{len=len+$2-$1}}END{print n"\t"len;}' | tail -n +2 | awk '{if(g != $1){print g"\t"t; g=$1; t=$2;maxlen=$3}else{if($3>maxlen){t=$2;maxlen=$3}}}END{print g"\t"t;}' | tail -n +2 > ${ResultsFolder}/Branchiostoma_lanceolatum.BraLan3_gene2tx_tmp.txt
fi
if [[ ! -s ${ResultsFolder}/FilteredProteomes/Branchiostoma_lanceolatum.BraLan3.fa ]]; then
	#	Extract protein sequences of selected trancripts 
	awk '{if(NR==FNR){a[">"$1]=$2;next} if($1 ~ />/){if(a[$1]){print $1;valid=1}else{valid=0}}else{if(valid==1){print $0}}}' Results/FilteringGeneSets/Branchiostoma_lanceolatum.BraLan3_gene2tx_tmp.txt <(zcat Data/Proteomes/Branchiostoma_lanceolatum.BraLan3.fa.gz) | sed 's/\.//g' > ${ResultsFolder}/FilteredProteomes/Branchiostoma_lanceolatum.BraLan3.fa
fi
if [[ ! -s ${ResultsFolder}/Branchiostoma_lanceolatum.BraLan3_gene2tx.txt ]]; then
	#	Check that only present sequences are in the gene2tx file
	awk '{if(NR==FNR){a[$1]=1;next}if(a[$1]){print $0}}' <(cat Results/FilteringGeneSets/FilteredProteomes/Branchiostoma_lanceolatum.BraLan3.fa | grep '>' | sed 's/>//g') Results/FilteringGeneSets/Branchiostoma_lanceolatum.BraLan3_gene2tx_tmp.txt > Results/FilteringGeneSets/Branchiostoma_lanceolatum.BraLan3_gene2tx.txt
fi
if [[ ! -s ${ResultsFolder}/FilteredGTFs/Branchiostoma_lanceolatum.BraLan3.gtf ]]; then
	#	Extract GTF
	awk -F'\t' '{if(NR==FNR){a[$1]=$2; next} n=split($9,b,"\x22"); for(i=1;i<n;i++){for(j=i+1;j<=n;j++){if((a[b[i]]==b[j] && a[b[i]] ~ sn) || (a[b[j]]==b[i] && a[b[j]] ~ sn)){print $0}}}}' <(cat ${ResultsFolder}/Branchiostoma_lanceolatum.BraLan3_gene2tx.txt) <(zcat ${GTFFolder}/Branchiostoma_lanceolatum.BraLan3.gtf.gz) | sort -u > ${ResultsFolder}/FilteredGTFs/Branchiostoma_lanceolatum.BraLan3.gtf
fi
	


################################################
### Species from ENSEMBL 103: 
#		- Homo sapiens
# 		- Mus musculus
# 		- Gallus gallus
# 		- Danio rerio
for Species in Homo_sapiens.GRCh38 Mus_musculus.GRCm39 Danio_rerio.GRCz11 Gallus_gallus.GRCg6a
do
	echo "	"${Species}
	if [[ ! -s ${ResultsFolder}/${Species}_gene2tx_tmp.txt ]]; then
		#	Filter for protein coding, "strong" evidence in any of the fields & print the longest transcript per gene
		zcat ${GTFFolder}/${Species}.103.gtf.gz | grep -v '^#' | grep 'gene_biotype "protein_coding"' | awk '{if($3 == "CDS"){print $0}}' | sed 's/"; /"\t/g' | cut -f4,5,9,11 | sed 's/\t[a-z_]\+ "/\t/g' | sed 's/"//g' | sort -k3,3 -k4,4 -k1,1V | awk '{if(n != $3"\t"$4){print n"\t"len; n=$3"\t"$4; len=0}else{len=len+$2-$1}}END{print n"\t"len;}' | tail -n +2 | awk '{if(g != $1){print g"\t"t; g=$1; t=$2;maxlen=$3}else{if($3>maxlen){t=$2;maxlen=$3}}}END{print g"\t"t;}' | tail -n +2 > ${ResultsFolder}/${Species}_gene2tx_tmp.txt
	fi
	if [[ ! -s ${ResultsFolder}/FilteredProteomes/${Species}.fa ]]; then
		#	Extract protein sequences of selected trancripts 
		awk '{if(NR==FNR){a[">"$1"\t"$2]=1;next} if($1 ~ />/){if(a[$1"\t"$2]){print $1;valid=1}else{valid=0}}else{if(valid==1){print $0}}}' Results/FilteringGeneSets/${Species}_gene2tx_tmp.txt <(zcat Data/Proteomes/${Species}.pep.all.fa.gz | awk -F' ' '{if($1 ~ />/){print seq; seq=""; print ">"$4"\t"$5}else{seq=seq$0}}END{print seq}' | tail -n +2 |  sed 's/gene:\([A-Z0-9]\+\)\.[0-9]\+\ttranscript:\([A-Z0-9]\+\)\.[0-9]\+/\1\t\2/g') > ${ResultsFolder}/FilteredProteomes/${Species}.fa
	fi
	if [[ ! -s ${ResultsFolder}/${Species}_gene2tx.txt ]]; then
		#	Check that only present sequences are in the gene2tx file
		awk '{if(NR==FNR){a[$1]=1;next}if(a[$1]){print $0}}' <(cat Results/FilteringGeneSets/FilteredProteomes/${Species}.fa | grep '>' | sed 's/>//g') Results/FilteringGeneSets/${Species}_gene2tx_tmp.txt > Results/FilteringGeneSets/${Species}_gene2tx.txt
	fi
	if [[ ! -s ${ResultsFolder}/FilteredGTFs/${Species}.gtf ]]; then
		#	Extract GTF
		awk -F'\t' '{if(NR==FNR){a[$1]=$2; next} n=split($9,b,"\x22"); for(i=1;i<n;i++){for(j=i+1;j<=n;j++){if(a[b[i]]==b[j] || a[b[j]]==b[i]){print $0}}}}' Results/FilteringGeneSets/${Species}_gene2tx.txt <(zcat ${GTFFolder}/${Species}.103.gtf.gz | grep -v '^#') | sort -u > ${ResultsFolder}/FilteredGTFs/${Species}.gtf
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
	echo "	"${Species}
	SpeciesShortName=$(echo ${Species^^} | awk -F '_' '{print substr($1, 1, 1)substr($2, 1, 3)}')
	if [[ ! -s ${ResultsFolder}/${Species}_gene2tx_tmp.txt ]]; then
		#	Filter for protein coding, "strong" evidence in any of the fields & print the longest transcript per gene
		zcat ${GTFFolder}/${Species}.gff.gz | grep -v '^#' | awk '{if($3=="CDS"){print $0}}' | sed 's/Parent=rna-/transcript_id=/g'  | sed 's/Parent=gene-/transcript_id=/g' | sed 's/ID.*;transcript_id=\([0-9a-zA-Z_\/\.\-]\+\);.*gene=\([0-9a-zA-Z_\/\.\-]\+\);.*protein_id=\([0-9a-zA-Z_\/\.\-]\+\).*/\2\t\1\t\3/g' | sed 's/\.[0-9]\+//g' | cut -f4,5,9,10,11 | awk -v sn=${SpeciesShortName} '{print $1"\t"$2"\t"sn""$3"\t"sn""$4";"sn""$5}' | sort -k3,3 -k4,4 -k1,1V | awk '{if(n != $3"\t"$4){print n"\t"len; n=$3"\t"$4; len=0}else{len=len+$2-$1}}END{print n"\t"len;}' | tail -n +2 | awk '{if(g != $1){print g"\t"t; g=$1; t=$2;maxlen=$3}else{if($3>maxlen){t=$2;maxlen=$3}}}END{print g"\t"t;}' | tail -n +2 > ${ResultsFolder}/${Species}_gene2tx_tmp.txt
	fi
	if [[ ! -s ${ResultsFolder}/FilteredProteomes/${Species}.fa ]]; then
		#	Extract protein sequences of selected trancripts 
		awk '{if(NR==FNR){a[">"$3]=$1;next} if($1 ~ />/){if(a[$1]){print ">"a[$1];valid=1}else{valid=0}}else{if(valid==1){print $0}}}' <(cat Results/FilteringGeneSets/${Species}_gene2tx_tmp.txt | sed 's/;/\t/g') <(zcat Data/Proteomes/${Species}.faa.gz | awk -F' ' '{if($1 ~ />/){print seq; seq=""; print $1}else{seq=seq$0}}END{print seq}' | tail -n +2 | sed 's/\.[0-9]\+//g' | sed "s/>/>${SpeciesShortName}/g")  > ${ResultsFolder}/FilteredProteomes/${Species}.fa
	fi
	if [[ ! -s ${ResultsFolder}/${Species}_gene2tx.txt ]]; then
		#	Check that only present sequences are in the gene2tx file
		awk '{if(NR==FNR){a[$1]=1;next}if(a[$1]){print $1"\t"$2}}' <(cat Results/FilteringGeneSets/FilteredProteomes/${Species}.fa | grep '>' | sed 's/>//g') Results/FilteringGeneSets/${Species}_gene2tx_tmp.txt > ${ResultsFolder}/${Species}_gene2tx.txt
	fi
	if [[ ! -s ${ResultsFolder}/FilteredGTFs/${Species}.gtf ]]; then
		#	Convert format, Extract GTF of the selected transcripts
		awk -v sn=${SpeciesShortName} -F'\t' '{if(NR==FNR){a[$1]=$2; next} n=split($9,b,"\x22"); for(i=1;i<n;i++){for(j=i+1;j<=n;j++){if((a[b[i]]==b[j] && a[b[i]] ~ sn) || (a[b[j]]==b[i] && a[b[j]] ~ sn)){print $0}}}}' <(cat ${ResultsFolder}/${Species}_gene2tx.txt | sed 's/;/\t/g') <(zcat ${GTFFolder}/${Species}.gff.gz | grep -v '^#' | sed 's/Parent=rna-/transcript_id=/g'  | sed 's/Parent=gene-/transcript_id=/g' | sed 's/=/ "/g' | sed 's/;/"; /g' | sed 's/$/";/g' | sed 's/ transcript_id \"[0-9a-zA-Z_\/\.\-]\+\";\(.* transcript_id \"[0-9a-zA-Z_\/\.\-]\+\";\)/\1/g' | sed 's/\.[0-9]\+//g' | sed "s/transcript_id \"/transcript_id \"${SpeciesShortName}/g" | sed "s/gene \"/gene \"${SpeciesShortName}/g" | sed "s/protein_id \"/protein_id \"${SpeciesShortName}/g" | sed 's/gene "/gene_id "/g') | sort -u > ${ResultsFolder}/FilteredGTFs/${Species}.gtf
	fi
done


################################################
###
### EXTRACTING DNA SEQUENCES
###
################################################
echo "### EXTRACTING DNA SEQUENCES"
for Species in Branchiostoma_lanceolatum.BraLan3 Homo_sapiens.GRCh38 Mus_musculus.GRCm39 Danio_rerio.GRCz11 Gallus_gallus.GRCg6a Branchiostoma_belcheri.Haploidv18h27 Branchiostoma_floridae.Bfl_VNyyK Strongylocentrotus_purpuratus.Spur5.0 Asterias_rubens.eAstRub1.3 Saccoglossus_kowalevskii.Skow1.1
#for Species in Branchiostoma_lanceolatum.BraLan3
do
	echo "	"${Species}
	GenomeFile=$(ls ${GenomesFolder}/${Species}*.fa)
	if [[ ! -s ${ResultsFolder}/${Species}_Extract_cDNA_tmp.txt  ]]; then
		cat ${ResultsFolder}/FilteredGTFs/${Species}.gtf | awk '{if($3 == "CDS"){print $0}}' | sed 's/\(.*\)\t\(.*\)\t\(.*\)\t\(.*\)\t\(.*\)\t\(.*\)\t\(.*\)\t\(.*\)\t.*gene_id "\([A-Z0-9_\.\-]\+\)".*/\1\t\2\t\3\t\4\t\5\t\6\t\7\t\8\t\9/g' | awk '{print $9"\t"$1"\t"$4"\t"$5"\t"$7}' | sort -k1,1 -k2,2V -k3,3V | uniq > ${ResultsFolder}/${Species}_Extract_cDNA_tmp.txt 
	fi
	if [[ ! -s ${ResultsFolder}/FilteredDNASequences/${Species}_DNA.fa ]]; then
		perl ${ExtractcDNAPerlScript} --pathCoordinates=${ResultsFolder}/${Species}_Extract_cDNA_tmp.txt  --pathGenomeSequence=${GenomeFile} --pathOutput=${ResultsFolder}/FilteredDNASequences/${Species}_DNA.fa > tmp/Perl_ExtractDNA.out
	fi
done

################################################
###
### CHECKING DNA - AA SEQUENCES CORRESPONDENCE
###
################################################
echo "### CHECKING DNA - AA SEQUENCES CORRESPONDENCE"
for Species in Branchiostoma_lanceolatum.BraLan3 Homo_sapiens.GRCh38 Mus_musculus.GRCm39 Danio_rerio.GRCz11 Gallus_gallus.GRCg6a Branchiostoma_belcheri.Haploidv18h27 Branchiostoma_floridae.Bfl_VNyyK Strongylocentrotus_purpuratus.Spur5.0 Asterias_rubens.eAstRub1.3 Saccoglossus_kowalevskii.Skow1.1
#for Species in Branchiostoma_lanceolatum.BraLan3
do
	echo "	"${Species}
	# Nuclear Genetic Code Backtranslation
	if [[ ! -s ${ResultsFolder}/BacktranslatedSequences/${Species}_bt.fa ]]; then
		echo "		Nuclear backtranslation"
		cat ${ResultsFolder}/FilteredDNASequences/${Species}_DNA.fa | tr '[:lower:]' '[:upper:]' | awk 'BEGIN{
			a["TTT"]="F"; 
			a["TTC"]="F"; 
			a["TTA"]="L"; 
			a["TTG"]="L"; 
			a["TCT"]="S"; 
			a["TCC"]="S"; 
			a["TCA"]="S"; 
			a["TCG"]="S"; 
			a["TAT"]="Y"; 
			a["TAC"]="Y"; 
			a["TAA"]="1"; 
			a["TAG"]="1"; 
			a["TGT"]="C"; 
			a["TGC"]="C"; 
			a["TGA"]="1"; 
			a["TGG"]="W"; 
			a["CTT"]="L"; 
			a["CTC"]="L"; 
			a["CTA"]="L"; 
			a["CTG"]="L"; 
			a["CCT"]="P"; 
			a["CCC"]="P"; 
			a["CCA"]="P"; 
			a["CCG"]="P"; 
			a["CAT"]="H"; 
			a["CAC"]="H"; 
			a["CAA"]="Q"; 
			a["CAG"]="Q"; 
			a["CGT"]="R"; 
			a["CGC"]="R"; 
			a["CGA"]="R"; 
			a["CGG"]="R"; 
			a["ATT"]="I"; 
			a["ATC"]="I"; 
			a["ATA"]="I"; 
			a["ATG"]="M"; 
			a["ACT"]="T"; 
			a["ACC"]="T"; 
			a["ACA"]="T"; 
			a["ACG"]="T"; 
			a["AAT"]="N"; 
			a["AAC"]="N"; 
			a["AAA"]="K"; 
			a["AAG"]="K"; 
			a["AGT"]="S"; 
			a["AGC"]="S"; 
			a["AGA"]="R"; 
			a["AGG"]="R"; 
			a["GTT"]="V"; 
			a["GTC"]="V"; 
			a["GTA"]="V"; 
			a["GTG"]="V"; 
			a["GCT"]="A"; 
			a["GCC"]="A"; 
			a["GCA"]="A"; 
			a["GCG"]="A"; 
			a["GAT"]="D"; 
			a["GAC"]="D"; 
			a["GAA"]="E"; 
			a["GAG"]="E"; 
			a["GGT"]="G"; 
			a["GGC"]="G"; 
			a["GGA"]="G"; 
			a["GGG"]="G"; 
		}{if($1 ~ />/){
			print $0; 
			btseq=""
		}else{
			for(i=1;i<=length($0);i=i+3){
				btseq=btseq""a[substr($0, i, 3)]
			}
			if(substr(btseq, length(btseq), 1)== "1"){
				btseq=substr(btseq, 1, length(btseq)-1)
			}
			print btseq
		}}' > ${ResultsFolder}/BacktranslatedSequences/${Species}_bt.fa
	fi
	# Mitochondrial Genetic Code Backtranslation
	if [[ ! -s ${ResultsFolder}/BacktranslatedSequences/${Species}_mbt.fa ]]; then
		echo "		Mitochondrial backtranslation"
		cat ${ResultsFolder}/FilteredDNASequences/${Species}_DNA.fa | tr '[:lower:]' '[:upper:]' | awk 'BEGIN{
			a["TTT"]="F"; 
			a["TTC"]="F"; 
			a["TTA"]="L"; 
			a["TTG"]="L"; 
			a["TCT"]="S"; 
			a["TCC"]="S"; 
			a["TCA"]="S"; 
			a["TCG"]="S"; 
			a["TAT"]="Y"; 
			a["TAC"]="Y"; 
			a["TAA"]="1"; 
			a["TAG"]="1"; 
			a["TGT"]="C"; 
			a["TGC"]="C"; 
			a["TGA"]="1"; 
			a["TGG"]="W"; 
			a["CTT"]="L"; 
			a["CTC"]="L"; 
			a["CTA"]="L"; 
			a["CTG"]="L"; 
			a["CCT"]="P"; 
			a["CCC"]="P"; 
			a["CCA"]="P"; 
			a["CCG"]="P"; 
			a["CAT"]="H"; 
			a["CAC"]="H"; 
			a["CAA"]="Q"; 
			a["CAG"]="Q"; 
			a["CGT"]="R"; 
			a["CGC"]="R"; 
			a["CGA"]="R"; 
			a["CGG"]="R"; 
			a["ATT"]="I"; 
			a["ATC"]="I"; 
			a["ATA"]="I"; 
			a["ATG"]="M"; 
			a["ACT"]="T"; 
			a["ACC"]="T"; 
			a["ACA"]="T"; 
			a["ACG"]="T"; 
			a["AAT"]="N"; 
			a["AAC"]="N"; 
			a["AAA"]="K"; 
			a["AAG"]="K"; 
			a["AGT"]="S"; 
			a["AGC"]="S"; 
			a["AGA"]="R"; 
			a["AGG"]="R"; 
			a["GTT"]="V"; 
			a["GTC"]="V"; 
			a["GTA"]="V"; 
			a["GTG"]="V"; 
			a["GCT"]="A"; 
			a["GCC"]="A"; 
			a["GCA"]="A"; 
			a["GCG"]="A"; 
			a["GAT"]="D"; 
			a["GAC"]="D"; 
			a["GAA"]="E"; 
			a["GAG"]="E"; 
			a["GGT"]="G"; 
			a["GGC"]="G"; 
			a["GGA"]="G"; 
			a["GGG"]="G"; 
			for(c in a){b[c]=a[c]};
			b["AGA"]="1"; 
			b["AGG"]="1"; 
			b["ATA"]="M"; 
			b["TGA"]="W"; 
		}{if($1 ~ />/){
			print $0; 
			btseq=""
		}else{
			for(i=1;i<=length($0);i=i+3){
				btseq=btseq""a[substr($0, i, 3)]
			}
			if(substr(btseq, length(btseq), 1)== "1"){
				btseq=substr(btseq, 1, length(btseq)-1)
			}
			print btseq
		}}' > ${ResultsFolder}/BacktranslatedSequences/${Species}_mbt.fa
	fi
	# Comparing the bt sequences with the actual AA sequence
	if [[ ! -s ${ResultsFolder}/Checking_DNA_AA_sequences/${Species}_checking_btDNA_AA.txt ]]; then
		echo "		Comparison"
		awk 'BEGIN{file=0}{if(FNR==1){file++} 
			# Reading File with nuclear backtranslations
			if(file==1){
				if($1 ~ />/){
					g=$1;next;
				}
				a[g]=$1;
				next;
			} 
			# Reading file with mitochondrial backtranslations
			if(file==2){
				if($1 ~ />/){
					g=$1;
					next;
				} 
				b[g]=$1;
				next;
			} 
			# Reading original proteomes
			if($1 ~ />/){
				g=$1;
				next;
			} 
			if(a[g]==$1){
				print g"\tIdentical\tNuclear"
			}else{
				if(b[g]==$1){
					print g"\tIdentical\tMitochondrial"
				}else{
					if(length(b[g])==length($1)){
						str = "";
						mismatch = 0;
						for(i=1;i<=length(b[g]);i++){
							if(substr(b[g],i,1) != substr($1,i,1)){
								str=str" "substr(b[g],i,1)"|"substr($1,i,1);
								mismatch++;
							}
						}
						print g"\tError\tSame length\t"length(b[g])"\t"mismatch"\t"str"\t"$1"\t"a[g]"\t"b[g]
					}else{
						print g"\tError\tDifferent length\t"length(b[g])"\t"length($1)"\t"$1"\t"a[g]"\t"b[g]						
					}
				}
			}}' ${ResultsFolder}/BacktranslatedSequences/${Species}_bt.fa ${ResultsFolder}/BacktranslatedSequences/${Species}_mbt.fa ${ResultsFolder}/FilteredProteomes/${Species}.fa > ${ResultsFolder}/Checking_DNA_AA_sequences/${Species}_checking_btDNA_AA.txt
	fi
done








################################################
###
### FINAL SET OF GENES
###
################################################
echo "### FINAL SET OF GENES"
for Species in Branchiostoma_lanceolatum.BraLan3 Homo_sapiens.GRCh38 Mus_musculus.GRCm39 Danio_rerio.GRCz11 Gallus_gallus.GRCg6a Branchiostoma_belcheri.Haploidv18h27 Branchiostoma_floridae.Bfl_VNyyK Strongylocentrotus_purpuratus.Spur5.0 Asterias_rubens.eAstRub1.3 Saccoglossus_kowalevskii.Skow1.1
#for Species in Branchiostoma_lanceolatum.BraLan3
do
	echo "	"${Species}
	# List of selected genes
	cat ${ResultsFolder}/Checking_DNA_AA_sequences/${Species}_checking_btDNA_AA.txt | grep 'Same' | awk -F'\t' '{if($5/$4 <= 0.1){print $1}}' | sed 's/>//g' >  ${ResultsFolder}/${Species}_FinalSet.list
	cat ${ResultsFolder}/Checking_DNA_AA_sequences/${Species}_checking_btDNA_AA.txt | grep 'Identical' | cut -f1 | sed 's/>//g' >> ${ResultsFolder}/${Species}_FinalSet.list

	if [[ ! -s ${ResultsFolder}/Proteomes/${Species}.fa ]]; then
		#	Extract protein sequences of selected genes
		awk '{if(NR==FNR){a[">"$1]=1;next} if($1 ~ />/){if(a[$1]){print $1;valid=1}else{valid=0}}else{if(valid==1){print $0}}}' ${ResultsFolder}/${Species}_FinalSet.list ${ResultsFolder}/FilteredProteomes/${Species}.fa >  ${ResultsFolder}/Proteomes/${Species}.fa
	fi
	if [[ ! -s ${ResultsFolder}/DNASequences/${Species}_DNA.fa ]]; then
		#	Extract DNA sequences of selected genes
		awk '{if(NR==FNR){a[">"$1]=1;next} if($1 ~ />/){if(a[$1]){print $1;valid=1}else{valid=0}}else{if(valid==1){print $0}}}' ${ResultsFolder}/${Species}_FinalSet.list ${ResultsFolder}/FilteredDNASequences/${Species}_DNA.fa >  ${ResultsFolder}/DNASequences/${Species}_DNA.fa
	fi
	if [[ ! -s ${ResultsFolder}/GTFs/${Species}.fa ]]; then
		#	Extract GTF of selected genes
		grep -w -f ${ResultsFolder}/${Species}_FinalSet.list ${ResultsFolder}/FilteredGTFs/${Species}.gtf > ${ResultsFolder}/GTFs/${Species}.gtf
	fi
done










