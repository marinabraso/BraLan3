#!/bin/bash

# Scripts
bedtools="../Software/bedtools.static.binary"

# Files & parameters
ProteomesFolder="Data/Proteomes"
TranscriptomesFolder="Data/Transcriptomes"
GenomesFolder="Data/Genomes"
ResultsFolder="Results/FindOrthologs"
mkdir -p ${ResultsFolder}/FilteredProteomes
mkdir -p ${ResultsFolder}/FilteredGTFs
mkdir -p ${ResultsFolder}/FilteredDNASequences
mkdir -p ${ResultsFolder}/Checking_DNA_AA_sequences
mkdir -p ${ResultsFolder}/CheckedProteomes
mkdir -p ${ResultsFolder}/CheckedGTFs
mkdir -p ${ResultsFolder}/CheckedDNASequences

# This script:
# - Filters the downloaded proteomes
# - Extracts the cDNA sequence of the filtered proteins
# - Checks correspondande of DNA - AA sequences
# - Filters again for the sequences that do correspond 

################################################
###
### FILTERING PROTEOMES
###  - Being annotated as protein coding
###  - Having CDS annotations in the corresponding GTF/GFF file
###  - Select one transcript per gene (longest)
###
################################################
echo "### FILTERING PROTEOMES"
################################################
### Branchiostoma lanceolatum
echo "Branchiostoma_lanceolatum.BraLan3"
if [[ $(ls ${ResultsFolder}/FilteredProteomes/Branchiostoma_lanceolatum.BraLan3_ProteinCoding_StrongEvidence.fa* 2> ~/null | wc -l) -lt 1 ]] || [[ $(ls ${ResultsFolder}/FilteredGTFs/Branchiostoma_lanceolatum.BraLan3_ProteinCoding_StrongEvidence.gtf* 2> ~/null | wc -l) -lt 1 ]]
then
	echo "Filtering proteome"
	# List of protein coding genes with "strong" evidence in any of the testedevidencies
	#zcat ${TranscriptomesFolder}/Branchiostoma_lanceolatum.BraLan3.gtf.gz | awk '{if($3 ~ /CDS/){print $0}}' | grep 'protein_coding' | cut -f9 | sed 's/gene_id \"\([A-Z0-9]\+\)\"; transcript_id.*bgee_evidence \"\([a-z]\+\)\"; uniprot_evidence \"\([a-z]\+\)\"; bflo_evidence \"\([a-z]\+\)\"; bbel_evidence \"\([a-z]\+\)\";/\1\t\2\t\3\t\4\t\5/g' |  sort | uniq | grep 'strong' | cut -f1 > ${ResultsFolder}/BraLan3_ProteinCoding_StrongEvidence.txt
	## Select the fasta sequences of those in the list & substitute . by * in fasta sequences (?)
	#awk '{if(NR==FNR){a[">"$1]=1;next;} if(a[$1]==1){valid=1;}else{if($1 ~ />/){valid=0;}} if(valid==1){print $0}}' ${ResultsFolder}/BraLan3_ProteinCoding_StrongEvidence.txt <(zcat ${ProteomesFolder}/Branchiostoma_lanceolatum.BraLan3.fa.gz) | sed 's/\./*/g' > ${ResultsFolder}/FilteredProteomes/Branchiostoma_lanceolatum.BraLan3_ProteinCoding_StrongEvidence.fa
	## Extract GTF info of the selected genes
	#awk '{if(NR==FNR){a[$1]=1;next;} if(a[$9]==1){print $0}}' <(cat ${ResultsFolder}/FilteredProteomes/Branchiostoma_lanceolatum.BraLan3_ProteinCoding_StrongEvidence.fa | grep '>' | sed 's/>//g') <(zcat ${TranscriptomesFolder}/Branchiostoma_lanceolatum.BraLan3.gtf.gz | awk '{if($3 == "CDS"){print $0}}' | sed 's/gene_id "\([A-Z0-9]\+\)".*/\1/g') > ${ResultsFolder}/FilteredGTFs/Branchiostoma_lanceolatum.BraLan3_ProteinCoding_StrongEvidence.gtf
fi

################################################
### Species from ENSEMBL 103: 
#		- Homo sapiens
# 		- Mus musculus
# 		- Gallus gallus
# 		- Danio rerio
#for Species in Homo_sapiens.GRCh38 Mus_musculus.GRCm39 Danio_rerio.GRCz11 Gallus_gallus.GRCg6a
for Species in Homo_sapiens.GRCh38
do
	echo ${Species}
	if [[ $(ls ${ResultsFolder}/FilteredProteomes/${Species}_ProteinCoding_LongestTx.fa* 2> ~/null | wc -l) -lt 1 ]] ##|| [[ $(ls ${ResultsFolder}/FilteredGTFs/${Species}_ProteinCoding_LongestTx_inGTF.gtf* 2> ~/null | wc -l) -lt 1 ]]
	then
		echo "Filtering proteome"
		# Extract GeneID and TrascriptID of protein coding genes with CDS entries in the GTF
		zcat ${TranscriptomesFolder}/${Species}.103.gtf.gz | awk '{if($3 ~ /CDS/){print $0}}' | grep 'protein_coding' | sed 's/.*gene_id "\([A-Z0-9]\+\)".*transcript_id "\([A-Z0-9]\+\)".*/\1\t\2/g' | sort -u > tmp/GT_gtf.txt
		# Extract GeneID and TrascriptID of protein coding genes in the proteome (fasta)
		zcat ${ProteomesFolder}/${Species}.pep.all.fa.gz | grep '>' | grep 'gene_biotype:protein_coding' | sed 's/.*gene:\([A-Z0-9]\+\).*transcript:\([A-Z0-9]\+\).*/\1\t\2/g' | sort -u > tmp/GT_fa.txt 
		# Intersection of the two above > common list
		awk 'NR==FNR { lines[$0]=1; next } $0 in lines' tmp/GT_gtf.txt tmp/GT_fa.txt > tmp/GT_comm.txt
		# Calculate the length of each protein entry
		zcat ${ProteomesFolder}/${Species}.pep.all.fa.gz | sed 's/.*gene:\([A-Z0-9]\+\).*transcript:\([A-Z0-9]\+\).*/>\1\t\2/g' | awk -F' ' '{if($1 ~ />/){print id"\t"len; id=$0; len=0}else{len=len+length($0)}}END{print id"\t"len;}' | tail -n +2 | sed 's/>//g' > tmp/GTL_fa.txt
		# Intersect common list with the length file and extract the lognest transcript per gene > Selected gene-transcript pairs
		awk 'NR==FNR { lines[$0]=1; next } $1"\t"$2 in lines' tmp/GT_comm.txt tmp/GTL_fa.txt | sort -k1,1 -k3,3V | awk          '{if(g!=$1){print line} line=$0;g=$1;}END{print $line}' | tail -n +2 > tmp/GTL_comm_longest.txt
		# Extract the fasta sequence of the selected gene-transcripts
		awk '{if(NR==FNR){l[">"$1"\t"$2]=1; next} if($1 ~ />/){if($1"\t"$2 in l){valid=1}else{valid=0}} if(valid==1){print $0}}' tmp/GTL_comm_longest.txt <(zcat ${ProteomesFolder}/${Species}.pep.all.fa.gz | sed 's/.*gene:\([A-Z0-9]\+\).*transcript:\([A-Z0-9]\+\).*/>\1\t\2/g') | sed 's/\t/|/g' > ${ResultsFolder}/FilteredProteomes/${Species}_ProteinCoding_LongestTx.fa
		# Extract the GTF CDS entries for the selected gene-transcripts 
		awk '{if(NR==FNR){l[$1"\t"$2]=1; next} if($9"\t"$10 in l){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"|"$10"\t"$9"\t"$10}}' tmp/GTL_comm_longest.txt <(zcat ${TranscriptomesFolder}/${Species}.103.gtf.gz | awk '{if($3 ~ /CDS/){print $0}}' | grep 'protein_coding' | sed 's/\tgene_id "\([A-Z0-9]\+\)".*transcript_id "\([A-Z0-9]\+\)".*/\t\1\t\2/g') > ${ResultsFolder}/FilteredGTFs/${Species}_ProteinCoding_LongestTx.gtf
		#rm tmp/GT*.txt
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
	if [[ $(ls ${ResultsFolder}/FilteredProteomes/${Species}_ProteinCoding_inGTF_CleanFormat_LongestProt.fa* 2> ~/null | wc -l) -lt 1 ]] || [[ $(ls ${ResultsFolder}/FilteredGTFs/${Species}_ProteinCoding_inGTF_CleanFormat_LongestProt.gtf* 2> ~/null | wc -l) -lt 1 ]]
	then
		echo "Filtering proteome"
		SpeciesShortName=$(echo ${Species^^} | awk -F '_' '{print substr($1, 1, 1)substr($2, 1, 3)}')
		# Extract info from GFF (filter for CDS & clean format of gene ID and protein ID)
		#zcat ${TranscriptomesFolder}/${Species}.gff.gz | grep -v '^#' | awk '{if($3=="CDS"){print $0}}' | sed 's/ID.*;gene=\([A-Za-z0-9\.\/\-]\+\);.*;protein_id=\([A-Za-z0-9_\.]\+\)$/\1\t\2/g' | sed 's/ID.*;gene=\([A-Za-z0-9\.\/\-]\+\);.*;protein_id=\([A-Za-z0-9_\.]\+\);.*/\1\t\2/g' | awk -v sn=${SpeciesShortName} '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"sn"_"$9"\t"sn"_"$10}' > ${ResultsFolder}/${Species}_ProteinCoding_inGTF_CleanFormat.txt
		## Extract protein sequences present in filtered gff and add geneID to protID
		#awk '{if(NR==FNR){a[">"$10]=1;b[">"$10]=$9;next} if($1 ~ /^>/){if(a[$1]==1){valid=1; print $1"|"b[$1]}else{valid=0}}else{if(valid==1){print $0}}}' ${ResultsFolder}/${Species}_ProteinCoding_inGTF_CleanFormat.txt <(zcat ${ProteomesFolder}/${Species}.faa.gz | awk -v sn=${SpeciesShortName} '{if($1 ~ /^>/){print ">"sn"_"substr($1, 2, length($1))}else{print $0}}') > ${ResultsFolder}/FilteredProteomes/${Species}_ProteinCoding_inGTF_CleanFormat.fa
		## Extract longest protein for each gene
		#awk '{if($1 ~ />/){print name"\t"len; name=$1; len=0; next;} len=len+length($1);}END{print name"\t"len;}' ${ResultsFolder}/FilteredProteomes/${Species}_ProteinCoding_inGTF_CleanFormat.fa | tail -n +2 | sed 's/>//g' | sed 's/|/\t/g' | sort -k2,2 -k3,3V | awk '{if(g==$2){p=$1}else{print g"\t"p;g=$2;p=$1}}END{print g"\t"p;}' | tail -n +2 > ${ResultsFolder}/${Species}_ProteinCoding_inGTF_CleanFormat_LongestProt.txt
		## Extract sequences of the longest protein per gene
		#awk '{if(NR==FNR){a[">"$2"|"$1]=1;next;} if($1 ~ /^>/){if(a[$1]==1){valid=1}else{valid=0}}if(valid==1){print $0}}' ${ResultsFolder}/${Species}_ProteinCoding_inGTF_CleanFormat_LongestProt.txt ${ResultsFolder}/FilteredProteomes/${Species}_ProteinCoding_inGTF_CleanFormat.fa > ${ResultsFolder}/FilteredProteomes/${Species}_ProteinCoding_inGTF_CleanFormat_LongestProt.fa
		## Remove unnecessary intermediate files
		#rm ${ResultsFolder}/FilteredProteomes/${Species}_ProteinCoding_inGTF_CleanFormat.fa
		## Extract GTF info of the selected genes
		#awk '{if(NR==FNR){a[$1]=1;next;} if(a[$9]==1){print $0}}' <(cat ${ResultsFolder}/FilteredProteomes/${Species}_ProteinCoding_inGTF_CleanFormat_LongestProt.fa | grep '>' | sed 's/>//g') <(cat ${ResultsFolder}/${Species}_ProteinCoding_inGTF_CleanFormat.txt | awk '{if($3 == "CDS"){print $0}}' | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$10"|"$9"\t"$9"\t"$10}') > ${ResultsFolder}/FilteredGTFs/${Species}_ProteinCoding_inGTF_CleanFormat_LongestProt.gtf
	fi
done

################################################
###
### EXTRACTING DNA SEQUENCES
###
################################################
echo "### EXTRACTING DNA SEQUENCES"
#for Species in Branchiostoma_lanceolatum.BraLan3 Homo_sapiens.GRCh38 Mus_musculus.GRCm39 Danio_rerio.GRCz11 Gallus_gallus.GRCg6a Branchiostoma_belcheri.Haploidv18h27 Branchiostoma_floridae.Bfl_VNyyK Strongylocentrotus_purpuratus.Spur5.0 Asterias_rubens.eAstRub1.3 Saccoglossus_kowalevskii.Skow1.1
#for Species in Homo_sapiens.GRCh38 Mus_musculus.GRCm39 Danio_rerio.GRCz11 Gallus_gallus.GRCg6a
for Species in Homo_sapiens.GRCh38
do
	echo ${Species}
	Genome=$(ls ${GenomesFolder}/${Species}* | grep -v '.fai')
	GTF=$(ls ${ResultsFolder}/FilteredGTFs/${Species}*)
	if [[ $(echo ${GTF} | sed 's/.*\.//g') =~ "gz" ]];	then gunzip ${GTF}; fi
	GTF=$(ls ${ResultsFolder}/FilteredGTFs/${Species}*.gtf)
	ProteomeTMP=$(ls ${ResultsFolder}/FilteredProteomes/${Species}*)
	if [[ $(echo ${ProteomeTMP} | sed 's/.*\.//g') =~ "gz" ]];	then gunzip ${ProteomeTMP}; fi
	ProteomeTMP=$(ls ${ResultsFolder}/FilteredProteomes/${Species}*.fa)

	# Extract DNA sequences of proteins from the genome reading GTF
	if [[ ! -s ${ResultsFolder}/FilteredDNASequences/${Species}.fa ]] && [[ ! -s ${ResultsFolder}/FilteredDNASequences/${Species}.fa.gz ]] 
	then
		echo "Extracting sequences"
		format=$(echo ${Genome} | sed 's/.*\.//g')
		if [[ ${format} =~ "gz" ]];	then gunzip ${Genome}; fi
		Genome=$(ls ${GenomesFolder}/${Species}* | grep -v '.fai')
		rm ${ResultsFolder}/FilteredDNASequences/${Species}.fa 2> ~/null
		for id in $(cat ${ProteomeTMP} | grep '>' | sed 's/>//g')
		do
			echo '>'${id} >> ${ResultsFolder}/FilteredDNASequences/${Species}.fa
			# Use bedtoos getfasta to extract sequences in a bed format (greped and formated from the GTF)
			PrimarySeq=$(${bedtools} getfasta -fi ${Genome} -bed <(grep -w ${id} <(cat ${GTF}) | awk '{if($3 ~ /CDS/){print $1"\t"$4-1"\t"$5}}' | sort -k2,2V) | grep -v '>' | awk '{seq=seq$0}END{print seq}')
			Strand=$(cat ${GTF} | grep -w ${id} | head -1 | cut -f7)
			if [[ ${Strand} == "+" ]]
			then
				echo ${PrimarySeq}	>> ${ResultsFolder}/FilteredDNASequences/${Species}.fa
			else
				# Reverse complement
				echo ${PrimarySeq} | rev | sed 's/A/Z/g' | sed 's/T/A/g' | sed 's/Z/T/g' | sed 's/C/Z/g' | sed 's/G/C/g' | sed 's/Z/G/g' >> ${ResultsFolder}/FilteredDNASequences/${Species}.fa
			fi
		done
	fi
done

################################################
###
### CHECKING DNA - AA SEQUENCES CORRESPONDENCE
###
################################################
echo "### CHECKING DNA - AA SEQUENCES CORRESPONDENCE"
#for Species in Branchiostoma_lanceolatum.BraLan3 Homo_sapiens.GRCh38 Mus_musculus.GRCm39 Danio_rerio.GRCz11 Gallus_gallus.GRCg6a Branchiostoma_belcheri.Haploidv18h27 Branchiostoma_floridae.Bfl_VNyyK Strongylocentrotus_purpuratus.Spur5.0 Asterias_rubens.eAstRub1.3 Saccoglossus_kowalevskii.Skow1.1
#for Species in Homo_sapiens.GRCh38 Mus_musculus.GRCm39 Danio_rerio.GRCz11 Gallus_gallus.GRCg6a
for Species in Homo_sapiens.GRCh38
do
	echo ${Species}
	GTF=$(ls ${ResultsFolder}/FilteredGTFs/${Species}*)
	if [[ $(echo ${GTF} | sed 's/.*\.//g') =~ "gz" ]];	then gunzip ${GTF}; fi
	GTF=$(ls ${ResultsFolder}/FilteredGTFs/${Species}*.gtf)
	ProteomeTMP=$(ls ${ResultsFolder}/FilteredProteomes/${Species}*)
	if [[ $(echo ${ProteomeTMP} | sed 's/.*\.//g') =~ "gz" ]];	then gunzip ${ProteomeTMP}; fi
	ProteomeTMP=$(ls ${ResultsFolder}/FilteredProteomes/${Species}*.fa)

	if [[ ! -s ${ResultsFolder}/Checking_DNA_AA_sequences/${Species}_DNA_AA_seq_check.txt ]]
	then
		# Check that extracted DNA sequences correspond to protein sequences
		echo "Comparing DNA and AA sequences"
		for id in $(cat ${ProteomeTMP} | grep '>' | sed 's/>//g')
		do
			#echo ${id}
			ProtSeq=$(awk -v id=${id} '{if($1 ~ />/){if($1 == ">"id){valid=1;next;}else{valid=0}} if(valid == 1){print $0}}' ${ProteomeTMP} | awk '{seq=seq""$1}END{print seq}')
			cDNASeq=$(awk -v id=${id} '{if($1 ~ />/){if($1 == ">"id){valid=1;next;}else{valid=0}} if(valid == 1){print $0}}' ${ResultsFolder}/FilteredDNASequences/${Species}.fa | awk '{seq=seq""$1}END{print seq}')
			chr=$(cat ${GTF} | grep -w ${id} | head -1 | cut -f1)
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
	fi
done



################################################
###
### Creating final proteome & CDS sequences files & GTF file 
###
################################################
echo "### Creating final proteome & CDS sequences files & GTF file"
#for Species in Branchiostoma_lanceolatum.BraLan3 Homo_sapiens.GRCh38 Mus_musculus.GRCm39 Danio_rerio.GRCz11 Gallus_gallus.GRCg6a Branchiostoma_belcheri.Haploidv18h27 Branchiostoma_floridae.Bfl_VNyyK Strongylocentrotus_purpuratus.Spur5.0 Asterias_rubens.eAstRub1.3 Saccoglossus_kowalevskii.Skow1.1
#for Species in Homo_sapiens.GRCh38 Mus_musculus.GRCm39 Danio_rerio.GRCz11 Gallus_gallus.GRCg6a
for Species in Homo_sapiens.GRCh38
do
	echo ${Species}
	GTFTMP=$(ls ${ResultsFolder}/FilteredGTFs/${Species}*)
	if [[ $(echo ${GTFTMP} | sed 's/.*\.//g') =~ "gz" ]];	then gunzip ${GTFTMP}; fi
	GTFTMP=$(ls ${ResultsFolder}/FilteredGTFs/${Species}*gtf)
	ProteomeTMP=$(ls ${ResultsFolder}/FilteredProteomes/${Species}*)
	if [[ $(echo ${ProteomeTMP} | sed 's/.*\.//g') =~ "gz" ]];	then gunzip ${ProteomeTMP}; fi
	ProteomeTMP=$(ls ${ResultsFolder}/FilteredProteomes/${Species}*fa)
	DNAseqTMP=$(ls ${ResultsFolder}/FilteredDNASequences/${Species}*)
	if [[ $(echo ${DNAseqTMP} | sed 's/.*\.//g') =~ "gz" ]];	then gunzip ${DNAseqTMP}; fi
	DNAseqTMP=$(ls ${ResultsFolder}/FilteredDNASequences/${Species}*fa)
	CheckFile=${ResultsFolder}/Checking_DNA_AA_sequences/${Species}_DNA_AA_seq_check.txt

	if [[ $(ls ${ResultsFolder}/CheckedProteomes/${Species}_AA.fa* 2> ~/null | wc -l) -lt 1 ]] || [[ $(ls ${ResultsFolder}/CheckedGTFs/${Species}.gtf* 2> ~/null | wc -l) -lt 1 ]] || [[ $(ls ${ResultsFolder}/CheckedDNASequences/${Species}_DNA.fa* 2> ~/null | wc -l) -lt 1 ]]
	then
		echo "Creating final files"
		# Extract the protein fasta sequence of the checked genes
		awk '{if(NR==FNR){l[">"$1]=1; next} if($1 ~ />/){if($1 in l){valid=1}else{valid=0}} if(valid==1){print $0}}' <(grep -v 'Error' ${CheckFile} | cut -f1) ${ProteomeTMP} > ${ResultsFolder}/CheckedProteomes/${Species}_AA.fa
		# Extract the GTF CDS entries for the selected gene-transcripts 
		awk '{if(NR==FNR){l[$1]=1; next} if($9 in l){print $0}}' <(grep -v 'Error' ${CheckFile} | cut -f1) ${GTFTMP} > ${ResultsFolder}/CheckedGTFs/${Species}.gtf
		# Extract the DNA fasta sequence of the checked genes
		awk '{if(NR==FNR){l[">"$1]=1; next} if($1 ~ />/){if($1 in l){valid=1}else{valid=0}} if(valid==1){print $0}}' <(grep -v 'Error' ${CheckFile} | cut -f1) ${DNAseqTMP} > ${ResultsFolder}/CheckedDNASequences/${Species}_DNA.fa
	fi
done



#rm -r ${ResultsFolder}/Filtered*
#gzip ${ResultsFolder}/CheckedProteomes/*.fa ${ResultsFolder}/CheckedDNASequences/*.fa ${ResultsFolder}/CheckedGTFs/*.gtf


 















