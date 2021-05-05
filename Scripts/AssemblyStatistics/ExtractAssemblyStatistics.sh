#!/bin/bash


module add UHTS/Analysis/busco/4.1.4
module add SequenceAnalysis/HMM-Profile/hmmer/3.1b2



pwd=$(pwd)
ResultsFolder=Results/AssemblyStatistics
ProteomesFolder=Data/Proteomes
GenomesFolder=Data/Genomes
Basename=$1
if [[ ${Basename} == "" ]]
then
	echo "Specify the genome basename"
	exit 0
fi
mkdir -p ${ResultsFolder}/${Basename}


GenomeFile=$(ls ${GenomesFolder}/${Basename}*fa | head -1)
if [[ ${GenomeFile} == "" ]]
then
	echo "Cannot find proteome file ${GenomeFile}"
	exit 0
fi


#############################################################
## BASIC Parameters
# Check for unrecognized parameters in the genome
#if [[ $(cat Data/Genomes/Branchiostoma_lanceolatum.BraLan3_genome.fa | grep -v '>' | sed 's/[ATCGatcgN]//g' | awk '{len=len+length($0)}END{print len}') -ne 0 ]]
#then
#	echo "Problem: unrecognized character in the genome (not ACGTacgtN)"
#	exit 0
#fi

# Build chr length files (gapped and ungapped)
if [[ ! -s ${ResultsFolder}/${Basename}/${Basename}_lengths.txt ]] || [[ ! -s ${ResultsFolder}/${Basename}/${Basename}_lengths_ungapped.txt ]]
then
	cat ${GenomeFile} | awk '{if($1 ~ />/){print name"\t"len; name=$1; len=0;next}len=len+length($1)}END{print name"\t"len}' | tail -n +2 | grep -v 'MT' | sort -k2,2Vr | awk '{alen=alen+$2; print $0"\t"NR"\t"alen;}' > ${ResultsFolder}/${Basename}/${Basename}_lengths.txt
	cat ${GenomeFile} | sed 's/N//g' | awk '{if($1 ~ />/){print name"\t"len; name=$1; len=0;next}len=len+length($1)}END{print name"\t"len}' | tail -n +2 | sort -k2,2Vr | awk '{alen=alen+$2; print $0"\t"NR"\t"alen;}' > ${ResultsFolder}/${Basename}/${Basename}_lengths_ungapped.txt
fi

#if [[ ! -s ${ResultsFolder}/${Basename}/${Basename}_statistics.txt ]]
#then
	# Total lengths & num of sequences
	echo "Total lengths & num of sequences"
	TotalLength=$(tail -1 ${ResultsFolder}/${Basename}/${Basename}_lengths.txt | cut -f4)
	echo "TotalLength	${TotalLength}" > ${ResultsFolder}/${Basename}/${Basename}_statistics.txt
	TotalLengthUngapped=$(tail -1 ${ResultsFolder}/${Basename}/${Basename}_lengths_ungapped.txt | cut -f4)
	echo "TotalLengthUngapped	${TotalLengthUngapped}" >> ${ResultsFolder}/${Basename}/${Basename}_statistics.txt
	NumSeq=$(cat ${ResultsFolder}/${Basename}/${Basename}_lengths.txt | wc -l)
	echo "NumSeq	${NumSeq}" >> ${ResultsFolder}/${Basename}/${Basename}_statistics.txt
	NumGaps=$(cat ${GenomeFile} | awk '{if($1 ~ />/){print seq; print $0; seq=""}else{seq=seq""$0}}END{print seq}' | grep -v '>' | sed 's/N\+/N/g' | sed 's/[ACGTacgt]\+//g' | awk '{seq=seq""$1}END{print length(seq)}')
	echo "NumGaps	${NumGaps}" >> ${ResultsFolder}/${Basename}/${Basename}_statistics.txt
	NumGaps1M=$(awk '{if(NR==FNR){a[$1]=1;next} if($1 ~ />/){if(a[$1]){valid=1}else{valid=0}}if(valid==1){print $0}}' <(cat ${ResultsFolder}/${Basename}/${Basename}_lengths.txt | awk '{if($2 >= 1000000){print $0}}' | cut -f1) ${GenomeFile} | awk '{if($1 ~ />/){print seq; print $0; seq=""}else{seq=seq""$0}}END{print seq}' | grep -v '>' | sed 's/N\+/N/g' | sed 's/[ACGTacgt]\+//g' | awk '{seq=seq""$1}END{print length(seq)}')
	echo "NumGaps1M	${NumGaps1M}" >> ${ResultsFolder}/${Basename}/${Basename}_statistics.txt

	# Scaffold N50, L50, N90, L90
	echo "Scaffold N50, L50, N90, L90"
	N50=$(awk '{if(NR==FNR){max=$1;next} if($4 > max/2){print $0; exit}}' <(tail -1 ${ResultsFolder}/${Basename}/${Basename}_lengths.txt | cut -f4) ${ResultsFolder}/${Basename}/${Basename}_lengths.txt | cut -f2)
	echo "N50	${N50}" >> ${ResultsFolder}/${Basename}/${Basename}_statistics.txt
	L50=$(awk '{if(NR==FNR){max=$1;next} if($4 > max/2){print $0; exit}}' <(tail -1 ${ResultsFolder}/${Basename}/${Basename}_lengths.txt | cut -f4) ${ResultsFolder}/${Basename}/${Basename}_lengths.txt | cut -f3)
	echo "L50	${L50}" >> ${ResultsFolder}/${Basename}/${Basename}_statistics.txt
	N90=$(awk '{if(NR==FNR){max=$1;next} if($4 > max*9/10){print $0; exit}}' <(tail -1 ${ResultsFolder}/${Basename}/${Basename}_lengths.txt | cut -f4) ${ResultsFolder}/${Basename}/${Basename}_lengths.txt | cut -f2)
	echo "N90	${N90}" >> ${ResultsFolder}/${Basename}/${Basename}_statistics.txt
	L90=$(awk '{if(NR==FNR){max=$1;next} if($4 > max*9/10){print $0; exit}}' <(tail -1 ${ResultsFolder}/${Basename}/${Basename}_lengths.txt | cut -f4) ${ResultsFolder}/${Basename}/${Basename}_lengths.txt | cut -f3)
	echo "L90	${L90}" >> ${ResultsFolder}/${Basename}/${Basename}_statistics.txt

	# Statistics of sequences >1M (putative chromosomes)
	echo "Statistics of sequences >1M (putative chromosomes)"
	num1M=$(awk '{if($2 > 1000000){line=$0;next}else{print line; exit}}' ${ResultsFolder}/${Basename}/${Basename}_lengths.txt | cut -f3)
	echo "num1M	${num1M}" >> ${ResultsFolder}/${Basename}/${Basename}_statistics.txt
	prop1M=$(awk '{if(NR==FNR){max=$1;next} if($2 > 1000000){accum=$4;next}else{print $4/max; exit}}' <(tail -1 ${ResultsFolder}/${Basename}/${Basename}_lengths.txt | cut -f4) ${ResultsFolder}/${Basename}/${Basename}_lengths.txt)
	echo "prop1M	${prop1M}" >> ${ResultsFolder}/${Basename}/${Basename}_statistics.txt

	## GC content and skew
	#echo "GC content and skew"
	#numG=$(cat ${GenomeFile} | grep -v '>' | tr -d -c 'Gg' | awk '{ print length}')
	#numC=$(cat ${GenomeFile} | grep -v '>' | tr -d -c 'Cc' | awk '{ print length}')
	#sumGC=$(( numG + numC ))
	#difGC=$(( numG - numC ))
	#GCcontent=$(echo "scale=4 ; ${sumGC} / ${TotalLengthUngapped}" | bc)
	#echo "GCcontent	${GCcontent}" >> ${ResultsFolder}/${Basename}/${Basename}_statistics.txt
	#GCskew=$(echo "scale=2 ; $difGC / $sumGC" | bc)
	#echo "GCskew	${GCskew}" >> ${ResultsFolder}/${Basename}/${Basename}_statistics.txt

	#############################################################
	## BUSCO score at the protein level with the metazoa database
	if [[ ! -s ${ResultsFolder}/${Basename}/short_summary.specific.metazoa_odb10.${Basename}.txt ]]
	then 
		echo "Running BUSCO score at the protein level with the metazoa database"
		gunzip ${ProteomesFolder}/${Basename}*fa.gz 2> ~/null
		ProteomeFile=$(ls ${ProteomesFolder}/${Basename}*fa | head -1)
		if [[ ${ProteomeFile} == "" ]]
		then
			echo "Cannot find proteome file ${ProteomeFile}"
			exit 0
		fi
		cd ${ResultsFolder}
		busco -f -m protein -l metazoa_odb10 -i ${pwd}/${ProteomeFile} -o ${Basename}
		cd ${pwd}
	fi
	BUSCOResults=$(cat ${ResultsFolder}/${Basename}/short_summary.specific.metazoa_odb10.${Basename}.txt | grep -P 'C:[0-9\.]+%\[S:' | sed 's/\t//g')
	echo "BUSCOproteinmetazoa	${BUSCOResults}" >> ${ResultsFolder}/${Basename}/${Basename}_statistics.txt
#fi



















