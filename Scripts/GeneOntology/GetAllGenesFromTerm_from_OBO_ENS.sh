#!/bin/bash

# Scripts

# Files & parameters
GO=$1
genesoutfile=$2
OBO=$3
ENSAnnot=$4

################################################
### 
ChildTermList=${GO}
newNumChildTerm=$(echo ${ChildTermList} | awk -F ' ' '{for(i=1;i<=NF;i++){print $i}}' | sort -u | wc -l)
oldNumChildTerm=0
while [  ${newNumChildTerm} -gt ${oldNumChildTerm} ]; do
	oldNumChildTerm=${newNumChildTerm}
	newChildTermList=$(awk -F '\t' '{if(NR==FNR){a[$1]=1;next} if($1 ~ /^id: /){id=substr($1, 5, length($1));next} if($1 ~ /^is_a: /){if(a[substr($1, 7, 10)]){print id}}}' <(echo ${ChildTermList} | awk -F ' ' '{for(i=1;i<=NF;i++){print $i}}' | sort -u) ${OBO})
	ChildTermList=$(echo ${ChildTermList}" "${newChildTermList} | awk -F ' ' '{for(i=1;i<=NF;i++){print $i}}' | sort -u)
	newNumChildTerm=$(echo ${ChildTermList} | awk -F ' ' '{for(i=1;i<=NF;i++){print $i}}' | sort -u | wc -l)
done
echo ${ChildTermList}
Genes=$(awk '{if(NR==FNR){a[$1]=1;next} if(a[$2]){print $1}}' <(echo ${ChildTermList} | awk -F ' ' '{for(i=1;i<=NF;i++){print $i}}') ${ENSAnnot} | sort | uniq)
numGenes=$(echo ${Genes} | awk -F ' ' '{for(i=1;i<=NF;i++){print $i}}' | sort -u | wc -l)
if [[ numGenes -gt 0 ]]; then
	echo ${Genes} | awk -F ' ' '{for(i=1;i<=NF;i++){print $i}}' | sort -u > ${genesoutfile}
fi
