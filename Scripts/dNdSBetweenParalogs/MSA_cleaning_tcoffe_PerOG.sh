#!/bin/bash

# Scripts
module add tcoffee/13.45.0.4846264
module add blast-plus/2.11.0

export MAX_N_PID_4_TCOFFEE=10000000

group=$1

t_coffee -infile=../MSA_mafft/${group}_AA.fa  -mode evaluate  -output=score_ascii
seqScore=$(cat ${group}_AA.score_ascii | grep 'cons' | grep -v ':' | sed 's/ \+/\t/g' | awk '{seq=seq$2}END{print seq}')
echo ${seqScore} > ${group}_AA.score
echo ${seqScore} | awk '{for(i=1;i<=length($0);i++){l=substr($0,i,1); seq=seq""l""l""l}}END{print seq}' > ${group}_DNA.score
awk '{if(NR==FNR){for(i=1;i<=length($0);i++){a[i]=substr($0,i,1)};next} if($1 ~ />/){print $0; seq=""}else{for(i=1;i<=length($0);i++){if(a[i]>5 || substr($0,i,1)=="-"){seq=seq""substr($0,i,1)}else{seq=seq"N"}}print seq}}' ${group}_DNA.score ../MSA_mafft/${group}_DNA.fa > ${group}_DNA_clean.fa
