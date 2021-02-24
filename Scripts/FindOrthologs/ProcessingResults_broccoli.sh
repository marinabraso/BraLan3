#!/bin/bash

# Scripts

# Files & parameters
ResultsFolder="Results/FindOrthologs"

################################################
### 
grep -w -f <(awk '{if($2>1){print $0}}' ${ResultsFolder}/broccoli/dir_step3/table_OGs_protein_counts.txt | cut -f1) ${ResultsFolder}/broccoli/dir_step3/table_OGs_protein_names.txt | cut -f2 | awk -F ' ' '{if(NR==1){print "#Paralogs\tGeneNames";next;} print NF"\t"$0}' > ${ResultsFolder}/Dre_Paralogs_FromOrthologusGroups.txt
grep -w -f <(awk '{if($3>1){print $0}}' ${ResultsFolder}/broccoli/dir_step3/table_OGs_protein_counts.txt | cut -f1) ${ResultsFolder}/broccoli/dir_step3/table_OGs_protein_names.txt | cut -f3 | awk -F ' ' '{if(NR==1){print "#Paralogs\tGeneNames";next;} print NF"\t"$0}' > ${ResultsFolder}/Bla_Paralogs_FromOrthologusGroups.txt
grep -w -f <(awk '{if($4>1){print $0}}' ${ResultsFolder}/broccoli/dir_step3/table_OGs_protein_counts.txt | cut -f1) ${ResultsFolder}/broccoli/dir_step3/table_OGs_protein_names.txt | cut -f4 | awk -F ' ' '{if(NR==1){print "#Paralogs\tGeneNames";next;} print NF"\t"$0}' > ${ResultsFolder}/Hsa_Paralogs_FromOrthologusGroups.txt
grep -w -f <(awk '{if($5>1){print $0}}' ${ResultsFolder}/broccoli/dir_step3/table_OGs_protein_counts.txt | cut -f1) ${ResultsFolder}/broccoli/dir_step3/table_OGs_protein_names.txt | cut -f5 | awk -F ' ' '{if(NR==1){print "#Paralogs\tGeneNames";next;} print NF"\t"$0}' > ${ResultsFolder}/Mmu_Paralogs_FromOrthologusGroups.txt


wc -l ${ResultsFolder}/*_Paralogs_FromOrthologusGroups.txt










