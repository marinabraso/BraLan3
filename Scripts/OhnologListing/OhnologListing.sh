#!/bin/bash


DataFolder="Data/Ohnologs"
GTFsFolder="Results/FilteringGeneSets/FilteredGTFs"
OG2genes="Results/FindOrthologs/AmphVerteb_broccoli/dir_step3/table_OGs_protein_names.txt"
ResultsFolder="Results/OhnologListing"
mkdir -p ${ResultsFolder}

# From original pairs files
cat ${DataFolder}/hsapiens.Pairs.Strict.2R.txt | tail -n +2 | awk '{print $1; print $2}' > ${ResultsFolder}/ListAllGenes2R_pairs.txt
cat ${DataFolder}/drerio.Pairs.Strict.2R.txt | tail -n +2 | awk '{print $1; print $2}' >> ${ResultsFolder}/ListAllGenes2R_pairs.txt
cat ${DataFolder}/ggallus.Pairs.Strict.2R.txt | tail -n +2 | awk '{print $1; print $2}' >> ${ResultsFolder}/ListAllGenes2R_pairs.txt
cat ${DataFolder}/mmusculus.Pairs.Strict.2R.txt | tail -n +2 | awk '{print $1; print $2}' >> ${ResultsFolder}/ListAllGenes2R_pairs.txt
cat ${DataFolder}/drerio.Pairs.Strict.3R.txt | tail -n +2 | awk '{print $1; print $2}' > ${ResultsFolder}/ListAllGenes3R_pairs.txt

awk '{if(FNR==NR){a[$1]=1;next;}if(a[$2]){print $1}}' ${ResultsFolder}/ListAllGenes2R_pairs.txt <(tail -n +2 ${OG2genes} | awk '{for(i=2;i<=NF;i++){print $1"\t"$i}}') | sort -u > ${ResultsFolder}/2R_Strict_OG_pairs.txt
awk '{if(FNR==NR){a[$1]=1;next;}if(a[$2]){print $1}}' ${ResultsFolder}/ListAllGenes3R_pairs.txt <(tail -n +2 ${OG2genes} | awk '{for(i=2;i<=NF;i++){print $1"\t"$i}}') | sort -u > ${ResultsFolder}/3R_Strict_OG_pairs.txt

# From original families files
cat ${GTFsFolder}/Danio_rerio.GRCz11.gtf | grep 'gene_name' | sed 's/.*gene_id "\(\S\+\)";.*gene_name "\(\S\+\)";.*/\1\t\2/g' | sort -u > ${ResultsFolder}/Drer_EnsID2names.list
cat ${GTFsFolder}/Gallus_gallus.GRCg6a.gtf | grep 'gene_name' | sed 's/.*gene_id "\(\S\+\)";.*gene_name "\(\S\+\)";.*/\1\t\2/g' | sort -u > ${ResultsFolder}/Ggal_EnsID2names.list
cat ${GTFsFolder}/Homo_sapiens.GRCh38.gtf | grep 'gene_name' | sed 's/.*gene_id "\(\S\+\)";.*gene_name "\(\S\+\)";.*/\1\t\2/g' | sort -u > ${ResultsFolder}/Hsap_EnsID2names.list
cat ${GTFsFolder}/Mus_musculus.GRCm39.gtf | grep 'gene_name' | sed 's/.*gene_id "\(\S\+\)";.*gene_name "\(\S\+\)";.*/\1\t\2/g' | sort -u > ${ResultsFolder}/Mmus_EnsID2names.list

cat ${DataFolder}/drerio.Families.Strict.3R.txt | sed 's/ \S /\t/g' | sed 's/ (1 to many)//g' | awk '{for(i=2;i<=NF;i++){print $i}}' > ${ResultsFolder}/Drer_3R_names_families.list
cat ${DataFolder}/drerio.Families.Strict.2R.txt | sed 's/ \S /\t/g' | sed 's/ (1 to many)//g' | awk '{for(i=2;i<=NF;i++){print $i}}' > ${ResultsFolder}/Drer_2R_names_families.list
cat ${DataFolder}/ggallus.Families.Strict.2R.txt | sed 's/ \S /\t/g' | sed 's/ (1 to many)//g' | awk '{for(i=2;i<=NF;i++){print $i}}' > ${ResultsFolder}/Ggal_2R_names_families.list
cat ${DataFolder}/hsapiens.Families.Strict.2R.txt | sed 's/ \S /\t/g' | sed 's/ (1 to many)//g' | awk '{for(i=2;i<=NF;i++){print $i}}' > ${ResultsFolder}/Hsap_2R_names_families.list
cat ${DataFolder}/mmusculus.Families.Strict.2R.txt | sed 's/ \S /\t/g' | sed 's/ (1 to many)//g' | awk '{for(i=2;i<=NF;i++){print $i}}' > ${ResultsFolder}/Mmus_2R_names_families.list

awk '{if(FNR==NR){a[$1]=1;next;}if(a[$2]){print $1}}' ${ResultsFolder}/Drer_3R_names_families.list ${ResultsFolder}/Drer_EnsID2names.list > ${ResultsFolder}/3R_EnsID_families.list
awk '{if(FNR==NR){a[$1]=1;next;}if(a[$2]){print $1}}' ${ResultsFolder}/Drer_2R_names_families.list ${ResultsFolder}/Drer_EnsID2names.list > ${ResultsFolder}/2R_EnsID_families.list
awk '{if(FNR==NR){a[$1]=1;next;}if(a[$2]){print $1}}' ${ResultsFolder}/Drer_2R_names_families.list ${ResultsFolder}/Ggal_EnsID2names.list >> ${ResultsFolder}/2R_EnsID_families.list
awk '{if(FNR==NR){a[$1]=1;next;}if(a[$2]){print $1}}' ${ResultsFolder}/Drer_2R_names_families.list ${ResultsFolder}/Hsap_EnsID2names.list >> ${ResultsFolder}/2R_EnsID_families.list
awk '{if(FNR==NR){a[$1]=1;next;}if(a[$2]){print $1}}' ${ResultsFolder}/Drer_2R_names_families.list ${ResultsFolder}/Mmus_EnsID2names.list >> ${ResultsFolder}/2R_EnsID_families.list

awk '{if(FNR==NR){a[$1]=1;next;}if(a[$2]){print $1}}' ${ResultsFolder}/3R_EnsID_families.list <(tail -n +2 ${OG2genes} | awk '{for(i=2;i<=NF;i++){print $1"\t"$i}}') | sort -u > ${ResultsFolder}/3R_Strict_OG_families.txt
awk '{if(FNR==NR){a[$1]=1;next;}if(a[$2]){print $1}}' ${ResultsFolder}/2R_EnsID_families.list <(tail -n +2 ${OG2genes} | awk '{for(i=2;i<=NF;i++){print $1"\t"$i}}') | sort -u > ${ResultsFolder}/2R_Strict_OG_families.txt



