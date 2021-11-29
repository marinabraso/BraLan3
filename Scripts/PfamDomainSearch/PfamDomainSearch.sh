#!/bin/bash


ProteomesFolder="Results/FilteringGeneSets/ProteomesAmphVerteb"
ResultsFolder="Results/PfamDomainSearch"
mkdir -p ${ResultsFolder}

Ethresh=0.001
PfamAHmm=../Software/Pfam-A/Pfam-A.hmm # from 16-Nov-2021
# Pfam release       : 35.0
# Pfam-A families    : 19632
# Date               : 2021-11
# Based on UniProtKB : 2021_03
hmmscan=/work/FAC/FBM/DEE/mrobinso/default/mbrasovi/Software/hmmer-3.3.2/src/hmmscan

#Species="Branchiostoma_lanceolatum.BraLan3 Gallus_gallus.GRCg6a Mus_musculus.GRCm39 Danio_rerio.GRCz11 Homo_sapiens.GRCh38"
Species="Danio_rerio.GRCz11"

for sp in $Species; do
    echo $sp
    ${hmmscan} -o ${ResultsFolder}/${sp}.out --noali --tblout ${ResultsFolder}/${sp}.tbl -E ${Ethresh} ${PfamAHmm} ${ProteomesFolder}/${sp}.fa
done



# sbatch -t 8:00:00 --mem=8000 -J PfamHmm -o tmp/PfamHmmd.out -e tmp/PfamHmmd.err Scripts/PfamDomainSearch/PfamDomainSearch.sh

