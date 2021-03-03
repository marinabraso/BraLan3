#!/bin/bash

# Scripts
godon="Scripts/dNdSBetweenParalogs/godon"
bedtools="/scratch/wally/FAC/FBM/DEE/mrobinso/default/mbrasovi/Software/bedtools.static.binary"

# Files & parameters
OrthologsFolder="Results/FindOrthologs"
ResultsFolder="Results/dNdSBetweenParalogs"
mkdir -p ${ResultsFolder}/Sequences
Genome="Data/Genomes/BraLan3/Bla_gDNA-FINAL.fa"
GTF="Data/Transcriptomes/BranchiostomaLanceolatum_BraLan3.gtf.gz"

################################################
###


#${godon} M0 --out-tree ${ResultsFolder}/M0_tree.nwk EMGT00050000000025.Drosophila.001.fst EMGT00050000000025.Drosophila.001.nwk













