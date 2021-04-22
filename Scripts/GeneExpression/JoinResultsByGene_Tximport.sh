#!/bin/bash


module add R/3.5.1;

# Scripts
TximportRScript=Scripts/GeneExpression/JoinResultsByGene_Tximport.R

# Files & parameters
ResultsFolder=Results/GeneExpression
Basename="Branchiostoma_lanceolatum.BraLan3"
AnnotationsTrasncToGene=${ResultsFolder}/${Basename}_tx2gene.txt


${TximportRScript} ${AnnotationsTrasncToGene} ${ResultsFolder}


