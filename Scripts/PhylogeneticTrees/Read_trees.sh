#!/bin/bash

# Scripts
RScript="Scripts/PhylogeneticTrees/Read_trees.R"

# Files & parameters
ResultsFolder="Results/PhylogeneticTrees"
TreeFolder=${ResultsFolder}/RAxML_Trees

echo "Small-scale & small-scale"
listSSDSSD=$(tail -n +2 Plots/DuplicatesOntologyExpression/OGInfo.txt | cut -f1,10,12 | sed 's/"//g' | grep -P 'Small.*Small' | cut -f1)
listSSDSSD=$(echo ${listSSDSSD} | sed 's/ /;/g' )
${RScript} ${TreeFolder} ${ResultsFolder} ${listSSDSSD} SSDSSD 2> ~/null

echo "Ohnologs & small-scale"
listOSSD=$(tail -n +2 Plots/DuplicatesOntologyExpression/OGInfo.txt | cut -f1,10,12 | sed 's/"//g' | grep -P 'Ohno.*Small' | cut -f1)
listOSSD=$(echo ${listOSSD} | sed 's/ /;/g' )
${RScript} ${TreeFolder} ${ResultsFolder} ${listOSSD} OSSD 2> ~/null

echo "Single-copy & single-copy"
listSCSC=$(tail -n +2 Plots/DuplicatesOntologyExpression/OGInfo.txt | cut -f1,10,12 | sed 's/"//g' | grep -P 'Sing.*Sing' | cut -f1)
listSCSC=$(echo ${listSCSC} | sed 's/ /;/g' )
${RScript} ${TreeFolder} ${ResultsFolder} ${listSCSC} SCSC 2> ~/null

