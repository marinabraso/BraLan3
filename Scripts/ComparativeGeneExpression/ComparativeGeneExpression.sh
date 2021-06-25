#!/bin/bash

module load curl/7.74.0


# Scripts
Rscript="Scripts/ComparativeGeneExpression/ComparativeGeneExpression.R"

# Files & parameters
Species1="Blan"
Species2="Drer"
ResultsFolder="Results/ComparativeGeneExpression"
mkdir -p ${ResultsFolder}

################################################
### 

${Rscript} ${ResultsFolder} ${Species1} ${Species2}



