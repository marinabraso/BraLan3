#!/bin/bash




# Scripts
Rscript="Analysis/ComparativeGeneExpression/ComparativeGeneExpression.R"

# Files & parameters
ResultsFolder="Plots/ComparativeGeneExpression"
mkdir -p ${ResultsFolder}

################################################
### 

${Rscript} ${ResultsFolder} ${Species1} ${Species2}



