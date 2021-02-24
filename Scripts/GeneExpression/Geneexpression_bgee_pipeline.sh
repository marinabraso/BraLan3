#!/bin/bash



# Scripts
${bgee_RNA_Seq}="bgee_pipeline/pipeline/RNA_Seq"
${prepare_GTF_script}=${bgee_RNA_Seq}"/0Before/prepare_GTF.R"


# Files & parameters
Basename="BranchiostomaLanceolatum_BraLan3"
ResultsFolder="Results/GeneExpression/"${Basename}
mkdir -p ${ResultsFolder}

OriginalGTF="Data/Transcriptomes/${Basename}.gtf.gz"

################################################
### Prepare GTF
R CMD BATCH --no-save --no-restore '--args gene_gtf_path=${OriginalGTF} output_gtf_path=${ResultsFolder}/${Basename}' ${prepare_GTF_script}











