#!/usr/bin/env Rscript


library(tximport)
library(rhdf5)


# Files, folders and parameters
args <- commandArgs(trailingOnly=TRUE)
AnnotationsTrasncToGene <- args[1]
ResultsFolder <- args[2]

gene.to.tx <- read.table(AnnotationsTrasncToGene, h=F)
tx2gene <- as.data.frame(cbind(as.character(gene.to.tx[,2]), as.character(gene.to.tx[,1])))
system_out <- system(paste("ls ", ResultsFolder,"/kallisto/*/abundance.h5", sep =""), intern=T)
tpmsamplefiles <- read.table(text=system_out, h=F)
tpmsamplefiles <- as.character(tpmsamplefiles[,1])
names <- strsplit(tpmsamplefiles, "/")
names(tpmsamplefiles) <- sapply(c(1:length(names)), function(x){return(unlist(names[x])[4])})

tx.imported <- tximport(files=tpmsamplefiles, type="kallisto", tx2gene=tx2gene, ignoreTxVersion=TRUE)
write.table(tx.imported$counts, file = paste(ResultsFolder,"/Gene_Counts_kallisto_tximport.tab", sep =""), quote = F, sep="\t", col.names = NA, row.names = TRUE)
write.table(tx.imported$abundance, file = paste(ResultsFolder,"/Gene_TPM_kallisto_tximport.tab", sep =""), quote = F, sep="\t", col.names = NA, row.names = TRUE)
write.table(tx.imported$length, file = paste(ResultsFolder,"/Gene_Length_kallisto_tximport.tab", sep =""), quote = F, sep="\t", col.names = NA, row.names = TRUE)








