#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")
library(ape)
args <- commandArgs(trailingOnly=TRUE)
TreeFolder <- args[1]
ResultsFolder <- args[2]
strlist <- args[3]
type <- args[4]
list <- unlist(strsplit(strlist, ";"))
print(length(list))

pdf(paste(ResultsFolder, "/VertebratesMonophyletic", type, ".pdf", sep=""), width=15, height=10)

VertMonophyletic <- c()
for(og in list){
	#print(og)
	if(file.exists(paste0(TreeFolder, "/RAxML_result.", og, "_PROTGAMMAAUTO"))){
		treetext <- readLines(paste0(TreeFolder, "/RAxML_result.", og, "_PROTGAMMAAUTO"))
		tree <- read.tree(text=treetext)
		VertMonophyletic <- c(VertMonophyletic, is.monophyletic(tree, tree$tip.label[grep("ENS", tree$tip.label)], plot=TRUE))
	}else{
		VertMonophyletic <- c(VertMonophyletic, NA)		
	}
}
t <- as.data.frame(table(VertMonophyletic))
cat("% of trees where vertebrates are monophyletic\n")
t$Freq[which(t$VertMonophyletic==TRUE)]/length(na.omit(VertMonophyletic))*100
cat("% of trees where vertebrates aren't monophyletic\n")
t$Freq[which(t$VertMonophyletic==FALSE)]/length(na.omit(VertMonophyletic))*100

cat("Total computed trees\n")
length(na.omit(VertMonophyletic))
cat("Total alignments without trees\n")
length(VertMonophyletic)-length(na.omit(VertMonophyletic))
cat("% of computed trees of the total OG in the category\n")
length(na.omit(VertMonophyletic))/length(VertMonophyletic)*100

dev.off()








