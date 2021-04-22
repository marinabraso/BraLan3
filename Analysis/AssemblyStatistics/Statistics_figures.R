#!/usr/bin/env Rscript



######################################################################
# Libraries & functions

`%out%` <- Negate(`%in%`)
script <- sub(".*=", "", commandArgs()[4])
source(paste(substr(script,1, nchar(script)-2), "_functions.R", sep=""))

######################################################################
# Files & folders

ResultsFolder <- "Plots/AssemblyStatistics"
StatisticsFile <- "Results/AssemblyStatistics/Joined_Statistics_t.txt"

######################################################################
# General parameters

Verteb <- c("Drer", "Ggal", "Mmus", "Hsap")
Amphi <- c("Bflo", "Bbel", "Blan3", "Blan2")
OutDeut <- c("Spur", "Arub", "Skow")
SortedSpecies <- c(Verteb, OutDeut, Amphi)
color <- "black"
###########################################################################
# Read data
Stat <- read.table(StatisticsFile, h=T, sep = "\t", row.names=1)
shortnames <- c()
for(name in rownames(Stat)){
	vecname <- unlist(strsplit(name, "_"))
	shortnames <- c(shortnames,  paste0(substr(vecname[1], 1, 1), substr(vecname[2], 1, 3)))
}
Stat$ShortName <- shortnames
Stat$ShortName[which(Stat$ShortName == "Blan")] <- "Blan3"
Stat$ShortName[which(Stat$ShortName == "BBra")] <- "Blan2"
Stat <- Stat[order(order(match(SortedSpecies, Stat$ShortName))),]
print(Stat$ShortName)
# Choose color
Stat$Color <- "grey50"
#Stat$Color[which(Stat$ShortName %in% Verteb)] <- "gold"
#Stat$Color[which(Stat$ShortName %in% OutDeut)] <- "black"
Stat$Color[which(Stat$ShortName == "Blan3")] <- "black"
Stat$Color[which(Stat$ShortName == "Blan2")] <- "black"

###########################################################################
###########################################################################
# Plotting
###########################################################################
###########################################################################
pdf(paste(ResultsFolder, "/AssembliesComparison.pdf", sep=""), width=20, height=10)
par(mar=c(10,10,5,5),oma=c(1,1,1,1))
layout(matrix(c(1,2),nrow=1,ncol=2,byrow=T), widths=c(1,1), heights=c(1), TRUE)


Scatter_plot(Stat$prop1M*100, Stat$BUSCOproteinmetazoa, c(0,100), c(70,100), Stat$ShortName, Stat$Color, "% of genome in >1M sequences", "BUSCO completeness score")
Scatter_plot(Stat$N50/1000000, Stat$BUSCOproteinmetazoa, c(0,150), c(70,100), Stat$ShortName, Stat$Color, "N50 (M)", "BUSCO completeness score")










