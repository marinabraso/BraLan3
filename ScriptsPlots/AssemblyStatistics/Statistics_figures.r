#!/usr/bin/env Rscript


######################################################################
# Libraries & functions

script <- sub(".*=", "", commandArgs()[4])
source(paste(substr(script,1, nchar(script)-2), "_functions.r", sep=""))

######################################################################
# Files & folders
ResultsFolder <- "Plots/AssemblyStatistics"
StatisticsFile <- "Results/AssemblyStatistics/Joined_Statistics_t.txt"

######################################################################
# General parameters
Verteb <- c("Drer", "Ggal", "Mmus", "Hsap")
Amphi <- c("Bflo", "Bbel", "Blan3", "Blan2")
OutDeut <- c("Spur", "Arub", "Skow")

###########################################################################
# Read data
Stat <- read.table(StatisticsFile, h=T, sep = "\t", row.names=1)
Stat$ShortName <- unlist(lapply(rownames(Stat), function(x){return(paste0(substr(unlist(strsplit(x, "_"))[1], 1, 1), substr(unlist(strsplit(x, "_"))[2], 1, 3)))}))
Stat$ShortName[which(Stat$ShortName == "Blan")] <- "Blan3"
Stat$ShortName[which(Stat$ShortName == "BBra")] <- "Blan2"
Stat <- Stat[which(Stat$ShortName %in% c(Verteb, Amphi)),]

# Choose color
Stat$Color <- "grey50"
Stat$Color[which(Stat$ShortName == "Blan3")] <- "forestgreen"
Stat$Color[which(Stat$ShortName == "Blan2")] <- "red3"

###########################################################################
###########################################################################
# Plotting
###########################################################################
###########################################################################
pdf(paste(ResultsFolder, "/AssembliesComparison.pdf", sep=""), width=10, height=5)
par(mar=c(6,6,3,3),oma=c(1,1,1,1))
layout(matrix(c(1,2),nrow=1,ncol=2,byrow=T), widths=c(1,1), heights=c(1), TRUE)

Scatter_plot(Stat$prop1M*100, Stat$BUSCO_C, c(50,100), c(70,100), 5, Stat$ShortName, Stat$Color, "% of genome in >1M sequences", "BUSCO completeness score")
Scatter_plot(Stat$N50/1000000, Stat$BUSCO_C, c(0,150), c(70,100), 5, Stat$ShortName, Stat$Color, "N50 (M)", "BUSCO completeness score")
Scatter_plot(Stat$BUSCO_D, Stat$BUSCO_C, c(0,100), c(70,100), 5, Stat$ShortName, Stat$Color, "BUSCO duplication score", "BUSCO completeness score")
Scatter_plot(Stat$NumGaps1M, Stat$BUSCO_C, c(0,max(Stat$NumGaps1M)), c(70,100), 5, Stat$ShortName, Stat$Color, "# gaps in >1M sequences", "BUSCO completeness score")
Scatter_plot(Stat$prop1M*100, Stat$NumGaps1M, c(50,100), c(0,60000), 6, Stat$ShortName, Stat$Color, "% of genome in >1M sequences", "# gaps in >1M sequences")

dev.off()

