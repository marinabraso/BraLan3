#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")

######################################################################
# Libraries & functions
pwd <- system("pwd", intern=TRUE)
script <- sub(".*=", "", commandArgs()[4])
source(paste0(dirname(script), "/", substr(basename(script),1, nchar(basename(script))-2), "_functions.r", sep=""))
# source("./Analysis/ComparativeGeneExpression/ComparativeGeneExpression_functions.r")

######################################################################
# Files & folders
RDataFolder <- "Plots"
ResultsFolder <- "Plots/ComparativeGeneExpression"
OrganOrthBgeeFile <- paste0(pwd, "/Data/AnatomicalOntology/Julien_AnatOrht_Bgee15.tsv")
CountsFileAV <- paste0(pwd, "/Results/FindOrthologs/AmphVerteb_broccoli/dir_step3/table_OGs_protein_counts.txt")
NamesFileAV <- paste0(pwd, "/Results/FindOrthologs/AmphVerteb_broccoli/dir_step3/table_OGs_protein_names.txt")
BlanGeneDataFile <- paste0(pwd, "/", RDataFolder, "/GeneDataProcessed.txt")
OhnologsSFile <- paste0(pwd, "/Results/FindOrthologs/OGOnhologs/OG_w_Onhologs.Strict.DR_GG_MM_HS.txt")
DrerGeneDataFile <- paste0(pwd, "/", ResultsFolder, "/DrerGeneData.txt")
DrerBgeeDataFile <- paste0(pwd, "/", ResultsFolder, "/DrerBgeeData.txt")

######################################################################
# General parameters
TPM.threshold <- 1
Species1 <- "Blan"
Species2 <- "Drer"

ShortSpeciesNames <- c("Drer", "Ggal", "Hsap", "Mmus", "Blan")
SpeciesNames <- c("Danio_rerio", "Gallus_gallus", "Homo_sapiens", "Mus_musculus", "Branchiostoma_lanceolatum")
SName1 <- SpeciesNames[which(ShortSpeciesNames==Species1)]
SName2 <- SpeciesNames[which(ShortSpeciesNames==Species2)]

BlanTissues <- c("cirri", "gills", "epidermis", "gut", "hepatic.diverticulum", "muscle", "neural.tube", "female.gonads", "male.gonads", "egg", "32cells", "Blastula", "MeanEmbr")
DrerTissues <- c("blastula", "embryo", "head", "tail", "granulocyte", "brain", "heart", "muscle tissue", "mesonephros", "intestine", "pharyngeal gill", "testis", "ovary", "bone element", "liver", "mature ovarian follicle", "zone of skin", "swim bladder", "head kidney", "spleen")
DrerTissuesIDs <- c("UBERON:0000307", "UBERON:0000922", "UBERON:0000033", "UBERON:0002415", "CL:0000094", "UBERON:0000955", "UBERON:0000948", "UBERON:0002385", "UBERON:0000080", "UBERON:0000160", "UBERON:0000206", "UBERON:0000473", "UBERON:0000992", "UBERON:0001474", "UBERON:0002107", "UBERON:0003982", "UBERON:0000014", "UBERON:0006860", "UBERON:0007132", "UBERON:0002106")

BlanMatchingTissues <- c("Blastula", "MeanEmbr", "male.gonads", "female.gonads", "muscle", "neural.tube", "gut", "gills", "hepatic.diverticulum", "epidermis")
DrerMatchingTissues <- c("blastula", "embryo", "testis", "ovary", "muscle tissue", "brain", "intestine", "pharyngeal gill", "liver", "zone of skin")
NamesMatchingTissues <- c("blastula", "embryo", "testis", "ovary", "muscle", "neural", "digestive", "gills", "hepatic", "skin")
colfunc <- colorRampPalette(c("forestgreen", "gold", "darkorange", "firebrick", "darkslateblue", "deepskyblue3"))
TissueColors <- colfunc(length(NamesMatchingTissues))

MatchingTissues <- as.data.frame(cbind(BlanMatchingTissues, DrerMatchingTissues, NamesMatchingTissues, TissueColors))
colnames(MatchingTissues) <- c("Blan", "Drer", "Name","Color")
MatchingTissues <- MatchingTissues[MatchingTissues$Name!="blastula" & MatchingTissues$Name!="ovary" & MatchingTissues$Name!="skin",]

###########################################################################
###########################################################################
# Read data

# List of onhologs
Ohnologs <- read.table(OhnologsSFile, h=F)[,1]

# Read Drer gene information previously processed or precess it and write itinto a file
DrerGeneData <- prepare_DrerGeneData(DrerGeneDataFile, DrerBgeeDataFile, SName2, MatchingTissues, ResultsFolder, pwd)
print(head(DrerGeneData))

# Read Blan gene data file 
BlanGeneData <- read.table(BlanGeneDataFile, h=T, sep = "\t", row.names=1)
colnames(BlanGeneData) <- sub("^X", "", colnames(BlanGeneData))

# Read orthologs counts information
CountsAV <- prepare_CountsAV(CountsFileAV, Ohnologs)
print(head(CountsAV))

# Prepare Blan - Drer gene pairs
GenePairs <- prepare_GenePairs(NamesFileAV, BlanGeneData, DrerGeneData, MatchingTissues, TPM.threshold)
print(head(GenePairs))
GenePairs$DifTissues <- lapply(c(1:length(GenePairs[,1])), difference_union_of_patterns_per_tissue_gp, GenePairs, MatchingTissues, TPM.threshold)

# Label genepairs that are ohnologs
GenePairs$Ohnologs <- rep(FALSE, length(GenePairs[,1]))
GenePairs$Ohnologs[which(GenePairs$OG %in% Ohnologs)] <- rep(TRUE, length(GenePairs[which(GenePairs$OG %in% Ohnologs),1]))

# Calc difference in the sum of the union of patterns un each species per OG
CountsAV$SumUnionDom1 <- unlist(lapply(CountsAV$OG, sum_union_of_patterns, 1, GenePairs, MatchingTissues, TPM.threshold))
CountsAV$SumUnionDom2 <- unlist(lapply(CountsAV$OG, sum_union_of_patterns, 2, GenePairs, MatchingTissues, TPM.threshold))
CountsAV$DiffUnionDom <- CountsAV$SumUnionDom1-CountsAV$SumUnionDom2

# Absolute difference between the union of patterns un each species per OG
CountsAV$DiffAbsUnion <- unlist(lapply(CountsAV$OG, absolute_difference_union_of_patterns, GenePairs, MatchingTissues, TPM.threshold))
CountsAV$DifTissues <- lapply(CountsAV$OG, difference_union_of_patterns_per_tissue_c, GenePairs, MatchingTissues, TPM.threshold)
print(head(CountsAV$DifTissues))
print(summary(CountsAV$DiffUnionDom))
print(summary(CountsAV$DiffAbsUnion))

# Union of paterns per tissue ans species (sum of TPM values) 
for(i in c(1:length(MatchingTissues[,1]))){
	CountsAV[, paste0(MatchingTissues$Name[i], 1)] <- unlist(lapply(CountsAV$OG, sum_TPM_per_species_tissue, 1, GenePairs, MatchingTissues$Name[i]))
	CountsAV[, paste0(MatchingTissues$Name[i], 2)] <- unlist(lapply(CountsAV$OG, sum_TPM_per_species_tissue, 2, GenePairs, MatchingTissues$Name[i]))
}

# Subdivisions of GenePairs & Counts (1 to 1 orthologs, 1 to many, many to 1)
GenePairs.121 <- GenePairs[which(GenePairs$OG %in% CountsAV$OG[which(CountsAV$Blan==1 & CountsAV$Drer==1)]),]
GenePairs.12m <- GenePairs[which(GenePairs$OG %in% CountsAV$OG[which(CountsAV$Blan==1 & CountsAV$Drer>=2)]),]
GenePairs.12mo <- GenePairs[which(GenePairs$OG %in% CountsAV$OG[which(CountsAV$Blan==1 & CountsAV$Drer>=2 & CountsAV$Ohnologs==TRUE)]),]
GenePairs.12ms <- GenePairs[which(GenePairs$OG %in% CountsAV$OG[which(CountsAV$Blan==1 & CountsAV$Drer>=2 & CountsAV$Ohnologs==FALSE)]),]
GenePairs.m21 <- GenePairs[which(GenePairs$OG %in% CountsAV$OG[which(CountsAV$Blan>=2 & CountsAV$Drer==1)]),]
CountsAV.121 <- CountsAV[which(CountsAV$Blan==1 & CountsAV$Drer==1),]
CountsAV.12m <- CountsAV[which(CountsAV$Blan==1 & CountsAV$Drer>=2),]
CountsAV.12mo <- CountsAV[which(CountsAV$Blan==1 & CountsAV$Drer>=2 & CountsAV$Ohnologs==TRUE),]
CountsAV.12ms <- CountsAV[which(CountsAV$Blan==1 & CountsAV$Drer>=2 & CountsAV$Ohnologs==FALSE),]
CountsAV.m21 <- CountsAV[which(CountsAV$Blan>=2 & CountsAV$Drer==1),]

print(table(CountsAV$Ohnologs))
print(head(CountsAV))

# Calc spearman correlations
for(i in c(1:length(MatchingTissues[,1]))){
	MatchingTissues$Corr[i] <- get_spearman_corr(GenePairs, MatchingTissues$Name[i])
	MatchingTissues$Corr.121[i] <- get_spearman_corr(GenePairs.121, MatchingTissues$Name[i])
	MatchingTissues$Corr.m21[i] <- get_spearman_corr(GenePairs.m21, MatchingTissues$Name[i])
	MatchingTissues$Corr.12m[i] <- get_spearman_corr(GenePairs.12m, MatchingTissues$Name[i])
	MatchingTissues$Corr.12mo[i] <- get_spearman_corr(GenePairs.12mo, MatchingTissues$Name[i])
	MatchingTissues$Corr.12ms[i] <- get_spearman_corr(GenePairs.12ms, MatchingTissues$Name[i])
	MatchingTissues$CorrPA[i] <- get_spearman_corr_presence_absence(GenePairs, MatchingTissues$Name[i], TPM.threshold)
	MatchingTissues$CorrPA.121[i] <- get_spearman_corr_presence_absence(GenePairs.121, MatchingTissues$Name[i], TPM.threshold)
	MatchingTissues$CorrPA.m21[i] <- get_spearman_corr_presence_absence(GenePairs.m21, MatchingTissues$Name[i], TPM.threshold)
	MatchingTissues$CorrPA.12m[i] <- get_spearman_corr_presence_absence(GenePairs.12m, MatchingTissues$Name[i], TPM.threshold)
	MatchingTissues$CorrPA.12mo[i] <- get_spearman_corr_presence_absence(GenePairs.12mo, MatchingTissues$Name[i], TPM.threshold)
	MatchingTissues$CorrPA.12ms[i] <- get_spearman_corr_presence_absence(GenePairs.12ms, MatchingTissues$Name[i], TPM.threshold)
	MatchingTissues$Corr.um21[i] <- get_spearman_corr(CountsAV.m21, MatchingTissues$Name[i])
	MatchingTissues$Corr.u12m[i] <- get_spearman_corr(CountsAV.12m, MatchingTissues$Name[i])
	MatchingTissues$Corr.u12mo[i] <- get_spearman_corr(CountsAV.12mo, MatchingTissues$Name[i])
	MatchingTissues$Corr.u12ms[i] <- get_spearman_corr(CountsAV.12ms, MatchingTissues$Name[i])
	MatchingTissues$CorrPA.um21[i] <- get_spearman_corr_presence_absence(CountsAV.m21, MatchingTissues$Name[i], TPM.threshold)
	MatchingTissues$CorrPA.u12m[i] <- get_spearman_corr_presence_absence(CountsAV.12m, MatchingTissues$Name[i], TPM.threshold)
	MatchingTissues$CorrPA.u12mo[i] <- get_spearman_corr_presence_absence(CountsAV.12mo, MatchingTissues$Name[i], TPM.threshold)
	MatchingTissues$CorrPA.u12ms[i] <- get_spearman_corr_presence_absence(CountsAV.12ms, MatchingTissues$Name[i], TPM.threshold)
}


pdf(paste0(ResultsFolder, "/ComparativeGeneExpression.pdf"), width=20, height=10)
par(mar=c(5,5,2,2),oma=c(1,1,1,1), yaxs='i', xaxs='i')
layout(matrix(c(1,2),nrow=1,ncol=2,byrow=T), widths=c(1), heights=c(1), TRUE)

breaks <- seq(-.5,length(MatchingTissues[,1])+.5,1)
#hist(GenePairs.121$SumDom1, breaks=breaks)
#hist(GenePairs.121$SumDom2, breaks=breaks)

layout(matrix(c(1,2,3,4,5,6,7,8),nrow=2,ncol=4,byrow=T), widths=c(1), heights=c(1), TRUE)

breaks <- seq(-.5,length(MatchingTissues[,1])+.5,1)
add_joined_relative_hist_line(breaks, c("black", "royalblue1", "royalblue4"), "many-to-1", "Absolute difference between patterns", GenePairs.121$DiffAbs, GenePairs.m21$DiffAbs, CountsAV.m21$DiffAbsUnion)
add_joined_relative_hist_line(breaks, c("black", "royalblue1", "royalblue4"), "1-to-many ohnologs", "Absolute difference between patterns", GenePairs.121$DiffAbs, GenePairs.12mo$DiffAbs, CountsAV.12mo$DiffAbsUnion)
add_joined_relative_hist_line(breaks, c("black", "royalblue1", "royalblue4"), "1-to-many small scale duplicates", "Absolute difference between patterns", GenePairs.121$DiffAbs, GenePairs.12ms$DiffAbs, CountsAV.12ms$DiffAbsUnion)
plot.new()

breaks <- seq(-length(MatchingTissues[,1])-.5,length(MatchingTissues[,1])+.5,1)
add_joined_relative_hist_line(breaks, c("black", "royalblue1", "royalblue4"), "many-to-1", "Blan domains - Drer domains", -GenePairs.121$DiffDom, -GenePairs.m21$DiffDom, -CountsAV.m21$DiffUnionDom)
add_joined_relative_hist_line(breaks, c("black", "royalblue1", "royalblue4"), "1-to-many ohnologs", "Blan domains - Drer domains", -GenePairs.121$DiffDom, -GenePairs.12mo$DiffDom, -CountsAV.12mo$DiffUnionDom)
add_joined_relative_hist_line(breaks, c("black", "royalblue1", "royalblue4"), "1-to-many small scale duplicates", "Blan domains - Drer domains", -GenePairs.121$DiffDom, -GenePairs.12ms$DiffDom, -CountsAV.12ms$DiffUnionDom)
plot.new()


plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,1), xlim=c(-.5,length(MatchingTissues[,1])+.5), col=NA)
mtext("Frequency", side = 2, line = 3, cex=1.5)
mtext("Differences between patterns", side = 1, line = 3, cex=1.5)
add_relative_freq_line(GenePairs.121$DiffAbs, breaks, "black")
#add_relative_freq_line(GenePairs.12m$DiffAbs, breaks, "royalblue1")
#add_relative_freq_line(CountsAV.12m$DiffAbsUnion, breaks, "royalblue3")
add_relative_freq_line(GenePairs.m21$DiffAbs, breaks, "orangered1")
add_relative_freq_line(CountsAV.m21$DiffAbsUnion, breaks, "orangered3")
add_relative_freq_line(GenePairs.12mo$DiffAbs, breaks, "olivedrab1")
add_relative_freq_line(CountsAV.12mo$DiffAbsUnion, breaks, "olivedrab3")
add_relative_freq_line(GenePairs.12ms$DiffAbs, breaks, "gold1")
add_relative_freq_line(CountsAV.12ms$DiffAbsUnion, breaks, "gold3")
axis(1, at = c(0:length(MatchingTissues[,1])), lwd.ticks=1, las=1, cex.axis=1)
axis(2, at = seq(0,1,.2), lwd.ticks=1, las=1, cex.axis=1)
legend("topright", c("Blan - Drer", "1-to-1", "1-to-many", "many-to-1", "1-to-many ohnologs", "1-to-many small scale duplicates", "Indvidual paterns", "Union of patterns"), pch=15, col=c("white", "black","royalblue1","orangered1", "olivedrab1", "gold1", "grey80", "grey40"), bty = "n", pt.cex=2, cex=1.5, xjust = 0, yjust = 0)
box()

breaks <- seq(-length(MatchingTissues[,1])-.5,length(MatchingTissues[,1])+.5,1)
plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,1), xlim=c(-length(MatchingTissues[,1])+.5,length(MatchingTissues[,1])+.5), col=NA)
mtext("Frequency", side = 2, line = 3, cex=1.5)
mtext("Differences between patterns", side = 1, line = 3, cex=1.5)
add_relative_freq_line(GenePairs.121$DiffDom, breaks, "black")
add_relative_freq_line(GenePairs.121$DiffDom, breaks, "black")
#add_relative_freq_line(GenePairs.12m$DiffDom, breaks, "royalblue1")
#add_relative_freq_line(CountsAV.12m$DiffUnionDom, breaks, "royalblue3")
add_relative_freq_line(GenePairs.m21$DiffDom, breaks, "orangered1")
add_relative_freq_line(CountsAV.m21$DiffUnionDom, breaks, "orangered3")
add_relative_freq_line(GenePairs.12mo$DiffDom, breaks, "olivedrab1")
add_relative_freq_line(CountsAV.12mo$DiffUnionDom, breaks, "olivedrab3")
add_relative_freq_line(GenePairs.12ms$DiffDom, breaks, "gold1")
add_relative_freq_line(CountsAV.12ms$DiffUnionDom, breaks, "gold3")
axis(1, at = c(0:length(MatchingTissues[,1])), lwd.ticks=1, las=1, cex.axis=1)
axis(2, at = seq(0,1,.2), lwd.ticks=1, las=1, cex.axis=1)
legend("topright", c("Blan - Drer", "1-to-1", "1-to-many", "many-to-1", "1-to-many ohnologs", "1-to-many small scale duplicates", "Indvidual paterns", "Union of patterns"), pch=15, col=c("white", "black","royalblue1","orangered1", "olivedrab1", "gold1", "grey80", "grey40"), bty = "n", pt.cex=2, cex=1.5, xjust = 0, yjust = 0)
box()

sep <- 0.1
layout(matrix(c(1,2,3,4),nrow=2,ncol=2,byrow=T), widths=c(2), heights=c(1), TRUE)
plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,1), xlim=c(0,length(MatchingTissues[,1])+1), col=NA)
mtext("Spearman correlation", side = 2, line = 5, cex=1.5)
points(c(1:length(MatchingTissues[,1]))-sep*1.5, MatchingTissues$Corr.121, col="black", pch=16, cex=2)
points(c(1:length(MatchingTissues[,1]))-sep*0.5, MatchingTissues$Corr.m21, col="orangered1", pch=16, cex=2)
points(c(1:length(MatchingTissues[,1]))-sep*0.5, MatchingTissues$Corr.um21, col="orangered1", pch=18, cex=2)
points(c(1:length(MatchingTissues[,1]))+sep*0.5, MatchingTissues$Corr.12mo, col="royalblue1", pch=16, cex=2)
points(c(1:length(MatchingTissues[,1]))+sep*0.5, MatchingTissues$Corr.u12mo, col="royalblue1", pch=18, cex=2)
points(c(1:length(MatchingTissues[,1]))+sep*1.5, MatchingTissues$Corr.12ms, col="gold1", pch=16, cex=2)
points(c(1:length(MatchingTissues[,1]))+sep*1.5, MatchingTissues$Corr.u12ms, col="gold1", pch=18, cex=2)
axis(1, at = c(1:length(MatchingTissues[,1])), labels=MatchingTissues$Name, lwd.ticks=1, las=1, cex.axis=1.2)
axis(2, at = seq(0,1,.2), lwd.ticks=1, las=1, cex.axis=1.2)

layout(matrix(c(1,2,3,4),nrow=2,ncol=2,byrow=T), widths=c(2), heights=c(1), TRUE)
plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,1), xlim=c(0,length(MatchingTissues[,1])+1), col=NA)
mtext("Spearman correlation", side = 2, line = 5, cex=1.5)
points(c(1:length(MatchingTissues[,1]))-sep*1.5, MatchingTissues$CorrPA.121, col="black", pch=16, cex=2)
points(c(1:length(MatchingTissues[,1]))-sep*0.5, MatchingTissues$CorrPA.m21, col="orangered1", pch=16, cex=2)
points(c(1:length(MatchingTissues[,1]))-sep*0.5, MatchingTissues$CorrPA.12mo, col="royalblue1", pch=16, cex=2)
points(c(1:length(MatchingTissues[,1]))+sep*0.5, MatchingTissues$CorrPA.12ms, col="gold1", pch=16, cex=2)
points(c(1:length(MatchingTissues[,1]))+sep*0.5, MatchingTissues$CorrPA.um21, col="orangered1", pch=18, cex=2)
points(c(1:length(MatchingTissues[,1]))+sep*1.5, MatchingTissues$CorrPA.u12mo, col="royalblue1", pch=18, cex=2)
points(c(1:length(MatchingTissues[,1]))+sep*1.5, MatchingTissues$CorrPA.u12ms, col="gold1", pch=18, cex=2)
axis(1, at = c(1:length(MatchingTissues[,1])), labels=MatchingTissues$Name, lwd.ticks=1, las=1, cex.axis=1.2)
axis(2, at = seq(0,1,.2), lwd.ticks=1, las=1, cex.axis=1.2)


layout(matrix(c(1,2,3,4),nrow=2,ncol=2,byrow=T), widths=c(2), heights=c(1), TRUE)
plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,.4), xlim=c(0,length(MatchingTissues[,1])+1), col=NA)
mtext("Number of differences", side = 2, line = 5, cex=1.5)
for(t in c(1:length(MatchingTissues[,1]))){
	points(t, length(grep(MatchingTissues$Name[t], unlist(GenePairs.121$DifTissues)))/length(na.omit(unlist(GenePairs.121$DifTissues))), col="black", pch=16, cex=2)
	points(t, length(grep(MatchingTissues$Name[t], unlist(GenePairs.m21$DifTissues)))/length(na.omit(unlist(GenePairs.m21$DifTissues))), col="orangered1", pch=16, cex=2)
	points(t, length(grep(MatchingTissues$Name[t], unlist(GenePairs.12mo$DifTissues)))/length(na.omit(unlist(GenePairs.12mo$DifTissues))), col="royalblue1", pch=16, cex=2)
	points(t, length(grep(MatchingTissues$Name[t], unlist(GenePairs.12ms$DifTissues)))/length(na.omit(unlist(GenePairs.12ms$DifTissues))), col="gold1", pch=16, cex=2)
	points(t, length(grep(MatchingTissues$Name[t], unlist(CountsAV.121$DifTissues)))/length(na.omit(unlist(CountsAV.121$DifTissues))), col="black", pch=18, cex=2)
	points(t, length(grep(MatchingTissues$Name[t], unlist(CountsAV.m21$DifTissues)))/length(na.omit(unlist(CountsAV.m21$DifTissues))), col="orangered1", pch=18, cex=2)
	points(t, length(grep(MatchingTissues$Name[t], unlist(CountsAV.12mo$DifTissues)))/length(na.omit(unlist(CountsAV.12mo$DifTissues))), col="royalblue1", pch=18, cex=2)
	points(t, length(grep(MatchingTissues$Name[t], unlist(CountsAV.12ms$DifTissues)))/length(na.omit(unlist(CountsAV.12ms$DifTissues))), col="gold1", pch=18, cex=2)
}
axis(1, at = c(1:length(MatchingTissues[,1])), labels=MatchingTissues$Name, lwd.ticks=1, las=1, cex.axis=1.2)
axis(2, at = seq(0,1,.2), lwd.ticks=1, las=1, cex.axis=1.2)

dev.off()






