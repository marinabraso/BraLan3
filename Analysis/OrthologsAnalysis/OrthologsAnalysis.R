#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")

######################################################################
# Libraries & functions

`%out%` <- Negate(`%in%`)
script <- sub(".*=", "", commandArgs()[4])
source(paste(substr(script,1, nchar(script)-2), "_functions.R", sep=""))
source(paste(dirname(script), "/Amphioxus_Vertebrate_categories.R", sep=""))
source(paste(dirname(script), "/GeneExpression.R", sep=""))
source(paste(dirname(script), "/PositiveSelection.R", sep=""))
source(paste(dirname(script), "/Synteny.R", sep=""))
source(paste(dirname(script), "/GeneFunction.R", sep=""))
library(qvalue)
library(ghibli)
library(alluvial)
library(viridis)

######################################################################
# Files & folders

RDataFolder <- "Plots"
ResultsFolder <- "Plots/OrthologsAnalysis"
GeneDataFile <- paste(RDataFolder, "/GeneDataProcessed.txt", sep ="")
OGCountsFile <- paste(RDataFolder, "/OGCountsProcessed.txt", sep ="")
OGNamesFileAV <- "Results/FindOrthologs/AmphVerteb_broccoli/dir_step3/table_OGs_protein_names.txt"
RNASeqMetadataFile <- paste(RDataFolder, "/RNASeqMetadataProcessed.txt", sep ="")
GOMFFile <- paste(RDataFolder, "/GeneOntologyMainMF_FromHsap2BlanProcessed.txt", sep ="")
ChrLengthFile <- "Results/AssemblyStatistics/Branchiostoma_lanceolatum.BraLan3/Branchiostoma_lanceolatum.BraLan3_lengths.txt"
ProteomesFolder <- "Results/FilteringGeneSets/Proteomes"

######################################################################
# General parameters
QvalT <- 0.1 # q-value threshold for positive selection
PQvalThreshold <- 0.01 # for hypergeometric test
TandemMaxDist <- 1000 

Species <- c("Blan", "Bflo", "Bbel", "Drer", "Ggal", "Mmus", "Hsap")
GenesBasenames <- c("BLA", "BFLO", "BBEL", "ENSDAR", "ENSGAL", "ENSMUS", "ENS")
SpeciesProteomes <- c("Branchiostoma_lanceolatum.BraLan3.fa", "Branchiostoma_floridae.Bfl_VNyyK.fa", "Branchiostoma_belcheri.Haploidv18h27.fa", "Danio_rerio.GRCz11.fa", "Gallus_gallus.GRCg6a.fa", "Mus_musculus.GRCm39.fa", "Homo_sapiens.GRCh38.fa")

Tissues <- c("cirri", "gills", "epidermis", "gut", "hepatic.diverticulum", "muscle", "neural.tube", "female.gonads", "male.gonads")
colfunc <- colorRampPalette(c("forestgreen", "gold", "darkorange", "firebrick", "darkslateblue", "deepskyblue3"))
TissueColors <- colfunc(length(Tissues))

EmbAges <- c("egg", "32cells", "Blastula", "7h", "8h", "10h", "11h", "15h", "18h", "21h", "24h", "27h", "36h", "50h", "60h", "Pre.metamorphic.larvae")
#EmbAges <- c("32cells.Blastula", "7h.8h", "10h.11h", "15h.18h", "21h.24h", "27h.36h", "50h.60h")
colfunc <- colorRampPalette(c("forestgreen", "gold", "darkorange", "firebrick", "darkslateblue", "deepskyblue3"))
EmbAgesColors <- colfunc(length(EmbAges))

VTypes <- c("Missing", "SingleCopy", "Ohnolog", "Duplicated")
BlanTypes <- c("SingleCopy", "Duplicated")
VTypeColors <- c("deepskyblue2", "gold", "olivedrab2", "olivedrab4")
VTypeColors <- c("#481567FF", "#FDE725FF", "#1F968BFF", "#39568CFF")
VTypeColors <- ghibli_palettes$MarnieMedium2[c(7,2,4,6)]
BaseColor <- ghibli_palettes$MarnieMedium2[5]
BlanType.pch <- c(16,18)
BlanType.lty <- c(1,2)

# Big molecular function terms
MFunctGO <- c("GO:0016887", "GO:0016209", "GO:0003824", "GO:0005198", "GO:0140110", "GO:0140223", "GO:0005215", "GO:0038024", "GO:0004857", "GO:0008047", "GO:0019207", "GO:0099106", "GO:0098772", "GO:0030546", "GO:0019888", "GO:0044092", "GO:0044093", "GO:0065009", "GO:0045182", "GO:0005488", "GO:0060090", "GO:0140104", "GO:0060089")
MFunctNames <- c("ATPase", "AntioxidantActivity", "CatalyticActivity", "StructuralMoleculeActivity", "TranscriptionRegulatorActivity", "GeneralTranscriptionInitiationFactorActivity", "TransporterActivity", "CargoReceptorActivity", "EnzymeInhibitorActivity", "EnzymeActivatorActivity", "KinaseRegulatorActivity", "IonChannelRegulatorActivity", "MolecularFunctionRegulator", "SignalingReceptorActivatorActivity", "ProteinPhosphataseRegulatorActivity", "NegativeRegulationOfMolecularFunction", "PositiveRegulationOfMolecularFunction", "RegulationOfMolecularFunction", "TranslationRegulatorActivity", "Binding", "MolecularAdaptorActivity", "MolecularCarrierActivity", "MolecularTransducerActivity")
GeneralMFNames <- c("Enzymes", "StructuralProteins", "TranscriptionFactors", "Transporters", "Regulators", "Binding", "Carriers", "Signaling")
GeneralMFNum <- c(1, 1, 1, 2, 3, 3, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 7, 8)
GeneralMFMidp <- c(2, 4, 5.5, 7.5, 14, 20.5, 22, 23)
GeneralMFColors <- sample(c(ghibli_palettes$MononokeMedium, "darkred"))
GeneralMFColors <- viridis(8)


###########################################################################
###########################################################################
# Read data

GeneData <- read.table(GeneDataFile, h=T, sep = "\t", row.names=1)
colnames(GeneData) <- sub("^X", "", colnames(GeneData))

OGData.AV <- read.table(OGCountsFile, h=T, sep = "\t", row.names=1)
OGData.AV$DupType <- rep(NA, length(OGData.AV[,1]))
OGData.AV$DupType[which(OGData.AV$NumChrs==1 & OGData.AV$BlanType=="Duplicated")] <- rep("Intra",sum(OGData.AV$NumChrs==1 & OGData.AV$BlanType=="Duplicated"))
OGData.AV$DupType[which(OGData.AV$NumChrs>1 & OGData.AV$BlanType=="Duplicated")] <- rep("Inter",sum(OGData.AV$NumChrs>1 & OGData.AV$BlanType=="Duplicated"))
OGData.AV$DupType[which(OGData.AV$NumChrs==1 & OGData.AV$MaxDist<=TandemMaxDist & OGData.AV$BlanType=="Duplicated")] <- rep("Tandem",sum(OGData.AV$NumChrs==1 & OGData.AV$MaxDist<=TandemMaxDist & OGData.AV$BlanType=="Duplicated"))
table(OGData.AV$DupType)
table(OGData.AV$BlanType)
head(OGData.AV)

MetaRNA <- read.table(RNASeqMetadataFile, h=T, sep = "\t", row.names=1)

ChrLengths <- read.table(ChrLengthFile, h=F, sep = "\t", row.names=1)
colnames(ChrLengths) <- c("Length", "Chr", "AccLength")
rownames(ChrLengths) <- sub("^>", "", rownames(ChrLengths))

GOMFData <- read.table(GOMFFile, h=T, sep = "\t")
print(head(GOMFData))

system_out <- system(paste0("cat ", OGNamesFileAV, " | tail -n +2 | awk '{for(i=2;i<=NF;i++){print $1\"\t\"$i}}' | sort -u | awk '{if(a[$2]){a[$2]=\"NA\"}else{a[$2]=$1}}END{for(i in a){print a[i]\"\t\"i}}'"), intern=T)
OG2gene.AV <- read.table(text=system_out, h=F, sep = "\t")
colnames(OG2gene.AV) <- c("OG", "Gene")
OG2gene.AV$Species <- rep(NA, length(OG2gene.AV[,1]))
for(sp in c(length(Species):1)){
	print(Species[sp])
	listGenesSp <- grep(GenesBasenames[sp], OG2gene.AV$Gene)
	OG2gene.AV$Species[listGenesSp] <- rep(Species[sp], length(listGenesSp))
}
table(OG2gene.AV$Species)

NumGenesSpecies <- rep(0, length(Species))
for(sp in c(1:length(Species))){
	system_out <- system(paste0("cat ", ProteomesFolder, "/", SpeciesProteomes[sp], " | grep '>' | sort | uniq | wc -l"), intern=T)
	NumGenesSpecies[sp] <- read.table(text=system_out, h=F, sep = "\t")
}
NumGenesSpecies <- unlist(NumGenesSpecies)

###########################################################################
###########################################################################
# Amphioxus - vertebrate shared OG numbers
paste("Shared V-A", length(OGData.AV[which(OGData.AV$SumVerteb>0 & OGData.AV$SumAmphi>0),1]))
paste("V", length(OGData.AV[which(OGData.AV$SumVerteb>0),1]))
paste("A", length(OGData.AV[which(OGData.AV$SumAmphi>0),1]))
paste("Shared/V*100", length(OGData.AV[which(OGData.AV$SumVerteb>0 & OGData.AV$SumAmphi>0),1])/length(OGData.AV[which(OGData.AV$SumVerteb>0),1])*100)
paste("Shared/A*100", length(OGData.AV[which(OGData.AV$SumVerteb>0 & OGData.AV$SumAmphi>0),1])/length(OGData.AV[which(OGData.AV$SumAmphi>0),1])*100)

# Amphioxus specfic OG numbers
OGData.A <- OGData.AV[which(OGData.AV$SumVerteb==0 & OGData.AV$SumAmphi>0),]
paste("A specific", length(OGData.A[,1]))
paste("Shared all A", length(OGData.A[which(OGData.A$Blan>0 & OGData.A$Bflo>0 & OGData.A$Bbel>0),1]))
paste("Shared all A/A specific*100", length(OGData.A[which(OGData.A$Blan>0 & OGData.A$Bflo>0 & OGData.A$Bbel>0),1])/length(OGData.A[,1])*100)
paste("blan-bbel specific", length(OGData.A[which(OGData.A$Blan>0 & OGData.A$Bflo==0 & OGData.A$Bbel>0),1]))
paste("blan-bbel specific/A specific*100", length(OGData.A[which(OGData.A$Blan>0 & OGData.A$Bflo==0 & OGData.A$Bbel>0),1])/length(OGData.A[,1])*100)
paste("blan-bflo specific", length(OGData.A[which(OGData.A$Blan>0 & OGData.A$Bflo>0 & OGData.A$Bbel==0),1]))
paste("blan-bflo specific/A specific*100", length(OGData.A[which(OGData.A$Blan>0 & OGData.A$Bflo>0 & OGData.A$Bbel==0),1])/length(OGData.A[,1])*100)

for(sp in c(1:length(Species))){
	OG2gene.sp <- OG2gene.AV[which(OG2gene.AV$Species==Species[sp]),]
	toprint <- c(1:7)
	toprint[1] <- Species[sp]
	toprint[2] <- sum(OGData.AV[,Species[sp]]>1)
	toprint[3] <- sum(OGData.AV[,Species[sp]]>1)/sum(OGData.AV[,Species[sp]]>0)*100
	toprint[4] <- length(unique(OG2gene.sp$Gene[which(OG2gene.sp$OG %in% rownames(OGData.AV[which(OGData.AV[,Species[sp]]>1),]))]))
	toprint[5] <- length(unique(OG2gene.sp$Gene[which(OG2gene.sp$OG %in% rownames(OGData.AV[which(OGData.AV[,Species[sp]]>1),]))]))/NumGenesSpecies[sp]*100
	toprint[6] <- NumGenesSpecies[sp]-length(unique(OG2gene.sp$Gene))
	toprint[7] <- (NumGenesSpecies[sp]-length(unique(OG2gene.sp$Gene)))/NumGenesSpecies[sp]*100
	toprint[8] <- length(unique(OG2gene.sp$Gene[which(OG2gene.sp$OG %in% rownames(OGData.AV[which(OGData.AV[,Species[sp]]>0 & OGData.AV$Sum==OGData.AV[,Species[sp]]),]))]))
	toprint[9] <- length(unique(OG2gene.sp$Gene[which(OG2gene.sp$OG %in% rownames(OGData.AV[which(OGData.AV[,Species[sp]]>0 & OGData.AV$Sum==OGData.AV[,Species[sp]]),]))]))/NumGenesSpecies[sp]*100
	print(paste(toprint, collapse=" "))
}
length(OGData.AV$Blan[which(OGData.AV$Sum==OGData.AV$Blan)])
length(OGData.AV$Blan[which(OGData.AV$Sum==OGData.AV$Blan)])

Amphioxus_Vertebrate_categories(OGData.AV)

#GeneExpression(GeneData)

#PositiveSelection(GeneData, OGData.AV)

#Synteny(GeneData, OGData.AV, ChrLengths)

#GeneFunction(GeneData, OGData.AV, GOMFData)

print("Nodal (Z linked?)")
print(GeneData[which(GeneData$Gene=="BLAG15000766"),])


