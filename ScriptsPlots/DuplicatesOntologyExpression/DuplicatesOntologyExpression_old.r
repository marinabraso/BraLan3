#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")

######################################################################
# Libraries & functions

#script <- sub(".*=", "", commandArgs()[4])
#source(paste(substr(script,1, nchar(script)-2), "_functions.r", sep=""))
source("ScriptsPlots/DuplicatesOntologyExpression/DuplicatesOntologyExpression_functions.r")
# source("ScriptsPlots/DuplicatesOntologyExpression/DuplicatesOntologyExpression.r")
library(ghibli)
library(viridis)
library(car)
library(MASS)

######################################################################
# Files & folders

ResultsFolder <- "Plots/DuplicatesOntologyExpression"
CountsFile <- "Results/FindOrthologs/AmphVerteb_broccoli/dir_step3/table_OGs_protein_counts.txt"
NamesFile <- "Results/FindOrthologs/AmphVerteb_broccoli/dir_step3/table_OGs_protein_names.txt"
OhnologsSFile <- "Results/FindOrthologs/OGOnhologs/OG_w_Onhologs.Strict.DR_GG_MM_HS.txt"
Ohnolog2ROG <- "Results/OhnologListing/2R_Strict_OG_pairs.txt"
Ohnolog3ROG <- "Results/OhnologListing/3R_Strict_OG_pairs.txt"
GOMFFile <- paste("Plots/GeneOntologyMainMF_FromHsap2BlanProcessed.txt", sep ="")
GOlistsFolder <- "Results/GeneOntology"
MainGOlistFile <- "Data/GeneOntology/MainGOMF_supclass.list"
GOlistFile <- "Results/GeneOntology/GOlist_MF.list"
GTFFolder <- "Results/FilteringGeneSets/FilteredGTFs"
BlanGeneDataFile <- paste0(ResultsFolder, "/BlanGeneData.txt")
DrerGeneDataFile <- paste0(ResultsFolder, "/DrerGeneData.txt")
DrerBgeeDataFile <- paste0(ResultsFolder, "/DrerBgeeData.txt")
RNASeqMetadataFile <- "Metadata/Marletaz2018_RNAseq_SRA.txt"
TPMFile <- "Results/GeneExpression/Gene_TPM_kallisto_tximport.tab"
ProteomesFolder <- "Results/FilteringGeneSets/FilteredProteomes"

######################################################################
# General parameters
PQvalThreshold <- 0.01 # for hypergeometric test
minGOsizeHsap <- 50
minGOsizeBlan <- 10
MaxTandemDist <- 10000
TPM.threshold <- 1


Species <- c("Blan", "Bflo", "Bbel", "Drer", "Ggal", "Mmus", "Hsap")
SpType <- c("AmphType", "AmphType", "AmphType", "VertType", "VertType", "VertType", "VertType")
Ampioxus <- c("Blan", "Bflo", "Bbel")
Vertebrates <- c("Drer", "Ggal", "Mmus", "Hsap")
GenesBasenames <- c("BLA", "BFLO", "BBEL", "ENSDAR", "ENSGAL", "ENSMUS", "ENS")
SpeciesLongNames <- c("Branchiostoma_lanceolatum.BraLan3", "Branchiostoma_floridae.Bfl_VNyyK", "Branchiostoma_belcheri.Haploidv18h27", "Danio_rerio.GRCz11", "Gallus_gallus.GRCg6a", "Mus_musculus.GRCm39", "Homo_sapiens.GRCh38")
BaseColor <- ghibli_palettes$MarnieMedium2[5]
VertTypes <- c("Single-copy", "Ohnologs", "Small-scale\nduplicates", "Missing")
VertTypes.abr <- c("SC", "O", "D", "M")
VertTypes.col <- c(ghibli_palettes$MarnieMedium2[c(2,4,7)], "royalblue2")
BlanTypes <- c("Single-copy", "Small-scale\nduplicates", "Missing")
BlanTypes.pch <- c(14, 16, 18)
BlanTypes.lty <- c(.5, 1, 2)

# Big molecular function terms
MFunctGO <- c("GO:0016887", "GO:0016209", "GO:0003824", "GO:0005198", "GO:0140110", "GO:0140223", "GO:0005215", "GO:0038024", "GO:0004857", "GO:0008047", "GO:0019207", "GO:0099106", "GO:0098772", "GO:0030546", "GO:0019888", "GO:0044092", "GO:0044093", "GO:0065009", "GO:0045182", "GO:0005488", "GO:0060090", "GO:0140104", "GO:0060089")
MFunctNames <- c("ATPase", "AntioxidantActivity", "CatalyticActivity", "StructuralMoleculeActivity", "TranscriptionRegulatorActivity", "GeneralTranscriptionInitiationFactorActivity", "TransporterActivity", "CargoReceptorActivity", "EnzymeInhibitorActivity", "EnzymeActivatorActivity", "KinaseRegulatorActivity", "IonChannelRegulatorActivity", "MolecularFunctionRegulator", "SignalingReceptorActivatorActivity", "ProteinPhosphataseRegulatorActivity", "NegativeRegulationOfMolecularFunction", "PositiveRegulationOfMolecularFunction", "RegulationOfMolecularFunction", "TranslationRegulatorActivity", "Binding", "MolecularAdaptorActivity", "MolecularCarrierActivity", "MolecularTransducerActivity")
GeneralMFNames <- c("Enzymes", "StructuralProteins", "TranscriptionFactors", "Transporters", "Regulators", "Binding", "Carriers", "Signaling")
GeneralMFNum <- c(1, 1, 1, 2, 3, 3, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 7, 8)
GeneralMFMidp <- c(2, 4, 5.5, 7.5, 14, 20.5, 22, 23)
GeneralMFColors <- viridis(8)

# Gene expression analysis
Species1 <- "Blan"
Species2 <- "Drer"
ShortSpeciesNames <- c("Drer", "Ggal", "Hsap", "Mmus", "Blan")
SpeciesNames <- c("Danio_rerio", "Gallus_gallus", "Homo_sapiens", "Mus_musculus", "Branchiostoma_lanceolatum")
SpeciesGRefs <- c("Danio_rerio.GRCz11", "Gallus_gallus.GRCg6a", "Homo_sapiens.GRCh38", "Mus_musculus.GRCm39", "Branchiostoma_lanceolatum.BraLan3", "Branchiostoma_floridae.Bfl_VNyyK")

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
OhnOG2R <- read.table(Ohnolog2ROG, h=F)[,1]
OhnOG3R <- read.table(Ohnolog3ROG, h=F)[,1]

# Orthologous groups counts
Counts <- read.table(CountsFile, h=F, sep = "\t", row.names=1)
system_out <- system(paste("head -1 ", CountsFile, " | cut -f2,3,4,5,6,7,8,9,10,11,12 | sed 's/\\([A-Z]\\)[a-z]\\+_\\([a-z][a-z][a-z]\\)[A-Za-z0-9\\._]\\+/\\1\\2/g'"), intern=T)
colnames(Counts) <- as.character(unlist(lapply(read.table(text=system_out, h=F, sep = "\t")[1,], as.character)))
Counts$OG <- rownames(Counts)
#Counts$VertType <- c("Missing", "Single-copy", "Small-scale\nduplicates")[unlist(apply(Counts[,Vertebrates], 1, function(x){r <- max(x); if(r>1){r=2}; return(r)}))+1]
#Counts$VertType[rownames(Counts) %in% OhnOG2R] <- "Ohnologs"
#table(Counts$VertType)
Counts$VertType <- c("Missing", "Single-copy", "Small-scale\nduplicates")[unlist(apply(Counts[,c("OG",Vertebrates)], 1, CalcVertType, Vertebrates, OhnOG3R))+1]
Counts$VertType[rownames(Counts) %in% OhnOG2R] <- "Ohnologs"
table(Counts$VertType)
Counts$AmphType <- c("Missing", "Single-copy", "Small-scale\nduplicates")[unlist(apply(Counts[,Ampioxus], 1, function(x){r <- max(x); if(r>1){r=2}; return(r)}))+1]
table(Counts$AmphType)
Counts$BlanType <- c("Missing", "Single-copy", "Small-scale\nduplicates")[unlist(lapply(Counts[,"Blan"], function(x){r <- x; if(r>1){r=2}; return(r)}))+1]
table(Counts$BlanType)
Counts$Sum <- apply(Counts[,Species], 1, sum)
Counts$MeanVert <- rowMeans(Counts[,Vertebrates])
Counts.BlanVert <- Counts[which(Counts$BlanType != "Missing" | Counts$VertType != "Missing"),]

# Ortogroup to gene for all species
system_out <- system(paste("cat ", NamesFile," | awk -F '\t' '{if(NR==1){split($0,h,\"\t\");next} split($0,a,\"\t\"); for(i = 2; i <= length(a); ++i){split(a[i],g,\" \"); for(j in g){print a[1]\"\t\"h[i]\"\t\"g[j]}}}'"), intern=T)
OG2Gene <- read.table(text=system_out, h=F, sep = "\t")
colnames(OG2Gene) <- c("OG", "Species", "Gene")
OG2Gene$Species <- Species[match(OG2Gene$Species, paste0(SpeciesLongNames, ".fa"))]

# Gene coordinates
Coord <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("Gene","Chr","Start", "End"))
for(s in c(1:length(Species))){
	system_out <- system(paste0("cat ",GTFFolder, "/", SpeciesLongNames[s], ".gtf | sed 's/[^\t]\\+\\?gene_id \"\\(\\S\\+\\)\";.*/\\1/g' | awk -F '\\t' '{print $9\"\t\"$1\"\t\"$4\"\t\"$5}' | sort -k1,1 -k2,2V -k3,3V | awk -F '\\t' '{if($1==gene){end=$4}else{print gene\"\t\"chr\"\t\"start\"\t\"end\"\t\"pos; if($2==chr){pos=pos+1}else{pos=1}; gene=$1; chr=$2; start=$3; end=$4}}END{print gene\"\t\"chr\"\t\"start\"\t\"end\"\t\"pos;}' | tail -n +2"), intern=T)
	tmp <- read.table(text=system_out, h=F, sep = "\t")
	print(head(tmp))
	Coord <- rbind(Coord, tmp)
}
colnames(Coord) <- c("Gene","Chr","Start", "End", "Position")
GeneInfo <- merge(Coord, OG2Gene, by = "Gene")
GeneInfo <- merge(GeneInfo, Counts, by = "OG")

# Calc parameters per ortogroup per species
for(s in c(1:length(Species))){
	supGeneInfo <- GeneInfo[which(GeneInfo$Species==Species[s]),]
	print(Species[s])
	# Number of chr per OG
	Counts[,paste0("NumChrs",Species[s])] <- unlist(lapply(Counts$OG, function(x){length(unique(supGeneInfo$Chr[which(supGeneInfo$OG==x)]))}))
	# Max distance between pairs of genes in OG
	Counts[,paste0("MaxDist",Species[s])]  <- unlist(lapply(Counts$OG, CalcMaxDist, supGeneInfo))
	# If genes in OG are consecutive
	Counts[,paste0("Consecutive",Species[s])]  <- unlist(lapply(Counts$OG, CalcIfConsecutive, supGeneInfo))
}

NumGenesSpecies <- rep(0, length(Species))
for(sp in c(1:length(Species))){
	system_out <- system(paste0("cat ", ProteomesFolder, "/", SpeciesLongNames[sp], ".fa | grep '>' | sort | uniq | wc -l"), intern=T)
	NumGenesSpecies[sp] <- read.table(text=system_out, h=F, sep = "\t")
}
NumGenesSpecies <- unlist(NumGenesSpecies)

# GO terms
GO <- read.delim(GOlistFile, h=F, sep = '\t')
colnames(GO) <- c("GO", "Class", "Name")
GO <- unique(GO[,c("GO", "Name")])
GO$Hsap <- unlist(lapply(GO$GO, GetNumberOfGOtermGenes, "Hsap", Counts, OG2Gene, GOlistsFolder))
GO$Blan <- unlist(lapply(GO$GO, GetNumberOfGOtermGenes, "Blan", Counts, OG2Gene, GOlistsFolder))
GO <- GO[which(GO$Hsap>=minGOsizeHsap & GO$Blan>=minGOsizeHsap),]
GO <- GO[order(GO$Hsap),]
GO$HsapD <- unlist(lapply(GO$GO, GetNumberOfGOtermGenes, "Hsap", Counts[which(Counts$VertType=="Small-scale\nduplicates"),], OG2Gene, GOlistsFolder))
GO$HsapS <- unlist(lapply(GO$GO, GetNumberOfGOtermGenes, "Hsap", Counts[which(Counts$VertType=="Single-copy"),], OG2Gene, GOlistsFolder))
GO$HsapO <- unlist(lapply(GO$GO, GetNumberOfGOtermGenes, "Hsap", Counts[which(Counts$VertType=="Ohnologs"),], OG2Gene, GOlistsFolder))
GO$BlanD <- unlist(lapply(GO$GO, GetNumberOfGOtermGenes, "Blan", Counts[which(Counts$BlanType=="Small-scale\nduplicates"),], OG2Gene, GOlistsFolder))
GO$BlanS <- unlist(lapply(GO$GO, GetNumberOfGOtermGenes, "Blan", Counts[which(Counts$BlanType=="Single-copy"),], OG2Gene, GOlistsFolder))
print(head(GO))

## Gene expression
# Read Drer & Blan gene information previously processed or process it and write itinto a file
DrerGeneData <- prepare_DrerGeneData(DrerGeneDataFile, DrerBgeeDataFile, MatchingTissues, ResultsFolder)
print(head(DrerGeneData))
BlanGeneData <- prepare_BlanGeneData(BlanGeneDataFile, RNASeqMetadataFile, TPMFile)
print(head(BlanGeneData))

# Read gene pairs and the gene copy number in each species
head <- unlist(strsplit(system(paste0("head -1 ", NamesFile, " | sed 's/#//g'"), intern=T), "\t"))
columnS1 <- which(head==paste0(SpeciesGRefs[which(ShortSpeciesNames==Species1)],".fa"))
columnS2 <- which(head==paste0(SpeciesGRefs[which(ShortSpeciesNames==Species2)],".fa"))
columnDrer <- which(head==paste0(SpeciesGRefs[which(ShortSpeciesNames=="Drer")],".fa"))
columnGgal <- which(head==paste0(SpeciesGRefs[which(ShortSpeciesNames=="Ggal")],".fa"))
columnMmus <- which(head==paste0(SpeciesGRefs[which(ShortSpeciesNames=="Mmus")],".fa"))
columnHsap <- which(head==paste0(SpeciesGRefs[which(ShortSpeciesNames=="Hsap")],".fa"))
system_out <- system(paste0("tail -n +2 ", NamesFile, " | awk -F '\t' '{print $1\"\t\"$",columnS1,"\"\t\"$",columnDrer,"\"\t\"$",columnGgal,"\"\t\"$",columnMmus,"\"\t\"$",columnHsap,"}' | awk -F '\t' '{split($2,a,\" \"); vcn=0; for(v = 3; v <= 6; ++v){split($v,b,\" \"); if(length(b)>vcn){vcn=length(b)}}; for(i in a){print $1\"\t\"a[i]\"\t\"length(a)\"\t\"vcn}}'"), intern=T)
GeneCN <- read.table(text=system_out, h=F, sep = "\t")
colnames(GeneCN) <- c("OG","Gene","CNBlan","CNVert")
GeneCN$TypeBlan <- c("SC", "D")[as.numeric(GeneCN$CNBlan>1)+1]
GeneCN$TypeVert <- c("SC", "D")[as.numeric(GeneCN$CNVert>1)+1]
GeneCN$TypeVert[which(GeneCN$CNVert==0)] <- "M"
GeneCN$TypeVert[which(GeneCN$OG %in% OhnOG2R)] <- "O"

# Construct gene pairs Blan - Drer
system_out <- system(paste0("tail -n +2 ", NamesFile, " | awk -F '\t' '{print $1\"\t\"$", columnS1, "\"\t\"$", columnS2,"}' | awk -F '\t' '{split($2,a,\" \");split($3,b,\" \"); for(i in a){for(j in b){print $1\"\t\"a[i]\"\t\"b[j]\"\t\"length(a)\"\t\"length(b)}}}'"), intern=T)
GenePairs <- read.table(text=system_out, h=F, sep = "\t")
colnames(GenePairs) <- c("OG","Gene1","Gene2","CopyNumber1","CopyNumber2")
GenePairs$Type1 <- c("SC", "D")[as.numeric(GenePairs$CopyNumber1>1)+1]
GenePairs$Type2 <- c("SC", "D")[as.numeric(GenePairs$CopyNumber2>1)+1]
GenePairs$Type2[which(GenePairs$OG %in% OhnOG2R)] <- "O"
print(head(GenePairs))

# Construct OG pxpression profile for species pair
OGProfile <- unique(GenePairs[,c("OG","Type1","Type2")])
for(t in c(1:length(MatchingTissues[,1]))){
	GenePairs[,paste0(MatchingTissues$Name[t],1)] <- BlanGeneData[match(GenePairs$Gene1, BlanGeneData$Gene),MatchingTissues$Blan[t]]>TPM.threshold
	GenePairs[,paste0(MatchingTissues$Name[t],2)] <- DrerGeneData[match(GenePairs$Gene2, DrerGeneData$Gene),paste0(MatchingTissues$Name[t],"Presence")]=="present"
	OGProfile[,paste0(MatchingTissues$Name[t],1)] <- unlist(lapply(OGProfile[,1], function(x){sum(GenePairs[which(GenePairs$OG==x),paste0(MatchingTissues$Name[t],1)])>0}))
	OGProfile[,paste0(MatchingTissues$Name[t],2)] <- unlist(lapply(OGProfile[,1], function(x){sum(GenePairs[which(GenePairs$OG==x),paste0(MatchingTissues$Name[t],2)])>0}))
}
GenePairs$BlanSum <- rowSums(GenePairs[,paste0(MatchingTissues$Name,1)])
GenePairs$DrerSum <- rowSums(GenePairs[,paste0(MatchingTissues$Name,2)])
GenePairs$DiffDom <- GenePairs$BlanSum-GenePairs$DrerSum
OGProfile$BlanSum <- rowSums(OGProfile[,paste0(MatchingTissues$Name,1)])
OGProfile$DrerSum <- rowSums(OGProfile[,paste0(MatchingTissues$Name,2)])
OGProfile$DiffDom <- OGProfile$BlanSum-OGProfile$DrerSum

OGProfile <- merge(OGProfile, Counts[,c("OG","NumChrsBlan","MaxDistBlan","ConsecutiveBlan","NumChrsDrer","MaxDistDrer","ConsecutiveDrer")], by = "OG")
GenePairs <- merge(GenePairs, Counts[,c("OG","NumChrsBlan","MaxDistBlan","ConsecutiveBlan","NumChrsDrer","MaxDistDrer","ConsecutiveDrer")], by = "OG")

GO$BlanD.AExp <- unlist(lapply(GO$GO, CalMeanExpGO, BlanGeneData, "MeanAdult", Counts$OG[which(Counts$BlanType=="Small-scale\nduplicates")], OG2Gene, GOlistsFolder))
GO$BlanD.EExp <- unlist(lapply(GO$GO, CalMeanExpGO, BlanGeneData, "MeanEmbr", Counts$OG[which(Counts$BlanType=="Small-scale\nduplicates")], OG2Gene, GOlistsFolder))
GO$BlanSC.AExp <- unlist(lapply(GO$GO, CalMeanExpGO, BlanGeneData, "MeanAdult", Counts$OG[which(Counts$BlanType=="Single-copy")], OG2Gene, GOlistsFolder))
GO$BlanSC.EExp <- unlist(lapply(GO$GO, CalMeanExpGO, BlanGeneData, "MeanEmbr", Counts$OG[which(Counts$BlanType=="Single-copy")], OG2Gene, GOlistsFolder))

# Pfam
PfamFolder <- "Results/PfamDomainSearch"
system_out <- system(paste0("cat ",PfamFolder, "/",SpeciesLongNames[which(Species=="Blan")],".tbl | grep -v '^#' | awk -F ' ' '{line=$1\"\t\"$2\"\t\"$3\"\t\"$5\"\t\"$6\"\t\"$19; for(i=20;i<=NF;i++){line=line\" \"$i} print line}'"), intern=T)
PfamBlan <- read.table(text=system_out, h=F, sep = "\t", quote = "")
colnames(PfamBlan) <- c("Taget","ID","Gene","Evalue","Score","Description")

system_out <- system(paste0("cat ",PfamFolder, "/",SpeciesLongNames[which(Species=="Mmus")],".tbl | grep -v '^#' | awk -F ' ' '{line=$1\"\t\"$2\"\t\"$3\"\t\"$5\"\t\"$6\"\t\"$19; for(i=20;i<=NF;i++){line=line\" \"$i} print line}'"), intern=T)
PfamHsap <- read.table(text=system_out, h=F, sep = "\t", quote = "")
colnames(PfamHsap) <- c("Taget","ID","Gene","Evalue","Score","Description")

Pfam <- as.data.frame(unique(rbind(PfamHsap[,c("Taget","ID")], PfamBlan[,c("Taget","ID")])))
Pfam$Blan <- unlist(lapply(Pfam$ID, function(x){}))
PfamBlan$Gene[which(PfamBlan$ID=="PF13927.9")]
OG2Gene$OG[which(OG2Gene$Gene %in% PfamBlan$Gene[which(PfamBlan$ID=="PF13927.9")])]












###########################################################################
###########################################################################
#### Printing statistics

# Amphioxus - vertebrate shared OG numbers
paste("Shared V-A", sum(Counts$VertType!="Missing" & Counts$AmphType!="Missing"))
paste("V", sum(Counts$VertType!="Missing"))
paste("A", sum(Counts$AmphType!="Missing"))
paste("Shared/V*100", sum(Counts$VertType!="Missing" & Counts$AmphType!="Missing")/sum(Counts$VertType!="Missing")*100)
paste("Shared/A*100", sum(Counts$VertType!="Missing" & Counts$AmphType!="Missing")/sum(Counts$AmphType!="Missing")*100)

# Amphioxus specfic OG numbers
Counts.A <- Counts[which(Counts$VertType=="Missing" & Counts$AmphType!="Missing"),]
paste("A specific", length(Counts.A[,1]))
paste("Shared all A", sum(unlist(apply(Counts.A[,Ampioxus], 1, min))>0))
paste("Shared all A/A specific*100", sum(unlist(apply(Counts.A[,Ampioxus], 1, min))>0)/length(Counts.A[,1])*100)
paste("blan-bbel specific", sum(Counts.A$Blan>0 & Counts.A$Bflo==0 & Counts.A$Bbel>0))
paste("blan-bbel specific/A specific*100", sum(Counts.A$Blan>0 & Counts.A$Bflo==0 & Counts.A$Bbel>0)/length(Counts.A[,1])*100)
paste("blan-bflo specific", sum(Counts.A$Blan>0 & Counts.A$Bflo>0 & Counts.A$Bbel==0))
paste("blan-bflo specific/A specific*100", sum(Counts.A$Blan>0 & Counts.A$Bflo>0 & Counts.A$Bbel==0)/length(Counts.A[,1])*100)

for(sp in c(1:length(Species))){
	OG2Gene.sp <- OG2Gene[which(OG2Gene$Species==Species[sp]),]
	toprint <- c(1:7)
	toprint[1] <- Species[sp]
	toprint[2] <- sum(Counts[,Species[sp]]>1)
	toprint[3] <- sum(Counts[,Species[sp]]>1)/sum(Counts[,Species[sp]]>0)*100
	toprint[4] <- length(unique(OG2Gene.sp$Gene[which(OG2Gene.sp$OG %in% rownames(Counts[which(Counts[,Species[sp]]>1),]))]))
	toprint[5] <- length(unique(OG2Gene.sp$Gene[which(OG2Gene.sp$OG %in% rownames(Counts[which(Counts[,Species[sp]]>1),]))]))/NumGenesSpecies[sp]*100
	toprint[6] <- NumGenesSpecies[sp]-length(unique(OG2Gene.sp$Gene))
	toprint[7] <- (NumGenesSpecies[sp]-length(unique(OG2Gene.sp$Gene)))/NumGenesSpecies[sp]*100
	print(paste(toprint, collapse=" "))
}

ContingencyTableCN(
	Counts$Blan[which(Counts$BlanType=="Small-scale\nduplicates" & Counts$VertType=="Small-scale\nduplicates")], 
	Counts$MeanVert[which(Counts$BlanType=="Small-scale\nduplicates" & Counts$VertType=="Small-scale\nduplicates")],
	paste0(ResultsFolder, "/ContingencyTableCN_SS.txt"))

ContingencyTableCN(
	Counts$Blan[which(Counts$BlanType=="Small-scale\nduplicates" & Counts$VertType=="Ohnologs")], 
	Counts$MeanVert[which(Counts$BlanType=="Small-scale\nduplicates" & Counts$VertType=="Ohnologs")],
	paste0(ResultsFolder, "/ContingencyTableCN_O.txt"))

ContingencyTableCN(
	Counts$Blan[which(Counts$BlanType=="Small-scale\nduplicates" & (Counts$VertType=="Small-scale\nduplicates" | Counts$VertType=="Ohnologs"))], 
	Counts$MeanVert[which(Counts$BlanType=="Small-scale\nduplicates" & (Counts$VertType=="Small-scale\nduplicates" | Counts$VertType=="Ohnologs"))],
	paste0(ResultsFolder, "/ContingencyTableCN_SSO.txt"))

GO$PercentHsapD = GO$HsapD/GO$Hsap*100
GO$PercentHsapO = GO$HsapO/GO$Hsap*100
GO$PercentBlanD = GO$BlanD/GO$Blan*100
write.table(GO[order(GO$PercentBlanD),], file = paste(ResultsFolder, "/GOterms_PropData.txt", sep =""), quote = F, sep="\t", col.names = TRUE, row.names = FALSE)




###########################################################################
###########################################################################
### Plotting
pdf(paste(ResultsFolder, "/DuplicatesOntologyExpression.pdf", sep=""), width=15, height=10)
par(mar=c(10,10,5,5),oma=c(1,1,1,1), yaxs='i', xaxs='i')

# Blan vs. vertebrate gene type 
layout(matrix(c(1,2,3,4,5,6),nrow=2,ncol=3,byrow=T), widths=c(1), heights=c(1), TRUE)
BarPlotVertTypesInBlanGenes(Counts.BlanVert[which(Counts.BlanVert$VertType!="Missing" & Counts.BlanVert$BlanType!="Missing"),], BlanTypes[which(BlanTypes!="Missing")], VertTypes[which(VertTypes!="Missing")], VertTypes.col[which(VertTypes!="Missing")])
PlotHypergeomTest_VertBlanTypes(Counts.BlanVert[which(Counts.BlanVert$VertType!="Missing" & Counts.BlanVert$BlanType!="Missing"),], VertTypes[which(VertTypes!="Missing")], BlanTypes[which(BlanTypes!="Missing")], VertTypes.col[which(VertTypes!="Missing")], PQvalThreshold, ResultsFolder)
plot.new()
legend("bottomright", c("In vertebrates", VertTypes[which(VertTypes!="Missing")]), pch=c(NA,15,15,15,15), text.col="black", col=c(NA, VertTypes.col[which(VertTypes!="Missing")]), bty = "n", pt.cex=2, cex=1.5, xjust = 0, yjust = 0)#
BarPlotVertTypesInBlanGenes(Counts.BlanVert, BlanTypes, VertTypes, VertTypes.col)
PlotHypergeomTest_VertBlanTypes(Counts.BlanVert, VertTypes, BlanTypes, VertTypes.col, PQvalThreshold, ResultsFolder)
plot.new()
legend("bottomright", c("In vertebrates", "Missing", "Single-copy", "Ohnologs", "Small-scale\nduplicates"), pch=c(NA,15,15,15,15), text.col="black", col=c(NA, VertTypes.col), bty = "n", pt.cex=2, cex=1.5, xjust = 0, yjust = 0)

# GO term in duplicates Hsap vs. Blan 
layout(matrix(c(1,2,3,4,5,6),nrow=2,ncol=3,byrow=T), widths=c(1), heights=c(1), TRUE)
ScatterPlotContour(GO$HsapD/GO$Hsap*100, GO$BlanD/GO$Blan*100, GO$Hsap, "black", "% of small-scale duplicates\nH. sapiens", "B. lanceolatum\n% of small-scale duplicates", c(0,100), c(0,100))
ScatterPlotContour(GO$HsapO/GO$Hsap*100, GO$BlanD/GO$Blan*100, GO$Hsap, "black", "% of ohnologs\nH. sapiens", "B. lanceolatum\n% of small-scale duplicates", c(0,100), c(0,100))
ScatterPlotContour((GO$HsapO+GO$HsapD)/GO$Hsap*100, GO$BlanD/GO$Blan*100, GO$Hsap, "black", "% of ohnologs\nH. sapiens", "B. lanceolatum\n% of small-scale duplicates", c(0,100), c(0,100))
ScatterPlotContour(GO$HsapS/GO$Hsap*100, GO$BlanS/GO$Blan*100, GO$Hsap, "black", "% of ohnologs\nH. sapiens", "B. lanceolatum\n% of small-scale duplicates", c(0,100), c(0,100))

layout(matrix(c(1,2),nrow=2,ncol=1,byrow=T), widths=c(15), heights=c(5), TRUE)
TandemIntraInterPerSpecies(Counts, Species, SpType, MaxTandemDist, "MaxDistance")
layout(matrix(c(1,2),nrow=2,ncol=1,byrow=T), widths=c(15), heights=c(5), TRUE)
TandemIntraInterPerSpecies(Counts, Species, SpType, MaxTandemDist, "Consecutive")

layout(matrix(c(1,2,5,3,4,6),nrow=2,ncol=3,byrow=T), widths=c(1.5,1.5), heights=c(1.5, 1.5), TRUE)
BoxPlot_BlanTypes_VertTypes(BlanGeneData, GeneCN, "MeanAdult", "Mean adult expression", c(0,100), VertTypes, VertTypes.abr, VertTypes.col)
BoxPlot_BlanTypes_VertTypes(BlanGeneData, GeneCN, "MeanEmbr", "Mean embrionic expression", c(0,100), VertTypes, VertTypes.abr, VertTypes.col)
BoxPlot_BlanTypes_VertTypes(BlanGeneData, GeneCN, "TauTissues", "Tau among tissues", c(0,1), VertTypes, VertTypes.abr, VertTypes.col)
BoxPlot_BlanTypes_VertTypes(BlanGeneData, GeneCN, "TauEmbAge", "Tau among developmental stages", c(0,1), VertTypes, VertTypes.abr, VertTypes.col)

layout(matrix(c(1,2,3,4,5,6),nrow=2,ncol=3,byrow=T), widths=c(1), heights=c(1), TRUE)
Hist_ExpressDomanis(GenePairs, OGProfile, "D", "SC", MatchingTissues, "B. lanceolatum donains - D. rerio domains", "Relative number of\npairwise comparisons", "B. lanceolatum specific\ngene duplicates")
Hist_ExpressDomanis(GenePairs, OGProfile, "SC", "D", MatchingTissues, "B. lanceolatum donains - D. rerio domains", "Relative number of\npairwise comparisons", "D. rerio specific\nsmall scale gene duplicates")
Hist_ExpressDomanis(GenePairs, OGProfile, "SC", "O", MatchingTissues, "B. lanceolatum donains - D. rerio domains", "Relative number of\npairwise comparisons", "D. rerio specific\nohnolog gene duplicates")

layout(matrix(c(1,2,3,4,5,6),nrow=2,ncol=3,byrow=T), widths=c(1), heights=c(1), TRUE)
Hist_ExpressDomanis(
	GenePairs[which(GenePairs$Type1=="SC" | (GenePairs$Type1=="D" & GenePairs$NumChrsBlan>1)),], 
	OGProfile[which(OGProfile$Type1=="SC" | (OGProfile$Type1=="D" & OGProfile$NumChrsBlan>1)),], "D", "SC", MatchingTissues, "B. lanceolatum donains - D. rerio domains", "Relative number of\npairwise comparisons", "B. lanceolatum specific\nmultichormosomal\ngene duplicates")
Hist_ExpressDomanis(
	GenePairs[which(GenePairs$Type1=="SC" | (GenePairs$Type1=="D" & GenePairs$NumChrsBlan==1 & GenePairs$ConsecutiveBlan==FALSE)),], 
	OGProfile[which(OGProfile$Type1=="SC" | (OGProfile$Type1=="D" & OGProfile$NumChrsBlan==1 & OGProfile$ConsecutiveBlan==FALSE)),], "D", "SC", MatchingTissues, "B. lanceolatum donains - D. rerio domains", "Relative number of\npairwise comparisons", "B. lanceolatum specific\ndistant monochromosomal\ngene duplicates")
Hist_ExpressDomanis(
	GenePairs[which(GenePairs$Type1=="SC" | (GenePairs$Type1=="D" & GenePairs$NumChrsBlan==1 & GenePairs$ConsecutiveBlan==TRUE)),], 
	OGProfile[which(OGProfile$Type1=="SC" | (OGProfile$Type1=="D" & OGProfile$NumChrsBlan==1 & OGProfile$ConsecutiveBlan==TRUE)),], "D", "SC", MatchingTissues, "B. lanceolatum donains - D. rerio domains", "Relative number of\npairwise comparisons", "B. lanceolatum specific\ntandem\ngene duplicates")

layout(matrix(c(1,2,3,4,5,6),nrow=2,ncol=3,byrow=T), widths=c(1), heights=c(1), TRUE)
DifferenceWithSC(GenePairs, OGProfile, MatchingTissues, "Branch specific gene duplicates", "Difference with\nsingle-copy genes distribution")
DifferenceWithSC_TandemIntraInter(GenePairs, OGProfile, MatchingTissues, "Branch specific gene duplicates", "Difference with\nsingle-copy genes distribution")

layout(matrix(c(1,2,3,4,5,6),nrow=2,ncol=3,byrow=T), widths=c(1), heights=c(1), TRUE)
ScatterPlot_GOExp(GO$BlanD.AExp, GO$BlanSC.AExp, GO$Blan, "black", "Mean adult expression\nof duplicate genes", "Mean adult expression\nof single-copy genes", c(0,120))
ScatterPlot_GOExp(GO$BlanD.EExp, GO$BlanSC.EExp, GO$Blan, "black", "Mean embrionic expression\nof duplicate genes", "Mean embrionic expression\nof single-copy genes", c(0,120))


dev.off()


###########################################################################
###########################################################################

