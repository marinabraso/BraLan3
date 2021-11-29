#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")

######################################################################
# Libraries & functions

#script <- sub(".*=", "", commandArgs()[4])
#source(paste(substr(script,1, nchar(script)-2), "_functions.r", sep=""))
source("ScriptsPlots/DuplicatesOntologyExpression/DuplicatesOntologyExpression_functions.r")
# source("ScriptsPlots/DuplicatesOntologyExpression/DuplicatesOntologyExpression.r")
library(ghibli)

######################################################################
# Files & folders

ResultsFolder <- "Plots/DuplicatesOntologyExpression"
CountsFile <- "Results/FindOrthologs/AmphVerteb_broccoli/dir_step3/table_OGs_protein_counts.txt"
NamesFile <- "Results/FindOrthologs/AmphVerteb_broccoli/dir_step3/table_OGs_protein_names.txt"
Ohnolog2ROG <- "Results/OhnologListing/2R_Strict_OG_pairs.txt"
Ohnolog3ROG <- "Results/OhnologListing/3R_Strict_OG_pairs.txt"

GOlistsFolder <- "Results/GeneOntology"
GOlistFile <- "Results/GeneOntology/GOlist_MF.list"

GTFFolder <- "Results/FilteringGeneSets/FilteredGTFs"
ProteomesFolder <- "Results/FilteringGeneSets/FilteredProteomes"

BlanGeneDataFile <- paste0(ResultsFolder, "/BlanGeneData.txt")
DrerGeneDataFile <- paste0(ResultsFolder, "/DrerGeneData.txt")
DrerBgeeDataFile <- paste0(ResultsFolder, "/DrerBgeeData.txt")
RNASeqMetadataFile <- "Metadata/Marletaz2018_RNAseq_SRA.txt"
TPMFile <- "Results/GeneExpression/Gene_TPM_kallisto_tximport.tab"

######################################################################
# General parameters
PQvalThreshold <- 0.01 # for hypergeometric test
minGOsize <- 50
MaxTandemDist <- 10000
TPM.threshold <- 1

# Species
Species <- c("Blan", "Bflo", "Bbel", "Drer", "Ggal", "Mmus", "Hsap")
Ampioxus <- c("Blan", "Bflo", "Bbel")
Vertebrates <- c("Drer", "Ggal", "Mmus", "Hsap")
SpeciesLongNames <- c("Branchiostoma_lanceolatum.BraLan3", "Branchiostoma_floridae.Bfl_VNyyK", "Branchiostoma_belcheri.Haploidv18h27", "Danio_rerio.GRCz11", "Gallus_gallus.GRCg6a", "Mus_musculus.GRCm39", "Homo_sapiens.GRCh38")

VertTypes <- c("Single-copy", "Ohnologs", "Small-scale\nduplicates", "Missing")
VertTypes.abr <- c("SC", "O", "D", "M")
VertTypes.col <- c(ghibli_palettes$MarnieMedium2[c(2,4,7)], "royalblue2")
BlanTypes <- c("Single-copy", "Small-scale\nduplicates", "Missing")
BlanTypes.pch <- c(14, 16, 18)
BlanTypes.lty <- c(.5, 1, 2)

# Gene expression analysis
Species1 <- "Blan"
Species2 <- "Drer"

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
OGInfo <- read.table(CountsFile, h=F, sep = "\t", row.names=1)
system_out <- system(paste("head -1 ", CountsFile, " | cut -f2,3,4,5,6,7,8,9,10,11,12 | sed 's/\\([A-Z]\\)[a-z]\\+_\\([a-z][a-z][a-z]\\)[A-Za-z0-9\\._]\\+/\\1\\2/g'"), intern=T)
colnames(OGInfo) <- read.table(text=system_out, h=F, sep = "\t")
OGInfo$OG <- rownames(OGInfo)
OGInfo$VertType <- c("Missing", "Single-copy", "Small-scale\nduplicates")[unlist(apply(OGInfo[,c("OG",Vertebrates)], 1, CalcVertType, Vertebrates, OhnOG3R))+1]
OGInfo$VertType[rownames(OGInfo) %in% OhnOG2R] <- "Ohnologs"
table(OGInfo$VertType)
OGInfo$AmphType <- c("Missing", "Single-copy", "Small-scale\nduplicates")[unlist(apply(OGInfo[,Ampioxus], 1, function(x){r <- max(x); if(r>1){r=2}; return(r)}))+1]
table(OGInfo$AmphType)
OGInfo$BlanType <- c("Missing", "Single-copy", "Small-scale\nduplicates")[unlist(lapply(OGInfo[,"Blan"], function(x){r <- x; if(r>1){r=2}; return(r)}))+1]
table(OGInfo$BlanType)
OGInfo$Sum <- apply(OGInfo[,Species], 1, sum)
OGInfo$MeanVert <- rowMeans(OGInfo[,Vertebrates])

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

# Unify gene coordinates and OGinfo into GeneInfo
GeneInfo <- merge(Coord, OG2Gene, by = "Gene")
GeneInfo <- merge(GeneInfo, OGInfo, by = "OG")

# Calc parameters per ortogroup per species
for(s in c(1:length(Species))){
	supGeneInfo <- GeneInfo[which(GeneInfo$Species==Species[s]),]
	# Number of chr per OG
	OGInfo[,paste0("NumChrs",Species[s])] <- unlist(lapply(OGInfo$OG, function(x){length(unique(supGeneInfo$Chr[which(supGeneInfo$OG==x)]))}))
	# Max distance between pairs of genes in OG
	OGInfo[,paste0("MaxDist",Species[s])]  <- unlist(lapply(OGInfo$OG, CalcMaxDist, supGeneInfo))
	# If genes in OG are consecutive
	OGInfo[,paste0("Consecutive",Species[s])]  <- unlist(lapply(OGInfo$OG, CalcIfConsecutive, supGeneInfo))
}
OGInfo$BlanLType <- NA
OGInfo$BlanLType[which(OGInfo$BlanType=="Small-scale\nduplicates" & OGInfo$NumChrsBlan>1)] <- "Inter"
OGInfo$BlanLType[which(OGInfo$BlanType=="Small-scale\nduplicates" & OGInfo$NumChrsBlan==1 & OGInfo$ConsecutiveBlan==FALSE)] <- "Intra"
OGInfo$BlanLType[which(OGInfo$BlanType=="Small-scale\nduplicates" & OGInfo$NumChrsBlan==1 & OGInfo$ConsecutiveBlan==TRUE)] <- "Tandem"


# Calc num genes per species
NumGenesSpecies <- rep(0, length(Species))
for(sp in c(1:length(Species))){
	system_out <- system(paste0("cat ", ProteomesFolder, "/", SpeciesLongNames[sp], ".fa | grep '>' | sort | uniq | wc -l"), intern=T)
	NumGenesSpecies[sp] <- read.table(text=system_out, h=F, sep = "\t")
}
NumGenesSpecies <- unlist(NumGenesSpecies)

## Gene expression
# Read Drer & Blan gene information previously processed or process it and write itinto a file
DrerGeneData <- prepare_DrerGeneData(DrerGeneDataFile, DrerBgeeDataFile, MatchingTissues, ResultsFolder)
print(head(DrerGeneData))
BlanGeneData <- prepare_BlanGeneData(BlanGeneDataFile, RNASeqMetadataFile, TPMFile)
print(head(BlanGeneData))

# Construct gene pairs Blan - Drer
head <- unlist(strsplit(system(paste0("head -1 ", NamesFile, " | sed 's/#//g'"), intern=T), "\t"))
columnS1 <- which(head==paste0(SpeciesLongNames[which(Species==Species1)],".fa"))
columnS2 <- which(head==paste0(SpeciesLongNames[which(Species==Species2)],".fa"))
system_out <- system(paste0("tail -n +2 ", NamesFile, " | awk -F '\t' '{print $1\"\t\"$", columnS1, "\"\t\"$", columnS2,"}' | awk -F '\t' '{split($2,a,\" \");split($3,b,\" \"); for(i in a){for(j in b){print $1\"\t\"a[i]\"\t\"b[j]\"\t\"length(a)\"\t\"length(b)}}}'"), intern=T)
GenePairs <- read.table(text=system_out, h=F, sep = "\t")
colnames(GenePairs) <- c("OG","Gene1","Gene2","CopyNumber1","CopyNumber2")
GenePairs$BlanType <- OGInfo$BlanType[match(GenePairs$OG, OGInfo$OG)]
GenePairs$VertType <- OGInfo$VertType[match(GenePairs$OG, OGInfo$OG)]
print(head(GenePairs))

# Construct OG pxpression profile for species pair
for(t in c(1:length(MatchingTissues[,1]))){
	GenePairs[,paste0(MatchingTissues$Name[t],1)] <- BlanGeneData[match(GenePairs$Gene1, BlanGeneData$Gene),MatchingTissues$Blan[t]]>TPM.threshold
	GenePairs[,paste0(MatchingTissues$Name[t],2)] <- DrerGeneData[match(GenePairs$Gene2, DrerGeneData$Gene),paste0(MatchingTissues$Name[t],"Presence")]=="present"
	OGInfo[,paste0(MatchingTissues$Name[t],1)] <- unlist(lapply(OGInfo$OG, function(x){sum(GenePairs[which(GenePairs$OG==x),paste0(MatchingTissues$Name[t],1)])>0}))
	OGInfo[,paste0(MatchingTissues$Name[t],2)] <- unlist(lapply(OGInfo$OG, function(x){sum(GenePairs[which(GenePairs$OG==x),paste0(MatchingTissues$Name[t],2)])>0}))
}
GenePairs$BlanSum <- rowSums(GenePairs[,paste0(MatchingTissues$Name,1)])
GenePairs$DrerSum <- rowSums(GenePairs[,paste0(MatchingTissues$Name,2)])
GenePairs$DiffDom <- GenePairs$BlanSum-GenePairs$DrerSum
OGInfo$BlanSum <- rowSums(OGInfo[,paste0(MatchingTissues$Name,1)])
OGInfo$DrerSum <- rowSums(OGInfo[,paste0(MatchingTissues$Name,2)])
OGInfo$DiffDom <- OGInfo$BlanSum-OGInfo$DrerSum

GenePairs <- merge(GenePairs, OGInfo[,c("OG","BlanLType")], by = "OG")


# GO terms
GOInfo <- read.delim(GOlistFile, h=F, sep = '\t')
colnames(GOInfo) <- c("GO", "Class", "Name")
GOInfo <- unique(GOInfo[,c("GO", "Name")])
GOInfo$Hsap <- unlist(lapply(GOInfo$GO, GetNumberOfGOtermGenes, "Hsap", OGInfo, OG2Gene, GOlistsFolder))
GOInfo$Blan <- unlist(lapply(GOInfo$GO, GetNumberOfGOtermGenes, "Blan", OGInfo, OG2Gene, GOlistsFolder))
GOInfo <- GOInfo[which(GOInfo$Hsap>=minGOsize & GOInfo$Blan>=minGOsize),]
GOInfo <- GOInfo[order(GOInfo$Hsap),]
GOInfo$HsapD <- unlist(lapply(GOInfo$GO, GetNumberOfGOtermGenes, "Hsap", OGInfo[which(OGInfo$VertType=="Small-scale\nduplicates"),], OG2Gene, GOlistsFolder))
GOInfo$HsapS <- unlist(lapply(GOInfo$GO, GetNumberOfGOtermGenes, "Hsap", OGInfo[which(OGInfo$VertType=="Single-copy"),], OG2Gene, GOlistsFolder))
GOInfo$HsapO <- unlist(lapply(GOInfo$GO, GetNumberOfGOtermGenes, "Hsap", OGInfo[which(OGInfo$VertType=="Ohnologs"),], OG2Gene, GOlistsFolder))
GOInfo$BlanD <- unlist(lapply(GOInfo$GO, GetNumberOfGOtermGenes, "Blan", OGInfo[which(OGInfo$BlanType=="Small-scale\nduplicates"),], OG2Gene, GOlistsFolder))
GOInfo$BlanS <- unlist(lapply(GOInfo$GO, GetNumberOfGOtermGenes, "Blan", OGInfo[which(OGInfo$BlanType=="Single-copy"),], OG2Gene, GOlistsFolder))
GOInfo$BlanD.AExp <- unlist(lapply(GOInfo$GO, CalMeanExpGO, BlanGeneData, "MeanAdult", OGInfo$OG[which(OGInfo$BlanType=="Small-scale\nduplicates")], OG2Gene, GOlistsFolder))
GOInfo$BlanD.EExp <- unlist(lapply(GOInfo$GO, CalMeanExpGO, BlanGeneData, "MeanEmbr", OGInfo$OG[which(OGInfo$BlanType=="Small-scale\nduplicates")], OG2Gene, GOlistsFolder))
GOInfo$BlanSC.AExp <- unlist(lapply(GOInfo$GO, CalMeanExpGO, BlanGeneData, "MeanAdult", OGInfo$OG[which(OGInfo$BlanType=="Single-copy")], OG2Gene, GOlistsFolder))
GOInfo$BlanSC.EExp <- unlist(lapply(GOInfo$GO, CalMeanExpGO, BlanGeneData, "MeanEmbr", OGInfo$OG[which(OGInfo$BlanType=="Single-copy")], OG2Gene, GOlistsFolder))
GOInfo$PercentHsapD <- GOInfo$HsapD/GOInfo$Hsap*100
GOInfo$PercentHsapO <- GOInfo$HsapO/GOInfo$Hsap*100
GOInfo$PercentBlanD <- GOInfo$BlanD/GOInfo$Blan*100


print(head(GOInfo))

# Pfam
PfamFolder <- "Results/PfamDomainSearch"
system_out <- system(paste0("cat ",PfamFolder, "/",SpeciesLongNames[which(Species=="Blan")],".tbl | grep -v '^#' | awk -F ' ' '{line=$1\"\t\"$2\"\t\"$3\"\t\"$5\"\t\"$6\"\t\"$19; for(i=20;i<=NF;i++){line=line\" \"$i} print line}'"), intern=T)
PfamBlan <- read.table(text=system_out, h=F, sep = "\t", quote = "")
colnames(PfamBlan) <- c("Taget","ID","Gene","Evalue","Score","Description")

system_out <- system(paste0("cat ",PfamFolder, "/",SpeciesLongNames[which(Species=="Mmus")],".tbl | grep -v '^#' | awk -F ' ' '{line=$1\"\t\"$2\"\t\"$3\"\t\"$5\"\t\"$6\"\t\"$19; for(i=20;i<=NF;i++){line=line\" \"$i} print line}'"), intern=T)
PfamHsap <- read.table(text=system_out, h=F, sep = "\t", quote = "")
colnames(PfamHsap) <- c("Taget","ID","Gene","Evalue","Score","Description")

Pfam <- as.data.frame(unique(rbind(PfamHsap[,c("Taget","ID")], PfamBlan[,c("Taget","ID")])))
Pfam$Blan <- unlist(lapply(Pfam$ID, function(x){length(unique(PfamBlan$Gene[which(PfamBlan$ID==x)]))}))
Pfam$BlanSC <- unlist(lapply(Pfam$ID, function(x){length(GeneInfo$BlanType[which(GeneInfo$Gene %in% unique(PfamBlan$Gene[which(PfamBlan$ID==x)]))]=="Single-copy")}))
Pfam$BlanD <- unlist(lapply(Pfam$ID, function(x){length(GeneInfo$BlanType[which(GeneInfo$Gene %in% unique(PfamBlan$Gene[which(PfamBlan$ID==x)]))]=="Small-scale\nduplicates")}))
Pfam$Hsap <- unlist(lapply(Pfam$ID, function(x){length(unique(PfamBlan$Gene[which(PfamBlan$ID==x)]))}))
Pfam$HsapSC <- unlist(lapply(Pfam$ID, function(x){length(GeneInfo$BlanType[which(GeneInfo$Gene %in% unique(PfamBlan$Gene[which(PfamBlan$ID==x)]))]=="Single-copy")}))
Pfam$HsapD <- unlist(lapply(Pfam$ID, function(x){length(GeneInfo$BlanType[which(GeneInfo$Gene %in% unique(PfamBlan$Gene[which(PfamBlan$ID==x)]))]=="Small-scale\nduplicates")}))


PfamBlan$Gene[which(PfamBlan$ID=="PF13927.9")]
GeneInfo$BlanType[which(GeneInfo$Gene %in% PfamBlan$Gene[which(PfamBlan$ID=="PF13927.9")])]



######################################################################
######################################################################
# Write files
cat("Write files\n")

write.table(OGInfo, file = paste(ResultsFolder, "/OGInfo", sep =""), quote = F, sep="\t", col.names = TRUE, row.names = TRUE)
write.table(GeneInfo, file = paste(ResultsFolder, "/GeneInfo.txt", sep =""), quote = F, sep="\t", col.names = TRUE, row.names = TRUE)
write.table(NumGenesSpecies, file = paste(ResultsFolder, "/NumGenesSpecies.txt", sep =""), quote = F, sep="\t", col.names = TRUE, row.names = TRUE)
write.table(GenePairs, file = paste(ResultsFolder, "/GenePairs.txt", sep =""), quote = F, sep="\t", col.names = TRUE, row.names = TRUE)
write.table(GOInfo, file = paste(ResultsFolder, "/GOInfo.txt", sep =""), quote = F, sep="\t", col.names = TRUE, row.names = TRUE)
write.table(GOInfo[order(GOInfo$PercentBlanD),], file = paste(ResultsFolder, "/GOterms_PropData.txt", sep =""), quote = F, sep="\t", col.names = TRUE, row.names = FALSE)










