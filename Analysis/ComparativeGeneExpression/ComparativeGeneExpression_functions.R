

add_hist_logYaxis <- function(h, col){
	#for(i in c(1:length(h$counts))){
	#	polygon(c(h$breaks[i],h$breaks[i+1],h$breaks[i+1],h$breaks[i]),c(0,0,log(h$counts[i]),log(h$counts[i])), col=modif_alpha(col,.2), border=col)
	#}
	lines(h$mids, log(h$counts), col=col, lwd=4)
	points(h$mids, log(h$counts), col=col, pch=16, cex=3)
}

add_relative_freq_line <- function(vec, b, col){
	h <- hist(vec, breaks=b, plot=FALSE)
	lines(h$mids, h$counts/length(vec), col=col, lwd=4)
}

add_joined_relative_hist_line <- function(b, cols, main, xlab, vec1, vec2, vec3=NULL){
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,1), xlim=c(min(b), max(b)), col=NA)
	mtext(main, side = 3, line = 1, cex=1)
	mtext("Relative number of pairwise comparisons", side = 2, line = 3, cex=1)
	mtext(xlab, side = 1, line = 3, cex=1)
	mids <- add_relative_hist(vec1, b, cols[1])
	add_relative_hist(vec2, b, cols[2])
	add_relative_hist(vec3, b, cols[3], dens=20)
	axis(1, at = mids, lwd.ticks=1, las=1, cex.axis=1)
	axis(2, at = seq(0,1,.2), lwd.ticks=1, las=1, cex.axis=1)
	box()
}

add_relative_hist <- function(vec, b, col, dens=NULL){
	h <- hist(vec, breaks=b, plot=FALSE)
	h$counts <- h$counts/length(vec)
	print(h$counts)
	for(i in c(1:length(h$counts))){
		polygon(c(h$breaks[i],h$breaks[i+1],h$breaks[i+1],h$breaks[i]),c(0,0,h$counts[i],h$counts[i]), col=modif_alpha(col,.3), border=col, density=dens, lwd=2)
	}
	return(h$mids)
}

get_spearman_corr <- function(pairs, tissue){
	vec1 <- pairs[,paste0(tissue,1)]
	vec2 <- pairs[,paste0(tissue,2)]
	tmp <- cbind(vec1[which(!is.na(vec1) & !is.na(vec2))], vec2[which(!is.na(vec1) & !is.na(vec2))])
	return(cor(tmp[,1], tmp[,2], method="spearman"))
}

get_spearman_corr_presence_absence <- function(pairs, tissue, thresh){
	vec1 <- pairs[,paste0(tissue,1)]>=thresh
	vec2 <- pairs[,paste0(tissue,2)]>=thresh
	tmp <- cbind(vec1[which(!is.na(vec1) & !is.na(vec2))], vec2[which(!is.na(vec1) & !is.na(vec2))])
	return(cor(tmp[,1], tmp[,2], method="spearman"))
}

modif_alpha <- function(col, alpha=.5){
  if(missing(col))
  stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha=alpha))  
}

sum_union_of_patterns <- function(OG, spnum, pairs, tissues, tpmthresh){
	return(sum(colSums(pairs[which(pairs$OG == OG),paste0(tissues$Name,spnum)]>tpmthresh)>0))
}

absolute_difference_union_of_patterns <- function(OG, pairs, tissues, tpmthresh){
	pattern1 <- colSums(pairs[which(pairs$OG == OG),paste0(tissues$Name,1)]>tpmthresh)>0
	pattern2 <- colSums(pairs[which(pairs$OG == OG),paste0(tissues$Name,2)]>tpmthresh)>0
	return(sum(abs(pattern1-pattern2)))
}

difference_union_of_patterns_per_tissue_gp <- function(x, pairs, tissues, tpmthresh){
	pattern1 <- pairs[x,paste0(tissues$Name, 1)]>tpmthresh
	pattern2 <- pairs[x,paste0(tissues$Name, 2)]>tpmthresh
	names(pattern1) <- tissues$Name
	DifTissues <- names(pattern1)[which(colSums(rbind(pattern1, pattern2))==1)]
	if(length(DifTissues)==0){DifTissues=NA}
	return(list(DifTissues))
}

difference_union_of_patterns_per_tissue_c <- function(OG, pairs, tissues, tpmthresh){
	pattern1 <- colSums(pairs[which(pairs$OG == OG),paste0(tissues$Name,1)]>tpmthresh)>0
	pattern2 <- colSums(pairs[which(pairs$OG == OG),paste0(tissues$Name,2)]>tpmthresh)>0
	names(pattern1) <- tissues$Name
	DifTissues <- names(pattern1)[which(colSums(rbind(pattern1, pattern2))==1)]
	if(length(DifTissues)==0){DifTissues=NA}
	return(list(DifTissues))
}

sum_TPM_per_species_tissue <- function(OG, spnum, pairs, tissue){
	if(length(pairs[which(pairs$OG == OG),1])==0){
		return(0) # Consider TPM=0 if no gene expression information 
	}else if(length(pairs[which(pairs$OG == OG),1])>1){
		return(sum(pairs[which(pairs$OG == OG),paste0(tissue,spnum)]))
	}else{
		return(pairs[which(pairs$OG == OG),paste0(tissue,spnum)])
	}
}

prepare_DrerGeneData <- function(file, bgeefile, sname, tissues, rfolder, pfolder){
	if(file.exists(file)){
		Data <- read.delim(file, h=T, stringsAsFactors=F)
	}else{
		library(BgeeDB)
		setwd(paste0("./", rfolder))
		DrerBgee <- Bgee$new(species=sname, dataType="rna_seq")
		DrerBgeeData <- getData(DrerBgee)
		write.table(DrerBgeeData, file = bgeefile, quote = F, sep="\t", col.names = TRUE, row.names = TRUE)
		Data <- as.data.frame(cbind(unique(DrerBgeeData$Gene.ID)))
		colnames(Data) <- c("Gene")
		for(i in c(1:length(tissues[,1]))){
			print(tissues$Name[i])
			system_out <- system(paste0("cat ", bgeefile, " | sed 's/\"//g' | awk -F'\\t' -v tissue=\"", tissues$Drer[i], "\" '{if($7==tissue){print $5\"\\t\"$13\"\\t\"$16}}' | sort -k1,1 | awk '{if(g!=$1){if(NR!=1){mTPM=mTPM/num;} print g\"\\t\"mTPM\"\\t\"pres; g=$1; num=1;mTPM=$2;pres=$3}else{num=num+1;mTPM=mTPM+$2;if($3==\"present\"){pres=$3}}}END{mTPM=mTPM/num; print g\"\\t\"mTPM\"\\t\"pres;}' | tail -n +2"), intern=T)
			tData <- read.table(text=system_out, h=F, sep = "\t")
			Data[,paste0(tissues$Name[i],"TPM")] <- tData[match(Data$Gene, tData[,1]),2]
			Data[,paste0(tissues$Name[i],"Presence")] <- tData[match(Data$Gene, tData[,1]),3]
		}
		write.table(Data, file = file, quote = F, sep="\t", col.names = TRUE, row.names = TRUE)
		setwd(pwd)
	}
	return(Data)
}


prepare_Data <- function(file, ohn){
	Data <- read.table(file, h=F, sep = "\t", row.names=NULL)
	system_out <- system(paste("head -1 ", file, " | cut -f2,3,4,5,6,7,8,9,10,11,12 | sed 's/\\([A-Z]\\)[a-z]\\+_\\([a-z][a-z][a-z]\\)[A-Za-z0-9\\._]\\+/\\1\\2/g'"), intern=T)
	Header <- read.table(text=system_out, h=F, sep = "\t")
	colnames(Data) <- c("OG", as.character(unlist(lapply(Header[1,], as.character))))
	print(head(Data))
	Data$Ohnologs <- rep(FALSE, length(Data[,1]))
	Data$Ohnologs[which(Data$OG %in% ohn)] <- rep(TRUE, length(Data[which(Data$OG %in% ohn),1]))
	return(Data)
}


prepare_GenePairs <- function(nfile, blangd, drergd, tissues, tpmthresh){
	system_out <- system(paste("cat ",nfile ," | tail -n +2 | awk -F '\\t' '{split($2,bl,\" \"); split($3,dr,\" \"); for(i in bl){for(j in dr){print $1\"\\t\"bl[i]\"\\t\"dr[j]}}}'"), intern=T)
	GenePairs <- read.table(text=system_out, h=F, sep = "\t")
	colnames(GenePairs) <- c("OG", "Gene1", "Gene2")
	GenePairs <- GenePairs[which(GenePairs$Gene1 %in% blangd$Gene & GenePairs$Gene2 %in% drergd$Gene),]

	# Getting tissue gene expression data for each pair 
	for(i in c(1:length(tissues[,1]))){
		GenePairs[,paste0(tissues$Name[i],1)] <- blangd[match(GenePairs$Gene1, blangd$Gene), tissues$Blan[i]]
		GenePairs[,paste0(tissues$Name[i],2)] <- drergd[match(GenePairs$Gene2, drergd$Gene), paste0(tissues$Name[i],"TPM")]
		GenePairs[,paste0(tissues$Name[i],"2p")] <- drergd[match(GenePairs$Gene2, drergd$Gene), paste0(tissues$Name[i],"Presence")]
	}

	# Calculate sum of expressed tissues (domains) in each species and its difference
	GenePairs$SumDom1 <- unlist(lapply(c(1:length(GenePairs[,1])), function(x){sum(GenePairs[x,paste0(tissues$Name, 1)]>tpmthresh)}))
	GenePairs$SumDom2 <- unlist(lapply(c(1:length(GenePairs[,1])), function(x){sum(GenePairs[x,paste0(tissues$Name, 2)]>tpmthresh)}))
	GenePairs$DiffDom <- GenePairs$SumDom1-GenePairs$SumDom2

	# Calculate absolute difference between species expression profiles
	GenePairs$DiffAbs <- unlist(lapply(c(1:length(GenePairs[,1])), function(x){sum(abs((GenePairs[x,paste0(tissues$Name, 1)]>tpmthresh)-(GenePairs[x,paste0(tissues$Name, 2)]>tpmthresh)))}))	
	return(GenePairs)
}

prepare_CountsAV <- function(file, ohn){
	Counts <- read.table(file, h=F, sep = "\t", row.names=NULL)
	system_out <- system(paste("head -1 ", file, " | cut -f2,3,4,5,6,7,8,9,10,11,12 | sed 's/\\([A-Z]\\)[a-z]\\+_\\([a-z][a-z][a-z]\\)[A-Za-z0-9\\._]\\+/\\1\\2/g'"), intern=T)
	Header <- read.table(text=system_out, h=F, sep = "\t")
	colnames(Counts) <- c("OG", as.character(unlist(lapply(Header[1,], as.character))))
	Counts$Ohnologs <- rep(FALSE, length(Counts[,1]))
	Counts$Ohnologs[which(Counts$OG %in% ohn)] <- rep(TRUE, length(Counts[which(Counts$OG %in% ohn),1]))
	return(Counts)
}

