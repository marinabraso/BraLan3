

add_hist_logYaxis <- function(h, col){
	#for(i in c(1:length(h$counts))){
	#	polygon(c(h$breaks[i],h$breaks[i+1],h$breaks[i+1],h$breaks[i]),c(0,0,log(h$counts[i]),log(h$counts[i])), col=modif_alpha(col,.2), border=col)
	#}
	lines(h$mids, log(h$counts), col=col, lwd=4)
	points(h$mids, log(h$counts), col=col, pch=16, cex=3)
}

get_spearman_corr <- function(pairs, tissue){
	vec1 <- pairs[,paste0(tissue,1)]
	vec2 <- pairs[,paste0(tissue,2)]
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
