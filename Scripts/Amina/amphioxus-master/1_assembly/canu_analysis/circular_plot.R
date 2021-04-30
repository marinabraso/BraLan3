checkFlag <- function(flags, flagIndex){
  flagToCheck <- 2^c(0:12)[flagIndex+1]
  for(i in 2^c(1:12)){
    if(i == flagToCheck){
      break;
    }
    flags <- flags - (flags %% i)
  }
  return((flags %% flagToCheck) != 0)
}

sam2LinkTable <- function(files){
  
  out_tab <- data.frame(ch1 = character(0), s1 = numeric(0), e1 = numeric(0), 
                        ch2 = character(0), s2 = numeric(0), e2 = numeric(0))
  
  for(i in 1:length(files)){
    samf <- read.csv(file = files[i], header = F, sep = '\t', skip = 3, comment.char = 'Y')
    samf$V10 <- c()
    samf$V11 <- c()
    
    sam_mp1 <- samf[seq(1, nrow(samf), by = 2),]
    sam_mp2 <- samf[seq(2, nrow(samf), by = 2),]
    
    # compute where mate pairs are RF or FR (in sets)
    rev1 <- checkFlag(sam_mp1$V2,5) & !checkFlag(sam_mp2$V2,5)
    rev2 <- !checkFlag(sam_mp1$V2,5) & checkFlag(sam_mp2$V2,5)
    
    # take only those where are going in opposite direction
    solid_mp <- (rev1 | rev2) & (abs(sam_mp1$V4 - sam_mp2$V4) > 1000)
    sam_mp1 <- sam_mp1[solid_mp,]
    sam_mp2 <- sam_mp2[solid_mp,]
    
    # recompute logical vectors for filtered tables
    rev1 <- checkFlag(sam_mp1$V2,5) & !checkFlag(sam_mp2$V2,5)
    rev2 <- !checkFlag(sam_mp1$V2,5) & checkFlag(sam_mp2$V2,5)
    
    # init vectors for coordinates
    s1 <- rep(0, nrow(sam_mp1))
    e1 <- rep(0, nrow(sam_mp1))
    s2 <- rep(0, nrow(sam_mp2))
    e2 <- rep(0, nrow(sam_mp2))
    
    sam_mp1$length <- unlist(lapply(lapply(strsplit(as.character(sam_mp1$V6), "[[:alpha:]]"), as.numeric), sum))
    sam_mp2$length <- unlist(lapply(lapply(strsplit(as.character(sam_mp2$V6), "[[:alpha:]]"), as.numeric), sum))
    
    s1[!rev1] <- sam_mp1$V4[!rev1]
    s1[rev1] <- sam_mp1$V4[rev1] + sam_mp1$length[rev1]
    e1[!rev1] <- sam_mp1$V4[!rev1] + sam_mp1$length[!rev1]
    e1[rev1] <- sam_mp1$V4[rev1]
    
    s2[!rev2] <- sam_mp2$V4[!rev2]
    s2[rev2] <- sam_mp2$V4[rev2] + sam_mp2$length[rev2]
    e2[!rev2] <- sam_mp2$V4[!rev2] + sam_mp2$length[!rev2]
    e2[rev2] <- sam_mp2$V4[rev2]
    
    out_tab <- rbind(out_tab,data.frame(ch1 = sam_mp1$V3, s1 = s1, e1 = e1, ch2 = sam_mp2$V3, s2 = s2, e2 = e2))
  }

  return(out_tab)
}

il_libs <- read.csv('/Volumes/dump/data/mapping/amphio/hox/il2pb/libraries.txt', header = F)$V1
il_libs <- as.integer(substr(il_libs, 30, 30)) * 1000

il2il_files <- paste0('/Volumes/dump/data/mapping/amphio/hox/il2il/', dir('/Volumes/dump/data/mapping/amphio/hox/il2il/', pattern = "lib"))
il2pb_files <- paste0('/Volumes/dump/data/mapping/amphio/hox/il2pb/', dir('/Volumes/dump/data/mapping/amphio/hox/il2pb/', pattern = "lib"))

out_tab <- sam2LinkTable(il2il_files)
dim(out_tab)
out_tab <- format(out_tab, scientific = FALSE)
write.table(x = out_tab, file = '../all_links_il2il', quote = F, row.names = F, col.names = F, sep = ' ')

files <- il2pb_files[1:6]
# pb_out_tab <- sam2LinkTable(files) function is not working, but reruning the same code without function wrapper makes it work. No idea what is wrong.
out_tab <- format(out_tab, scientific = FALSE)
write.table(x = out_tab, file = '../all_links_il2pn', quote = F, row.names = F, col.names = F, sep = ' ')

# setwd('/Volumes/dump/projects/PacBio/hox/')
# genes_pb_localisation <- read.csv('pb_bflo_hox_genes.out', sep = '\t', header = F)
genes_il_localisation <- read.csv('/Volumes/dump/projects/amphio/hox/il_bflo_hox_genes.out', sep = '\t', header = F)
genes_pb_localisation <- read.csv('/Volumes/dump/projects/amphio/hox/pb_bflo_hox_genes.out', sep = '\t', header = F)

colnames(genes_il_localisation)[c(1,2,7:10,13)] <- c('gene','contig','g_from','g_to','c_from','c_to','g_size')
colnames(genes_pb_localisation)[c(1,2,7:10,13)] <- c('gene','contig','g_from','g_to','c_from','c_to','g_size')

genes_il_localisation <-  shrink(genes_il_localisation) # functino shrink is in blast_hox.R
genes_pb_localisation <-  shrink(genes_pb_localisation) # functino shrink is in blast_hox.R

colnames(genes_il_localisation) <- c('gene','contig','c_from','c_to','q_size')
colnames(genes_pb_localisation) <- c('gene','contig','c_from','c_to','q_size')

write.table(x = genes_il_localisation[,c(2,3,4)], file = '../hox_genes_on_il', quote = F, row.names = F, col.names = F, sep = ' ')
write.table(x = genes_pb_localisation[,c(2,3,4)], file = '../hox_genes_on_pb', quote = F, row.names = F, col.names = F, sep = ' ')

# sam_mp1_rev <- sam_mp1[rev1,1:5]
# sam_mp2_fw <- sam_mp2[rev1,1:5]
# 
# head(sam_mp1_rev[sam_mp1_rev[,4] > sam_mp2_fw[,4],])
# head(sam_mp2_fw[sam_mp1_rev[,4] > sam_mp2_fw[,4],])
# sum(sam_mp1[rev1,4] > sam_mp2[sam_mp2$V1 %in% sam_mp1[rev1,1],4])
# 
# np <- checkFlag(sam_mp1$V2,2)
# hist(log10(abs(sam_mp1$V4[!np] - sam_mp2$V4[!np])), breaks = 40, main = 'Histogram of insert sizes', xlab = 'log10( is )')
# hist(log10(abs(sam_mp1$V4[np] - sam_mp2$V4[np])), add = T , col = 'red')
# hist(log10(abs(sam_mp1$V4[!np] - sam_mp2$V4[!np])), breaks = 40, add = T, col = 'blue')
# legend('topright', legend = c('T','F'), pch = 20, col = c('red','blue'), title = 'flag bit 2')
# sum(np)
# 
# np <- !checkFlag(sam_mp1$V2,2) & ((abs(sam_mp1$V4 - sam_mp2$V4)) < 1000)
# sum(np)
# 
# head(sam_mp1[np,1:5])
# head(sam_mp2[np,1:5])
# 
# hist(log10(sam_mp1_rev[sam_mp1_rev[,4] > sam_mp2_fw[,4],4] - sam_mp2_fw[sam_mp1_rev[,4] > sam_mp2_fw[,4],4]))
# 
# sum(checkFlag(sam_mp2[sam_mp2$V1 %in% sam_mp1[rev1,1],2],5))





#as.integer(intToBits(129))

#is_same <- sam_mp1$V4[checkFlag(sam_mp1$V2,5) == checkFlag(sam_mp2$V2,5)] - sam_mp2$V4[checkFlag(sam_mp1$V2,5) == checkFlag(sam_mp2$V2,5)]

#write.table(x = chlinks, file = '../fake_links.txt', quote = F, row.names = F, col.names = F, sep = ' ')