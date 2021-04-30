graphicsPATH <- '/Volumes/dump/pictures/pool/ngs/ampiho/'
sbired <- rgb(255 / 255,36 / 255,31 / 255)
sbiredT <- rgb(255 / 255,36 / 255,31 / 255, alpha = 0.5)
illublueT <- rgb(6 / 255,243 / 255,255 / 255, alpha = 0.5)

getN50 <- function(cl, genomeSize = F, order = F){
  cl <- as.numeric(cl)
  if(genomeSize == F & order == F){
    return(sort(cl, decreasing = T)[max(which(cumsum(sort(cl, decreasing = T)) < (sum(cl) / 2))) + 1])
  }
  if(order == F){
    return(sort(cl, decreasing = T)[max(which(cumsum(sort(cl, decreasing = T, na.last = T)) < (genomeSize / 2)), na.rm = T) + 1])
  }
  return(max(which(cumsum(sort(cl, decreasing = T, na.last = T)) < (sum(cl) / 2)), na.rm = T) + 1)
}

minil <- read.csv('/Volumes/dump/projects/PacBio/miniassembly/12G_sam_file/lengths.csv', header = F)$V1

cl <- read.csv('/Volumes/dump/data/sequences/draft_genomes/Amphioxus/canu_e025/amphio_contig_lengths.csv', header = F)$V1
cl_ref <- read.csv('/Volumes/dump/data/sequences/draft_genomes/Amphioxus/canu_e025/amphio_contigs_lengths_ref_haplotype.csv', header = F)$V1
cl_alt <- read.csv('/Volumes/dump/data/sequences/draft_genomes/Amphioxus/canu_e025/amphio_contigs_lengths_alt_haplotype.csv', header = F)$V1

rrl <- read.csv('/Volumes/dump/data/sequences/raw_sequences/amhio/raw_reads_lengths.csv', header = F)
rl <- read.csv('/Volumes/dump/data/sequences/raw_sequences/amhio/amphio.correctedReads.length', header = F, sep = '\t')
# rlc <- read.csv('/Volumes/dump/data/sequences/draft_genomes/Amphioxus/canu_e025/amphio_contigs_lengths_ref_haplotype.csv', header = F)$V1

sl_il <- read.csv('/Volumes/dump/data/sequences/draft_genomes/Amphioxus/Bl71/lengths_illumina.csv', header = F)$V1

cumcl <- cumsum(cl_ref)

rrl <- rrl$V1
sum(as.numeric(rrl))
length(as.numeric(rrl))
80028931172 / 520000000
getN50(rrl)
sum(as.numeric(rrl[rrl > 13000])) / 520000000


pdf(paste(graphicsPATH,'raw_reads_histogram.pdf',sep = ''))
hist(rrl, xlab = 'length', main = '', col = sbired)
dev.off()

pdf(paste(graphicsPATH,'miniasm_histogram.pdf',sep = ''))
hist(log(minil,10), xlab = 'length', main = '', col = sbired)
dev.off()

PBmin
PBmax
min(which(cumcl > PBmin))
min(which(cumcl > PBmax))

pdf(paste(graphicsPATH,'canu_ctg_distr_e025.pdf',sep = ''))
  hist(log(cl,10), xlab = expression(paste(log[10],' contig length')), ylab = 'Number of contigs', main = 'Distribution of contig lengths')
dev.off()

pdf(paste(graphicsPATH,'canu_ctg_distr_ref.pdf',sep = ''))
  hist(log(cl_ref,10), xlab = expression(paste(log[10],' contig length')), ylab = 'Number of contigs', main = 'Distribution of contig lengths', breaks = 16)
dev.off()

bulish <- rgb(28, 44, 255, max = 255, alpha = 125, names = "blue50")
yellowish <- rgb(220, 160, 140, max = 255, alpha = 125, names = "yellow50")

pdf(paste(graphicsPATH,'haplomerger_efect_on_distr.pdf',sep = ''))
  hist(log(cl,10), xlab = expression(paste(log[10],' contig length')), ylab = 'Number of contigs', main = 'Distribution of contig lengths', freq = F, col = bulish)
  hist(log(cl_ref,10), xlab = expression(paste(log[10],' contig length')), ylab = 'Number of contigs', main = 'Distribution of contig lengths', breaks = 16, freq = F, add = T, col = yellowish)
  legend('topright', col = c(bulish,yellowish), pch = c(20,20), legend = c('diploid mix assembly',"haplotype assembly"))
dev.off()

sum(cl)
sum(cl_ref)

length(cl)
length(cl_ref)

max(cl)
max(cl_ref)
# computed N50
getN50(cl)
getN50(cl_ref)

pdf(paste(graphicsPATH,'raw_contigs_histogram.pdf',sep = ''))
hist(log(cl,10), xlab = expression(paste(log[10],' length')), main = '', col = sbired)
dev.off()

pdf(paste(graphicsPATH,'merged_contigs_histogram.pdf',sep = ''))
hist(log(cl_ref,10), xlab = expression(paste(log[10],' length')), main = '', col = sbired, breaks = 16)
dev.off()

getN50(sl_il)

# remove xlim, ylim for non zoomed
pdf(paste(graphicsPATH,'lengths_PB_vs_IL_ZOOMED.pdf',sep = ''))
hist(log(cl_ref,10), xlab = expression(paste(log[10],' length')), main = '', col = sbired, breaks = 16, xlim = c(6.3,8), ylim = c(0, 20))
hist(log(sl_il,10), col = illublueT, breaks = 16, add = T)
dev.off()



rlc <- rl$V3[!is.na(rl$V3)]
getN50(rlc)

pdf(paste(graphicsPATH,'corrected_reads_histogram.pdf',sep = ''))
  hist(rlc, xlab = 'length', main = '', col = sbired)
dev.off()

getN50(rrl$V1)
# theoretical N50 (approx. from known gneome size)
getN50(cl,520000000)
sort(cl, decreasing = T)[1:10]
length(cl)
length(cl_ref)
mean(cl)
median(cl)

fcl <- cl[cl > 100000]
length(fcl)
getN50(fcl)

fcl <- cl_ref[cl_ref > 2000000]
length(fcl)
getN50(fcl)

# readToTig showes which read is in which contig at what position
# tigInfo showes information about contigs...
rpos <- read.csv('/Volumes/dump/data/sequences/draft_genomes/Amphioxus/canu/amphio.layout.readToTig', header = F, sep = '\t')
cpos <- read.csv('/Volumes/dump/data/sequences/draft_genomes/Amphioxus/canu/amphio.layout.tinInfo', header = F, sep = '\t')

summary(cpos)
hist(log(cpos$V2,10))
hist(log(cpos$V4,10))
plot(cpos$V9)

# obviously there were more contigs before some filternig
mean(sort(cpos$V2, decreasing = T)[1:5505])
empirical_genomeSize <- sum(cpos$V2) / 2

# N50 for separated haplotypes
getN50(cl, sum(cpos$V2))
