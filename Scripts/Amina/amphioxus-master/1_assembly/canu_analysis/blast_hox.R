library('RColorBrewer')


shrink <- function(toShrink){
  shrinked <- data.frame()
  for(gene in unique(as.character(toShrink$gene))){
    gene_table <- toShrink[toShrink$gene == gene,]
    if(length(levels(factor(gene_table$contig))) > 1){
      contigs = length(levels(gene_table$contig))
      contig_alignments = c()
      for(con in 1:contigs){
        align_len <- sum((subset(gene_table, contig == levels(gene_table$contig)[con])$c_to) -
                           (subset(gene_table, contig == levels(gene_table$contig)[con])$c_from))
        contig_alignments <- c(contig_alignments, align_len)
      }
      gene_table <- subset(gene_table, contig == levels(gene_table$contig)[which.max(contig_alignments)])
    }
    form_sbset <- gene_table[which.min(gene_table$g_from),]
    to_sbset <- gene_table[which.min(abs(gene_table$g_size - gene_table$g_to)),]
    new_line <- data.frame('gene' = gene, 
                           'contig' = gene_table$contig[1], 
                           'from' = form_sbset$c_from[which.max(form_sbset$V3)], 
                           'to' = to_sbset$c_to[which.max(to_sbset$V3)], 
                           'size' = gene_table$g_size[1])
    shrinked <- rbind(shrinked, new_line)
  }
  return(shrinked)
}
# BLAST
setwd('../../../data/hox/')
genes_pb_localisation <- read.csv('pb_bflo_hox_genes.out', sep = '\t', header = F)
genes_il_localisation <- read.csv('il_bflo_hox_genes.out', sep = '\t', header = F)


colnames(genes_pb_localisation)[c(1,2,7:10,13)] <- c('gene','contig','g_from','g_to','c_from','c_to','g_size')
colnames(genes_il_localisation)[c(1,2,7:10,13)] <- c('gene','contig','g_from','g_to','c_from','c_to','g_size')

# manual_pb_pos_table <- data.frame(
#   gene = c('Hox1','Hox2','Hox3','Hox4','Hox5'),
#   from = c(1280414,1272393,1265464,1218514,1164911),
#   to = c(1281707,1272616,1267370,1219703,1168120),
#   length = c(867,618,1236,828,852)
# )


genes_pb_localisation <- shrink(genes_pb_localisation)
genes_il_localisation <-  shrink(genes_il_localisation)

# genes_pb_localisation <- genes_pb_localisation[-which(genes_pb_localisation$V4 < 60),]

# genes_pb_localisation$nc_from <- genes_pb_localisation$c_from - apply(data.frame(from = genes_pb_localisation$c_from, to = genes_pb_localisation$c_to), 1, min)
# genes_pb_localisation$nc_to <- genes_pb_localisation$c_to - apply(data.frame(from = genes_pb_localisation$c_from, to = genes_pb_localisation$c_to), 1, min)

colnames(genes_pb_localisation) <- c('gene','contig','c_from','c_to','q_size')
colnames(genes_il_localisation) <- c('gene','contig','c_from','c_to','q_size')

max(genes_il_localisation$c_to) - min(genes_il_localisation$c_from)
max(genes_pb_localisation$c_to) - min(genes_pb_localisation$c_from)

plot(c(), 
     xlim = c(min(c(genes_il_localisation$c_from, genes_pb_localisation$c_from)),max(c(genes_il_localisation$c_to, genes_pb_localisation$c_to))), 
     ylim = c(0,1), 
     xlab = 'position', 
     ylab = 'assembly', yaxt='n')

bpal <- brewer.pal(12, 'Paired')
bpal <- c(bpal, 'black', 'blue')

colour = 1
shift <- -0.14
# rnorm(1, mean = 0, sd = 0.02)

for(i in 1:(nrow(genes_pb_localisation) - 1)){
  PB <- c(genes_pb_localisation$c_from[i],genes_pb_localisation$c_to[i])
  arrows(PB[1], 0.25 + shift, PB[2] , 0.25 + shift, lwd = 5, col = bpal[colour])
  if(genes_pb_localisation$gene[i] != genes_pb_localisation$gene[i+1]){
    shift <- shift + 0.02
    colour = colour + 1
  }
}
PB <- c(genes_pb_localisation$c_from[i],genes_pb_localisation$c_to[i])
arrows(PB[1], 0.25 + shift, PB[2] , 0.25 + shift, lwd = 5, col = bpal[colour])

colour = 1
shift <- -0.14
# rnorm(1, mean = 0, sd = 0.02)

for(i in 1:(nrow(genes_il_localisation) - 1)){
  IL <- c(genes_il_localisation$c_from[i],genes_il_localisation$c_to[i])
  arrows(IL[1], 0.75 + shift, IL[2] , 0.75 + shift, lwd = 5, col = bpal[colour])
  if(genes_il_localisation$gene[i] != genes_il_localisation$gene[i+1]){
    shift <- shift + 0.02
    colour = colour + 1
  }
}
IL <- c(genes_il_localisation$c_from[i],genes_il_localisation$c_to[i])
arrows(IL[1], 0.75 + shift, IL[2] , 0.75 + shift, lwd = 5, col = bpal[colour])

axis(2, c(0.25,0.75), c('PB','IL'), las=1)
legend(1250000,1.05,as.character(unique(genes_il_localisation$gene)) , pch = 20, cex = 0.5, col = bpal)

max(genes_pb_localisation$c_to) - min(genes_pb_localisation$c_from)
max(genes_il_localisation$c_to) - min(genes_il_localisation$c_from)

PBdists <- genes_pb_localisation$c_to[1:13] - genes_pb_localisation$c_from[2:14]
Ildists <- genes_il_localisation$c_to[1:13] - genes_il_localisation$c_from[2:14]
PBdists
Ildists

ref_il <- read.csv('Illu_hox_Pb_asse_blast.out', sep = '\t', header = F)
alt_il <- read.csv('Ilhox_vs_altb.out', sep = '\t', header = F, skip = 5)
hapl_il <- read.csv('Ilhox_vs_haploidAssembly.out', sep = '\t', header = F)

alt_il <- alt_il[-which('# BLAST processed 1 queries' == alt_il$V1),]

summary(hapl_il)
summary(hapl_il[1:11443,])

colnames(ref_il)[c(1,2,7:10,13)] <- c('gene','contig','g_from','g_to','c_from','c_to','g_size')
colnames(alt_il)[c(1,2,7:10)] <- c('gene','contig','g_from','g_to','c_from','c_to')
colnames(hapl_il)[c(1,2,7:10)] <- c('gene','contig','g_from','g_to','c_from','c_to')

hapl_il_source <- hapl_il

summary(ref_il)

ref_il <- subset(ref_il,ref_il$contig == 'Sc0000087')
alt_il <- subset(alt_il,alt_il$contig == 'Sc0000087')

# plot(ref_il$V4)

max(ref_il$g_to) - min(ref_il$g_from)
max(ref_il$c_to) - min(ref_il$c_from)

# REF PB st
pdf('/Volumes/dump/pictures/pool/ngs/ampiho/HOX_blast.pdf')
  plot(c(), c(), ylim = c(min(ref_il[,c(7,8)]),max(ref_il[,c(7,8)])), 
       xlim = c(min(ref_il[,c(9,10)]),max(ref_il[,c(9,10)])), 
       ylab = 'Illumina HOX cluster', xlab = 'PacBio contig 87')
  for(i in 1:nrow(ref_il)){
    PB <- abs(ref_il[i,c(7,8)])
    lines(ref_il[i, c(9,10)], PB, lwd = 1)
  }
  
  for(i in 1:14){
    xs <- c(genes_il_localisation$from[i], genes_il_localisation$to[i])
    ys <- c(genes_pb_localisation$from[i],genes_pb_localisation$to[i])
    lines(xs, ys, col = 'red', lwd = 3)
  }
  
dev.off()


zoomPBToEdge <- function(x1, x2, y1, y2){
  plot(c(), c(), ylim = c(y1,y2), 
       xlim = c(x1,x2), 
       ylab = 'Illumina HOX cluster', xlab = 'PacBio contig 87')
  for(i in 1:nrow(ref_il)){
    PB <- abs(ref_il[i,c(7,8)])
    lines(ref_il[i, c(9,10)], PB, lwd = 1)
  }
  
  for(i in 1:14){
    xs <- c(genes_il_localisation$from[i], genes_il_localisation$to[i])
    ys <- c(genes_pb_localisation$from[i],genes_pb_localisation$to[i])
    lines(xs, ys, col = 'red', lwd = 3)
  }
}

zoomPBToEdge(0,2000000,0,2000000)
 zoomPBToEdge(500000,700000,0,1500000)
  zoomPBToEdge(610000,625000,250000,1500000)
   zoomPBToEdge(618080,618085,250000,1500000)
    zoomPBToEdge(618080,618085,494500,495500)
   zoomPBToEdge(621000,625000,250000,1500000)

# alt PB
pdf('/Volumes/dump/pictures/pool/ngs/ampiho/HOX_alt_blast.pdf')
plot(c(), c(), ylim = c(min(alt_il[,c(7,8)]),max(alt_il[,c(7,8)])), 
     xlim = c(min(alt_il[,c(9,10)]),max(alt_il[,c(9,10)])), 
     ylab = 'Illumina HOX cluster', xlab = 'PacBio alt contig 87')
for(i in 1:nrow(alt_il)){
  PB <- alt_il[i,c(7,8)]
  lines(alt_il[i, c(9,10)], PB, lwd = 1)
}
xs <- c(0, max(hapl_il$c_to))
ys <- c(min(genes_il_localisation$c_from),min(genes_il_localisation$c_from))
lines(xs, ys, col = 'red')
ys <- c(max(genes_il_localisation$c_to),max(genes_il_localisation$c_to))
lines(xs, ys, col = 'red')
dev.off()


#  hapl PB
hapl_il <- subset(hapl_il_source,hapl_il_source$contig == 'tig00010474')

pdf('/Volumes/dump/pictures/pool/ngs/ampiho/HOX_hapl_blast_c1.pdf')
plot(c(), c(), ylim = c(min(hapl_il[,c(7,8)]),max(hapl_il[,c(7,8)])), 
     xlim = c(min(hapl_il[,c(9,10)]),max(hapl_il[,c(9,10)])), 
     ylab = 'Illumina HOX cluster', xlab = 'PacBio Canu_e025 contig tig00010474')
for(i in 1:nrow(hapl_il)){
  PB <- hapl_il[i,c(7,8)]
  lines(hapl_il[i, c(9,10)], PB, lwd = 1)
}
xs <- c(0, max(hapl_il$c_to))
ys <- c(min(genes_il_localisation$c_from),min(genes_il_localisation$c_from))
lines(xs, ys, col = 'red')
ys <- c(max(genes_il_localisation$c_to),max(genes_il_localisation$c_to))
lines(xs, ys, col = 'red')
dev.off()

# summary(hapl_il_source[1:11444,])

pdf('/Volumes/dump/pictures/pool/ngs/ampiho/HOX_hapl_blast_c2.pdf')
hapl_il <- subset(hapl_il_source,hapl_il_source$contig == 'tig00003013')
plot(c(), c(), ylim = c(min(hapl_il[,c(7,8)]),max(hapl_il[,c(7,8)])), 
     xlim = c(min(hapl_il[,c(9,10)]),max(hapl_il[,c(9,10)])), 
     ylab = 'Illumina HOX cluster', xlab = 'PacBio Canu_e025 contig tig00003013')
for(i in 1:nrow(hapl_il)){
  PB <- hapl_il[i,c(7,8)]
  lines(hapl_il[i, c(9,10)], PB, lwd = 1)
}
xs <- c(0, max(hapl_il$c_to))
ys <- c(min(genes_il_localisation$c_from),min(genes_il_localisation$c_from))
lines(xs, ys, col = 'red')
ys <- c(max(genes_il_localisation$c_to),max(genes_il_localisation$c_to))
lines(xs, ys, col = 'red')
dev.off()

pdf('/Volumes/dump/pictures/pool/ngs/ampiho/HOX_hapl_blast_c3.pdf')
# summary(hapl_il_source[1:16444,])
hapl_il <- subset(hapl_il_source,hapl_il_source$contig == 'tig00003417')
plot(c(), c(), ylim = c(min(hapl_il[,c(7,8)]),max(hapl_il[,c(7,8)])), 
     xlim = c(min(hapl_il[,c(9,10)]),max(hapl_il[,c(9,10)])), 
     ylab = 'Illumina HOX cluster', xlab = 'PacBio Canu_e025 contig tig00003417')
for(i in 1:nrow(hapl_il)){
  PB <- hapl_il[i,c(7,8)]
  lines(hapl_il[i, c(9,10)], PB, lwd = 1)
}
xs <- c(0, max(hapl_il$c_to))
ys <- c(min(genes_il_localisation$from),min(genes_il_localisation$from))
lines(xs, ys, col = 'red')
ys <- c(max(genes_il_localisation$to),max(genes_il_localisation$to))
lines(xs, ys, col = 'red')
dev.off()

pdf('/Volumes/dump/pictures/pool/ngs/ampiho/HOX_hapl_blast_c4.pdf')
# summary(hapl_il_source[1:20000,])
hapl_il <- subset(hapl_il_source,hapl_il_source$contig == 'tig00001117')
plot(c(), c(), ylim = c(min(hapl_il[,c(7,8)]),max(hapl_il[,c(7,8)])), 
     xlim = c(min(hapl_il[,c(9,10)]),max(hapl_il[,c(9,10)])), 
     ylab = 'Illumina HOX cluster', xlab = 'PacBio Canu_e025 contig tig00001117')
for(i in 1:nrow(hapl_il)){
  PB <- hapl_il[i,c(7,8)]
  lines(hapl_il[i, c(9,10)], PB, lwd = 1)
}
xs <- c(0, max(hapl_il$c_to))
ys <- c(min(genes_il_localisation$c_from),min(genes_il_localisation$c_from))
lines(xs, ys, col = 'red')
ys <- c(max(genes_il_localisation$c_to),max(genes_il_localisation$c_to))
lines(xs, ys, col = 'red')
dev.off()

pdf('/Volumes/dump/pictures/pool/ngs/ampiho/HOX_hapl_blast_c5.pdf')
# summary(hapl_il_source[1:15000,])
hapl_il <- subset(hapl_il_source,hapl_il_source$contig == 'tig00003042')
plot(c(), c(), ylim = c(min(hapl_il[,c(7,8)]),max(hapl_il[,c(7,8)])), 
     xlim = c(min(hapl_il[,c(9,10)]),max(hapl_il[,c(9,10)])), 
     ylab = 'Illumina HOX cluster', xlab = 'PacBio Canu_e025 contig tig00003042')
for(i in 1:nrow(hapl_il)){
  PB <- hapl_il[i,c(7,8)]
  lines(hapl_il[i, c(9,10)], PB, lwd = 1)
}
xs <- c(0, max(hapl_il$c_to))
ys <- c(min(genes_il_localisation$from),min(genes_il_localisation$from))
lines(xs, ys, col = 'red')
ys <- c(max(genes_il_localisation$to),max(genes_il_localisation$to))
lines(xs, ys, col = 'red')

dev.off()


fileName <- "./scaffold37_amphio_blast.out"
conn <- file(fileName,open="r")
linn <-readLines(conn)

paraHOX <- read.csv('./scaffold37_amphio_blast.out', sep = '\t', header = F, skip = 5, col.names = strsplit(linn[4], ',')[[1]])
paraHOX <- paraHOX[!paraHOX$X..Fields..query.id == '# BLAST processed 1 queries',]
summary(paraHOX)

Sc81_ph <- paraHOX[paraHOX$X.subject.id == 'Sc0000081',]
sum(Sc81_ph$X.alignment.length)

pdf('/Volumes/dump/pictures/pool/ngs/ampiho/paraHox_inPBassembly.pdf')
plot(c(), c(), ylim = c(min(c(Sc81_ph$X.s..start, Sc81_ph$X.s..end)),max(c(Sc81_ph$X.s..start, Sc81_ph$X.s..end))), 
     xlim = c(min(c(Sc81_ph$X.q..start, Sc81_ph$X.q..start)),max(c(Sc81_ph$X.q..start, Sc81_ph$X.q..start))), 
     ylab = 'PacBio assembly, Sc0000081', xlab = 'Scaffold 37 - paraHOX')
for(i in 1:nrow(Sc81_ph)){
  PB <- Sc81_ph[i,c(7,8)]
  lines(PB, Sc81_ph[i, c(9,10)], lwd = 1)
}
dev.off()

# total length 1643025



hapl_il <- subset(hapl_il_source,hapl_il_source$contig == 'tig00010474')

plot(c(), c(), xlim = c(min(hapl_il[,c(7,8)]),max(hapl_il[,c(7,8)])), 
     ylim = c(min(hapl_il[,c(9,10)]),max(hapl_il[,c(9,10)])), 
     xlab = 'Illumina HOX cluster', ylab = 'PacBio Canu_e025 contig tig00010474')
for(i in 1:nrow(hapl_il)){
  PB <- hapl_il[i,c(7,8)]
  lines(PB,hapl_il[i, c(9,10)], lwd = 1)
}
xs <- c(0, max(hapl_il$c_to))
ys <- c(min(genes_il_localisation$c_from),min(genes_il_localisation$c_from))
lines(ys,xs, col = 'red')
ys <- c(max(genes_il_localisation$c_to),max(genes_il_localisation$c_to))
lines(ys,xs, col = 'red')

nrow(hapl_il)
nrow(hapl_il[hapl_il$g_from > 450000 & hapl_il$g_to < 550000 & hapl_il$c_to < 1500000,])
nrow(hapl_il[hapl_il$g_from > 450000 & hapl_il$g_to < 550000 & hapl_il$c_to < 1500000 & hapl_il$c_from > 1000000,])

for(i in 1:nrow(taget_seq)){
  PB <- taget_seq[i,c(7,8)]
  lines(PB,taget_seq[i, c(9,10)], lwd = 1, col = 'blue')
}

taget_seq <- hapl_il[hapl_il$g_from > 1450000 & hapl_il$g_to < 1550000,]
for(i in 1:nrow(taget_seq)){
  PB <- taget_seq[i,c(7,8)]
  lines(PB,taget_seq[i, c(9,10)], lwd = 1, col = 'green')
}

taget_seq <- hapl_il[hapl_il$g_from > 450000 & hapl_il$g_to < 550000 & hapl_il$c_to < 1250000 & hapl_il$c_from > 1150000,]
hist(taget_seq$c_to)
hist(taget_seq$g_to)

taget_seq <- hapl_il[hapl_il$g_from > 1450000 & hapl_il$g_to < 1550000,]
hist(taget_seq$g_to)

############
## seqinr ##
############

library('seqinr')

IL_hox <- unlist(read.fasta(file = 'seq/Illumina_hox.fa'))
PB_alt_hox <- unlist(read.fasta(file = 'seq/PB_alt_hox.fa'))

zoomPBToEdge(617500,618500,494500,495500)
dotPlot(PB_alt_hox[617500:618500], IL_hox[494500:495500], nmatch = 6, wsize = 10, wstep = 5)
dotPlot(PB_alt_hox[494500:495500], IL_hox[617500:618500], nmatch = 15, wsize = 20, wstep = 5)
length(PB_alt_hox)
length(IL_hox)

dotPlot(paste0(PB_alt_hox[578381:578405]),paste0(IL_hox[482558:482582]), nmatch = 2, wsize = 2, wstep = 1)
dotPlot(paste0(PB_alt_hox[578810:578948]),paste0(IL_hox[482943:483043]), nmatch = 3, wsize = 5, wstep = 1)

