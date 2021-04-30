setwd('/Volumes/dump/projects/PacBio/canu/')
backbone <- read.csv('mauve_out.xmfa.backbone', header = T, sep = '\t')

setwd('/Volumes/dump/projects/PacBio/Canu_report/')
backbone <- read.csv('PBrefxIL_mauve.backbone', header = T, sep = '\t')

clIL <- read.csv('/Volumes/dump/data/sequences/draft_genomes/Amphioxus/Bl71/lengths_illumina.csv', header = F)$V1
clPB <- read.csv('/Volumes/dump/data/sequences/draft_genomes/Amphioxus/canu_e025/amphio_contigs_lengths_ref_haplotype.csv', header = F)$V1



contig_87 <- backbone[which(backbone$seq0_leftend > cumclPB[87]  & backbone$seq0_leftend < cumclPB[88]),]
contig_87 <- contig_87[contig_87$seq1_leftend != 0 & contig_87$seq1_leftend < 20000000,]
plot(contig_87$seq1_leftend ~ contig_87$seq0_leftend)


sum(clPB)
max(backbone[,2],backbone[,1])

sum(clIL)
max(backbone[,4],backbone[,3])

# PB st
# pdf('/Volumes/dump/pictures/pool/ngs/ampiho/Illu_Canu_converved_seq_tl.pdf')
#   plot(c(), c(), xlim = c(0,max(backbone[,c(1,2)])), ylim = c(0,max(backbone[,c(3,4)])), xlab = 'PacBio assembly', ylab = 'Illumina assembly')
#   for(i in 1:nrow(backbone)){
#     if(any(backbone[i,] == 0)){
#       next
#     }
#     PB <- abs(backbone[i,c(1,2)])
#     lines(PB, backbone[i, c(3,4)], lwd = 0.1, col = mean(backbone[i,c(1,2)] < 0) + 1)
#   }
# dev.off()

notmapped_il <- backbone[which(backbone$seq1_leftend == 0),]
notmapped_pb <- backbone[which(backbone$seqo_leftend == 0),]
mapped <- backbone[which(backbone$seq0_leftend != 0 & backbone$seq1_leftend != 0),]

length(backbone$seq0_leftend != 0 & backbone$seq1_leftend != 0)

pb_total <- sum(abs(backbone[,2] - backbone[,1]))
il_total <- sum(abs(backbone[,4] - backbone[,3]))

# poswise <- rep(0, max(max(backbone[,c(1,2)])))
# for(i in 1:nrow(backbone)){
#   y1 <- abs(backbone[i,1])
#   y2 <- abs(backbone[i,2])
#   poswise[y1:y2] <- poswise[y1:y2]+1
# }
# 
# backbone[1:10,]
# 0 or 1

#pdf('/Volumes/dump/pictures/pool/ngs/ampiho/Illu_Canu_mauve_long.pdf')
# plot(c(), c(), ylim = c(0,max(backbone[,c(3,4)])), xlim = c(0,max(backbone[,c(1,2)])), ylab = 'PacBio assembly', xlab = 'Illumina assembly')
# for(i in 1:nrow(backbone)){
#   if(abs(backbone[i,c(1)] - backbone[i,c(2)]) < 2500 | abs(backbone[i,c(3)] - backbone[i,c(4)]) < 2500){
#     next
#   }
#   yc <- abs(backbone[i,c(1,2)])
#   lines(yc, backbone[i,c(3,4)], lwd = 0.4, col = mean(backbone[i,c(1,2)] < 0) + 1)
# }

pdf('/Volumes/dump/pictures/pool/ngs/ampiho/Illu_Canu_ref_long.pdf')
plot(c(), c(), xlim = c(0,max(backbone[,c(1,2)])), ylim = c(0,max(backbone[,c(3,4)])), xlab = 'PacBio assembly', ylab = 'Illumina assembly')
for(i in 1:nrow(backbone)){
  if(abs(backbone[i,c(1)] - backbone[i,c(2)]) < 2500 | abs(backbone[i,c(3)] - backbone[i,c(4)]) < 2500){
    next
  }
  PB <- abs(backbone[i,c(1,2)])
  lines(PB,backbone[i, c(3,4)], lwd = 0.1, col = mean(backbone[i,c(1,2)] < 0) + 1)
}

dev.off()

pdf('/Volumes/dump/pictures/pool/ngs/ampiho/Illu_Canu_ref_1000.pdf')
plot(c(), c(), ylim = c(0,max(backbone[,c(3,4)])), xlim = c(0,max(backbone[,c(1,2)])), ylab = 'Illumina assembly', xlab = 'PacBio assembly')
for(i in 1:nrow(backbone)){
  if(abs(backbone[i,c(1)] - backbone[i,c(2)]) < 1000 | abs(backbone[i,c(3)] - backbone[i,c(4)]) < 1000){
    next
  }
  yc <- abs(backbone[i,c(1,2)])
  lines(yc, backbone[i,c(3,4)], lwd = 0.4, col = mean(backbone[i,c(1,2)] < 0) + 1)
}
dev.off()

pdf('/Volumes/dump/pictures/pool/ngs/ampiho/Illu_Canu_ref_all.pdf')
plot(c(), c(), ylim = c(0,max(backbone[,c(3,4)])), xlim = c(0,max(backbone[,c(1,2)])), ylab = 'Illumina assembly', xlab = 'PacBio assembly')
for(i in 1:nrow(backbone)){
  if(any(backbone[i,] == 0)){
    next
  }
  yc <- abs(backbone[i,c(1,2)])
  lines(yc, backbone[i,c(3,4)], lwd = 0.4, col = mean(backbone[i,c(1,2)] < 0) + 1)
}
dev.off()

pdf('/Volumes/dump/pictures/pool/ngs/ampiho/Illu_Canu_ref_LITTLE_ZOOMED.pdf')
  
  
ILmax <- max(mapped[,c(3,4)]) / 2
PBmax <- max(mapped[,c(1,2)]) / 2

plot(c(), c(), ylim = c(0, ILmax), xlim = c(0,PBmax), ylab = 'Illumina assembly', xlab = 'PacBio assembly')
for(i in 1:nrow(mapped)){
  if(any(abs(mapped[i,c(1,2)]) > PBmax | abs(mapped[i,c(3,4)]) > ILmax)){
    next
  }
  yc <- abs(mapped[i,c(1,2)])
  lines(yc, mapped[i,c(3,4)], lwd = 0.4, col = mean(mapped[i,c(1,2)] < 0) + 1)
}

cumclIL <- cumsum(clIL)
cumclPB <- cumsum(clPB)

for(i in 1:78){
 lines(c(1, PBmax),c(cumclIL[i],cumclIL[i]), lwd = 0.3)
}

for(i in 1:88){
  lines(c(cumclPB[i],cumclPB[i]),c(1, ILmax), lwd = 0.3) 
}

dev.off()

pdf('/Volumes/dump/pictures/pool/ngs/ampiho/Illu_Canu_ref_ZOOMED.pdf')
ILmax <- max(mapped[,c(3,4)]) / 4
PBmax <- max(mapped[,c(1,2)]) / 4
BESTmapped <- mapped[mapped$seq0_leftend < PBmax & mapped$seq0_leftend < PBmax,]
BESTmapped <- BESTmapped[BESTmapped$seq1_leftend < ILmax & BESTmapped$seq1_rightend < ILmax,]
plot(c(), c(), ylim = c(0, ILmax), xlim = c(0,PBmax), ylab = 'Illumina assembly', xlab = 'PacBio assembly')
for(i in 1:nrow(BESTmapped)){
  yc <- abs(BESTmapped[i,c(1,2)])
  lines(yc, BESTmapped[i,c(3,4)], lwd = 0.4, col = mean(BESTmapped[i,c(1,2)] < 0) + 1)
}

cumclIL <- cumsum(clIL)
cumclPB <- cumsum(clPB)

for(i in 1:18){
  lines(c(1, PBmax),c(cumclIL[i],cumclIL[i]), lwd = 0.15)
}

for(i in 1:27){
  lines(c(cumclPB[i],cumclPB[i]),c(1, ILmax), lwd = 0.15) 
}

dev.off()


long_maps <- backbone[(abs(backbone[,c(1)] - backbone[,c(2)]) > 1000 & abs(backbone[,c(3)] - backbone[,c(4)]) > 1000),]

print(c("Pb mapped:", sum(abs(mapped[,4] - mapped[,3])) / pb_total))
print(c("Pb long mapped:", sum(abs(long_maps[,4] - long_maps[,3])) / pb_total))
print(c("Pb notmapped:", sum(abs(notmapped_il[,4] - notmapped_il[,3])) / pb_total))
print(c("Il mapped:", sum(abs(mapped[,2] - mapped[,1])) / il_total))
print(c("Il long mapped:", sum(abs(long_maps[,2] - long_maps[,1])) / il_total))
print(c("Il notmapped:", sum(abs(notmapped_pb[,2] - notmapped_pb[,1])) / il_total))

