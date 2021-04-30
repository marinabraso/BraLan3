


setwd('/Volumes/dump/data/sequences/draft_genomes/Amphioxus/HOX_Blan/')
backbone <- read.csv('HOX_pb.backbone', header = T, sep = '\t')
clPB <- read.csv('/Volumes/dump/data/sequences/draft_genomes/Amphioxus/canu_e025/amphio_contigs_lengths_ref_haplotype.csv', header = F)$V1

sum(clPB)
max(backbone[,2],backbone[,1])
max(backbone[,4],backbone[,3])

mapped <- backbone[which(backbone$seq0_rightend != 0 & backbone$seq1_rightend != 0),]
NONmapped <- backbone[which(backbone$seq0_rightend == 0 & backbone$seq1_rightend != 0),]

dim(NONmapped)
summary(NONmapped)
sum(NONmapped$seq1_rightend - NONmapped$seq1_leftend)

ALLmapped <- backbone[which(backbone$seq0_rightend != 0 & backbone$seq1_rightend != 0),]
LEFTmapped <- backbone[which(backbone$seq0_rightend != 0 & backbone$seq1_rightend != 0 & abs(backbone$seq0_rightend) < 200000000),]
RIGHTmapped <- backbone[which(backbone$seq0_rightend != 0 & backbone$seq1_rightend != 0 & abs(backbone$seq0_rightend) > 230000000 & abs(backbone$seq0_rightend) < 260000000),]

mapped <- ALLmapped
PBmin <- min(abs(mapped[,c(1,2)]))
PBmax <- max(abs(mapped[,c(1,2)]))
HOXmin <- min(abs(mapped[,c(3,4)]))
HOXmax <- max(abs(mapped[,c(3,4)]))

# pdf('/Volumes/dump/pictures/pool/ngs/ampiho/Hox_mapping.pdf')
# plot(c(), c(), ylim = c(HOXmin,HOXmax), xlim = c(PBmin,PBmax), ylab = 'HOX of Illumina assembly', xlab = 'PacBio assembly')
# for(i in 1:nrow(mapped)){
#   if(mapped[i,c(1)] - mapped[i,c(2)] == 0 | abs(mapped[i,c(3)] - mapped[i,c(4)]) == 0){
#     next
#   }
#   PB <- abs(mapped[i,c(1,2)])
#   lines(PB, mapped[i,c(3,4)], lwd = 0.4, col = mean(mapped[i,c(1,2)] < 0) + 1)
# }
# dev.off()

# mapped <- LEFTmapped[-c(1:4),]
# PBmin <- min(abs(mapped[,c(1,2)]))
# PBmax <- max(abs(mapped[,c(1,2)]))
# HOXmin <- min(abs(mapped[,c(3,4)]))
# HOXmax <- max(abs(mapped[,c(3,4)]))
# PBmax - PBmin
# 
# pdf('/Volumes/dump/pictures/pool/ngs/ampiho/Hox_mapping_leftZOOM_pos.pdf')
# plot(c(), c(), ylim = c(HOXmin,HOXmax), xlim = c(PBmin,PBmax), ylab = 'HOX of Illumina assembly', xlab = 'PacBio assembly')
# for(i in 1:nrow(mapped)){
#   if(mapped[i,c(1)] - mapped[i,c(2)] == 0 | abs(mapped[i,c(3)] - mapped[i,c(4)]) == 0){
#     next
#   }
#   PB <- mapped[i,c(1,2)]
#   lines(PB, mapped[i,c(3,4)], lwd = 0.4, col = mean(mapped[i,c(1,2)] < 0) + 1)
# }
# dev.off()

mapped <- LEFTmapped[-c(1:4),]
PBmin <- min(abs(mapped[,c(1,2)]))
PBmax <- max(abs(mapped[,c(1,2)]))
HOXmin <- min(abs(mapped[,c(3,4)]))
HOXmax <- max(abs(mapped[,c(3,4)]))
# PBmax - PBmin
# 
# pdf('/Volumes/dump/pictures/pool/ngs/ampiho/Hox_mapping_leftZOOM_neg.pdf')
plot(c(), c(), ylim = c(HOXmin,HOXmax), xlim = c(PBmin,PBmax), ylab = 'HOX of Illumina assembly', xlab = 'PacBio assembly')
for(i in 1:nrow(mapped)){
  if(mapped[i,c(1)] - mapped[i,c(2)] == 0 | abs(mapped[i,c(3)] - mapped[i,c(4)]) == 0){
    next
  }
  PB <- mapped[i,c(1,2)]
  lines(PB, mapped[i,c(3,4)], lwd = 0.4, col = mean(mapped[i,c(1,2)] < 0) + 1)
}
# dev.off()

mapped <- RIGHTmapped[-c(1:4),]
PBmin <- min(abs(mapped[,c(1,2)]))
PBmax <- max(abs(mapped[,c(1,2)]))
HOXmin <- min(abs(mapped[,c(3,4)]))
HOXmax <- max(abs(mapped[,c(3,4)]))
PBmax - PBmin

# pdf('/Volumes/dump/pictures/pool/ngs/ampiho/Hox_mapping_rightZOOM_neg.pdf')
plot(c(), c(), ylim = c(HOXmin,HOXmax), xlim = c(PBmin,PBmax), ylab = 'HOX of Illumina assembly', xlab = 'PacBio assembly')
for(i in 1:nrow(mapped)){
  if(mapped[i,c(1)] - mapped[i,c(2)] == 0 | abs(mapped[i,c(3)] - mapped[i,c(4)]) == 0){
    next
  }
  PB <- abs(mapped[i,c(1,2)])
  lines(PB, mapped[i,c(3,4)], lwd = 0.4, col = mean(mapped[i,c(1,2)] < 0) + 1)
}
# dev.off()

# mapped <- RIGHTmapped
# PBmin <- min(abs(mapped[,c(1,2)]))
# PBmax <- max(abs(mapped[,c(1,2)]))
# HOXmin <- min(abs(mapped[,c(3,4)]))
# HOXmax <- max(abs(mapped[,c(3,4)]))
# PBmax - PBmin
# 
# pdf('/Volumes/dump/pictures/pool/ngs/ampiho/Hox_mapping_rightZOOM.pdf')
# plot(c(), c(), ylim = c(HOXmin,HOXmax), xlim = c(PBmin,PBmax), ylab = 'HOX of Illumina assembly', xlab = 'PacBio assembly')
# for(i in 1:nrow(mapped)){
#   if(mapped[i,c(1)] - mapped[i,c(2)] == 0 | abs(mapped[i,c(3)] - mapped[i,c(4)]) == 0){
#     next
#   }
#   PB <- abs(mapped[i,c(1,2)])
#   lines(PB, mapped[i,c(3,4)], lwd = 0.4, col = mean(mapped[i,c(1,2)] < 0) + 1)
# }
# dev.off()
