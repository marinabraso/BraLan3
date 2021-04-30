setwd('/Volumes/dump/projects/amphio/hox/mapping_PBreads2ILhox/')

mr <- read.csv(header = F, file = './PBreads2ILhox_stats.csv', sep = ' ')
colnames(mr) <- c('read','flag','pos','legth')

lmr <- subset(mr, pos > 450000 & pos < 550000 )
rmr <- subset(mr, pos > 1450000 & pos < 1550000 )

dicontmr <- lmr[which(lmr$read %in% rmr$read),]
dicontmr <- rbind(dicontmr, rmr[which(rmr$read %in% lmr$read),])

dicontmr <- dicontmr[order(dicontmr$read),]

length(unique(dicontmr$read))
