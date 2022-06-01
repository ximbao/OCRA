### Hichip summary for OCRA supp table ###


## prep the objects in GRanges mode
ft246.hichip <- fread("~/Documents/Hichip_calls/FitHiChIP.interactions_FT246_FitHiC_Q0.01.bed")
kura.hichip <- fread("~/Documents/Hichip_calls/FitHiChIP.interactions_FitHiC_Q0.01_KURA_NARROWP_20KB2MB.bed")
ft33.hichip <- fread("~/FitHiChIP/FT33/FitHiChIP_Peak2ALL_b5000_L20000_U2000000/P2PBckgr_1/Coverage_Bias/FitHiC_BiasCorr/FitHiChIP_FT33.interactions_FitHiC_Q0.01.bed")
uwb.hichip <- fread("~/Documents/Hichip_calls/FitHiChIP.interactions_FitHiC_Q0.01_UWB1289_NARROWP_20KB_2MB.bed")

ft246.hichip$s1 <- as.numeric(ft246.hichip$s1)
ft246.hichip$s2 <- as.numeric(ft246.hichip$s2)
kura.hichip$V2 <- as.numeric(kura.hichip$V2)
kura.hichip$V5 <- as.numeric(kura.hichip$V5)
ft33.hichip$s1 <- as.numeric(ft33.hichip$s1)
ft33.hichip$s2 <- as.numeric(ft33.hichip$s2)
uwb.hichip$s1 <- as.numeric(uwb.hichip$s1)
uwb.hichip$s2 <- as.numeric(uwb.hichip$s2)


colnames(ft246.hichip)[1:3] <- c("chr", "start", "end")
colnames(kura.hichip)[1:3] <- c("chr", "start", "end")
colnames(ft33.hichip)[1:3] <- c("chr", "start", "end")
colnames(uwb.hichip)[1:3] <- c("chr", "start", "end")
colnames(ft246.hichip)[4:6] <- c("chr", "start", "end")
colnames(kura.hichip)[4:6] <- c("chr", "start", "end")
colnames(ft33.hichip)[4:6] <- c("chr", "start", "end")
colnames(uwb.hichip)[4:6] <- c("chr", "start", "end")

ft246.hichip1 <- makeGRangesFromDataFrame(na.omit(ft246.hichip[, 1:3]))
kura.hichip1 <- makeGRangesFromDataFrame(na.omit(kura.hichip[, 1:3]))
ft33.hichip1 <- makeGRangesFromDataFrame(na.omit(ft33.hichip[, 1:3]))
uwb.hichip1 <- makeGRangesFromDataFrame(na.omit(uwb.hichip[, 1:3]))

ft246.hichip2 <- makeGRangesFromDataFrame(na.omit(ft246.hichip[, 4:6]))
kura.hichip2 <- makeGRangesFromDataFrame(na.omit(kura.hichip[, 4:6]))
ft33.hichip2 <- makeGRangesFromDataFrame(na.omit(ft33.hichip[, 4:6]))
uwb.hichip2 <- makeGRangesFromDataFrame(na.omit(uwb.hichip[, 4:6]))


## Number of loops in genes 
a <- queryHits(findOverlaps(c(uwb.hichip1), genes(TxDb.Hsapiens.UCSC.hg19.knownGene)))
b <- queryHits(findOverlaps(c(uwb.hichip2), genes(TxDb.Hsapiens.UCSC.hg19.knownGene)))
length(unique(c(a,b)))

a <- queryHits(findOverlaps(c(kura.hichip1), genes(TxDb.Hsapiens.UCSC.hg19.knownGene)))
b <- queryHits(findOverlaps(c(kura.hichip2), genes(TxDb.Hsapiens.UCSC.hg19.knownGene)))
length(unique(c(a,b)))

a <- queryHits(findOverlaps(c(ft246.hichip1), genes(TxDb.Hsapiens.UCSC.hg19.knownGene)))
b <- queryHits(findOverlaps(c(ft246.hichip2), genes(TxDb.Hsapiens.UCSC.hg19.knownGene)))
length(unique(c(a,b)))

a <- queryHits(findOverlaps(c(ft33.hichip1), genes(TxDb.Hsapiens.UCSC.hg19.knownGene)))
b <- queryHits(findOverlaps(c(ft33.hichip2), genes(TxDb.Hsapiens.UCSC.hg19.knownGene)))
length(unique(c(a,b)))

## Number of loops in promoters 
promot <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
start(promot) <- start(promot)-500
end(promot) <- start(promot)+500

a <- queryHits(findOverlaps(c(uwb.hichip1), promot))
b <- queryHits(findOverlaps(c(uwb.hichip2), promot))
length(unique(c(a,b)))

a <- queryHits(findOverlaps(c(kura.hichip1), promot))
b <- queryHits(findOverlaps(c(kura.hichip2), promot))
length(unique(c(a,b)))

a <- queryHits(findOverlaps(c(ft246.hichip1), promot))
b <- queryHits(findOverlaps(c(ft246.hichip2), promot))
length(unique(c(a,b)))

a <- queryHits(findOverlaps(c(ft33.hichip1), promot))
b <- queryHits(findOverlaps(c(ft33.hichip2), promot))
length(unique(c(a,b)))



