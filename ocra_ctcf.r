###### OCRA Peaks #######
load("~/pacoquita_bkp/ocra_ctcf_diffbind.Rda")

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(ReactomePA)
library(biomaRt)


txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

ocra.files <- read.delim2("~/ocra_peaks.txt", sep = "\t", header = F, stringsAsFactors = F)
names(ocra.files) <- c("path", "line", "mark")

# read.narrowPeak <- function(filename){
#   extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")
#   gr <- import(con = filename, format = "BED", extraCols = extraCols_narrowPeak)
#   return(gr)
# }

ctcf <- ocra.files$path[c(4,8,12,16,19,23,26,30,34,38,41,45,49,53,57,60,64,68)]
names(ctcf) <- ocra.files$line[c(4,8,12,16,19,23,26,30,34,38,41,45,49,53,57,60,64,68)]
peak <- readPeakFile(ctcf[1])
covplot(peak, weightCol="V5")

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(peak, windows=promoter)
tagMatrixList <- lapply(ctcf, getTagMatrix, windows=promoter)
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.95, resample=500, facet="row")
tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL)
peakAnnoList <- lapply(ctcf, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE)
plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList)


#### DiffBind CTCF ####
library(DiffBind)
manifest <- read.delim("~/Documents/ocra_ctcf.txt", sep = "\t", header = F, stringsAsFactors = F)
names(manifest) <- c("peaks", "cellLine", "mark", "bam_rep1", "bam_rep2")
aux <- manifest[, c(2,3,1,4)]
aux$Replicate <- "1"
names(aux)[4] <- "bamReads" 
aux2 <- manifest[, c(2,3,1,5)]
aux2$Replicate <- "2"
names(aux2)[4] <- "bamReads"
manifest <- rbind(aux,aux2)
manifest$Tissue <- c(rep("OvCancer", 10), rep("Normal", 5), rep("OvCancer", 13), rep("Normal", 5), rep("OvCancer", 3))
manifest$SubType <- c(rep("CCOC",3), "Endometriosis", rep("HGSOC", 5), "LGSOC", rep("FT",3), rep("OSE",2), rep("MUC", 3),
                      rep("CCOC",3), "Endometriosis", rep("HGSOC", 5), "LGSOC", rep("FT",3), rep("OSE",2), rep("MUC", 3))

manifest$PeakCaller <- "narrowPeak"
names(manifest) <- c("SampleID", "Factor1", "Peaks", "bamReads", "Replicate", "Tissue", "Factor", "PeakCaller")

ctcf.db <- dba(sampleSheet = manifest[1:18,], minOverlap = 1)
ctcf.db.c <- dba.count(ctcf.db, bParallel = T)
dba.plotPCA(ctcf.db.c, label = DBA_ID)
dba.plotHeatmap(ctcf.db.c, colScheme = "Greys")


### Making PCA without diffbind function - PRETTY PCA
aux <- dba.peakset(ctcf.db.c, bRetrieve = T) 
aux <- mcols(aux)
aux <- as.data.frame(aux)
pca1 <- prcomp(t(log2(aux+1)))
pca <- as.data.frame(pca1$x[, 1:2])
pca$anno <- manifest$Factor[1:18]
pca$shape <- manifest$Tissue[1:18]
pca$shape[4] <- "Normal"
rownames(pca)[which(!rownames(pca) %in% ocra.colors$CellLine)] <- c("HEY", "UWB1", "kuramochi", "IOSE4", "IOSE11", "GFTR230")

pdf("~/Documents/OCRA_Files/CTCF_Peaks_PCA_updated.pdf", 10,8)
ggplot(pca, aes(PC1, PC2, label = rownames(pca), color = ocra.colors[match(rownames(pca), ocra.colors$CellLine),]$Histotype)) + 
  geom_point(aes(shape=shape, size = 5)) + geom_text_repel(point.padding = 0.5, size = 5, segment.color = NA) + 
  theme_bw() +
  scale_color_manual(values = c("#D4126A", "#FD9933", "#97E17A", "#34CCCC", "#7F7F7F", "#6D69B3",
                                "#FF6699")) +
  coord_fixed() + 
  xlab(paste0("PC1 (", prettyNum(summary(pca1)$importance[2,1]*100, digits =2),"%)")) +
  ylab(paste0("PC2 (", prettyNum(summary(pca1)$importance[2,2]*100, digits =2),"%)")) +
  labs(color = "Histotype", shape = "Group") + guides(size = FALSE) +
  ggtitle("OCRA - CTCF Peaks")
dev.off()


## Normal vs Cancer 
contrast <- dba.contrast(ctcf.db.c, group1 = ctcf.db.c$masks$Normal, group2 = ctcf.db.c$masks$OvCancer)
contrast <- dba.analyze(contrast, bBlacklist = F)
plot(contrast, contrast=1)
ctcf.nVc <- dba.report(contrast)
save(ctcf.db, ctcf.db.c,ctcf.nVc, file="~/Documents/ocra_ctcf_diffbind.Rda")


tagMatrix <- getTagMatrix(ctcf.nVc, windows=promoter)
peakAnnoList <- annotatePeak(ctcf.nVc, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE)
plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList)
genes = (substr(as.data.frame(peakAnnoList)$transcriptId, 1, 15))
genes2 <- as.data.frame(peakAnnoList)
hmart=useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
bm <- getBM(attributes=c('ensembl_transcript_id', 'hgnc_symbol'),
                        filters = 'chromosome_name',
                        values = c(1:22,'X'),
                        mart = hmart)
bm <- bm[bm$ensembl_transcript_id %in% genes, ]
genes2$symbol <- bm$hgnc_symbol[match(bm$ensembl_transcript_id, substr(genes2$transcriptId, 1, 15))]




## FTSEC vs HGSOC 
## FTSEC vs HGSOC 
ctcf.db.c$masks$HGSOC[c("HeyA8", "OAW42")] <- FALSE
contrast2 <- dba.contrast(ctcf.db.c, group1 = ctcf.db.c$masks$FT, group2 = ctcf.db.c$masks$HGSOC)
contrast2 <- dba.analyze(contrast2)
pdf("~/Documents/OCRA_Files/DiffBind_CorHeatmap_CTCF_FTSECvsHGSOC.pdf", 10,8)
dba.plotHeatmap(contrast2, contrast=1, colScheme = "Greys",colSideCols = c("darkolivegreen4", "deepskyblue4"), 
                main="CTCF Peaks - FTSEC vs HGSOC")
dev.off()
ctcf.ftVc <- dba.report(contrast2)
ctcf.ftVc <- as.data.frame(ctcf.ftVc)
ctcf.ftVc <- ctcf.ftVc[ctcf.ftVc$Fold != 0 & ctcf.ftVc$FDR < 0.01, ]
# for remap divide by FT and HGSOC specific
write.table(ctcf.ftVc[ctcf.ftVc$Fold > 2, 1:3], sep= "\t", row.names = F, col.names = F, quote = F, 
            file = "~/Documents/OCRA_Files/ctcf_FTSEC.bed")
write.table(ctcf.ftVc[ctcf.ftVc$Fold < -2, 1:3], sep= "\t", row.names = F, col.names = F, quote = F, 
            file = "~/Documents/OCRA_Files/ctcf_HGSOC.bed")


#####################################################
k27.nVc$Mark <- "H3k27ac" 
k27.nVc <- as.data.frame(k27.nVc)
k4me1.nVc$Mark <- "H3k4me1"
k4me1.nVc <- as.data.frame(k4me1.nVc)
k4me3.nVc$Mark <- "H3k4me3"
k4me3.nVc <- as.data.frame(k4me3.nVc)
ctcf.nVc <- "CTCF"
ctcf.nVc <- as.data.frame(ctcf.nVc)
nVc.marks <- rbind(k27.nVc, k4me1.nVc, k4me3.nVc, ctcf.nVc)
























