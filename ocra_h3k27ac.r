#### OCRA H3k27ac peaks ####\

load("~/pacoquita_bkp/ocra_h3k27ac_diffbind.Rda")

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(ReactomePA)
library(biomaRt)
library(DiffBind)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


ocra.27ac   <- read.delim("~/pacoquita_bkp/ocra_h3k27ac_peaks.txt", header = F, sep = "\t", stringsAsFactors = F)
ocra.27ac <- ocra.27ac[, -5]
names(ocra.27ac) <- c("SampleID", "mark", "Peaks", "bamReads")
ocra.27ac$Replicate <- "1"
ocra.27ac$Tissue <- c(rep("OvCancer", 10), rep("Normal", 4), rep("OvCancer", 2), "Normal", "OvCancer")
ocra.27ac$Factor <- c(rep("CCOC",3), "Endometrioid", rep("HGSOC",5), "LGSOC", rep("FTSEC",2),
                       rep("OSEC",2), rep("Mucinous", 2), "FTSEC", "Mucinous")
ocra.27ac$PeakCaller <- "narrowPeak"

ocra.27ac$Tissue[4] <- "Normal"

k27.db <- dba(sampleSheet = ocra.27ac, minOverlap = 1)
k27.db.c <- dba.count(k27.db, bParallel = T)
dba.plotPCA(k27.db.c, label = DBA_ID)
dba.plotHeatmap(k27.db.c, colScheme = "Greys")

## Normal vs Cancer 
contrast <- dba.contrast(k27.db.c, group1 = k27.db.c$masks$Normal, group2 = k27.db.c$masks$OvCancer)
contrast <- dba.analyze(contrast)
plot(contrast, contrast=1)
k27.nVc <- dba.report(contrast)
# save(k27.db, k27.db.c, k27.nVc, file = "~/Documents/ocra_h3k27ac_diffbind.Rda")

tagMatrix <- getTagMatrix(k27.nVc, windows=promoter)
peakAnnoList <- annotatePeak(k27.nVc, TxDb=txdb,
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

### Making PCA without diffbind function - PRETTY PCA
aux <- dba.peakset(k27.db.c, bRetrieve = T) 
aux <- mcols(aux)
aux <- as.data.frame(aux)
pca1 <- prcomp(t(log2(aux+1)))
pca <- as.data.frame(pca1$x[, 1:3])
pca$anno <- ocra.27ac$Factor
pca$shape <- ocra.27ac$Tissue

pdf("~/Documents/OCRA_Files/OCRA_h3k27ac_PCA_updated.pdf", 10,8)
ggplot(pca, aes(PC1, PC2, label = rownames(pca), color = shape)) + 
  geom_point(aes(size= 5,shape=anno)) + geom_text_repel(point.padding = 0.5, size = 5, segment.color = NA) + 
  theme_bw() +
  scale_color_manual(values = c("darkolivegreen3", "deepskyblue3")) +
  scale_shape_manual(values=c(15,18,22,16,3,2,4,2)) +
  coord_fixed() + 
  xlab(paste0("PC1 (", prettyNum(summary(pca1)$importance[2,1]*100, digits =2),"%)")) +
  ylab(paste0("PC2 (", prettyNum(summary(pca1)$importance[2,2]*100, digits =2),"%)")) +
  labs(color = "Histotype", shape = "Group") + guides(size = FALSE) +
  ggtitle("OCRA - H3K27ac Peaks")
dev.off()

library(plot3D)
plot3D::scatter3D(pca1$x[, 1], pca1$x[, 2], pca1$x[, 3], colvar = NULL, pch = 20, bty="g",
                  cex = 2, phi=20, col = kelly.colours[2:dim(pca1$x)[1]])
text3D(pca1$x[, 1], pca1$x[, 2], pca1$x[, 3], labels=rownames(pca1$x), add = T, cex = 0.6)


## FTSEC vs HGSOC 
contrast2 <- dba.contrast(k27.db.c, group1 = k27.db.c$masks$FT, group2 = k27.db.c$masks$HGSOC)
contrast2 <- dba.analyze(contrast2, bBlacklist = F)
pdf("~/Documents/OCRA_Files/DiffBind_CorHeatmap_H3k27ac_FTSECvsHGSOC_updated.pdf", 10,8)
dba.plotHeatmap(contrast2, contrast=1, colScheme = "Greys",colSideCols = c("darkolivegreen3", "deepskyblue3"), 
                main="H3k27ac Peaks - FTSEC vs HGSOC")
dev.off()

pdf("~/Documents/OCRA_Files/volc_k27ac_FTSEC_vs_HGSOC.pdf", 8,10)
dba.plotVolcano(contrast2)
dev.off()

k27.ftVc <- dba.report(contrast2)
k27.ftVc <- as.data.frame(k27.ftVc)
k27.ftVc <- k27.ftVc[k27.ftVc$Fold != 0 & k27.ftVc$FDR < 0.01, ]
# for remap divide by FT and HGSOC specific
write.table(k27.ftVc[k27.ftVc$Fold > 2, 1:3], sep= "\t", row.names = F, col.names = F, quote = F, 
            file = "~/Documents/OCRA_Files/k27ac_DiffBind_FTSEC.bed")
write.table(k27.ftVc[k27.ftVc$Fold < -2, 1:3], sep= "\t", row.names = F, col.names = F, quote = F, 
            file = "~/Documents/OCRA_Files/k27ac_DiffBindHGSOC.bed")

### Plot remap enrichment ###
enrichment.diffbind.ft <- read.csv("~/Documents/OCRA_Files/Result_enrichmentDiffBind_FTSEC_csv.csv")
pdf("~/Documents/OCRA_Files/TFBS_H3k27ac_DiffBind_FTSEC.pdf", 8, 10)
enrichmentDotPlot(enrichment.diffbind.ft, top=20, main = "TFBS - FTSEC Peaks differentially bound FTSEC vs HGSOC")
dev.off()

enrichment.diffbind.hgsoc <- read.csv("~/Documents/OCRA_Files/Result_enrichmentDiffBind_HGSOC_csv.csv")
pdf("~/Documents/OCRA_Files/TFBS_H3k27ac_DiffBind_HGSOC.pdf", 8, 10)
enrichmentDotPlot(enrichment.diffbind.hgsoc, top=20, main = "TFBS - HGSOC Peaks differentially bound FTSEC vs HGSOC")
dev.off()



### Get genes/nearest gene on diffbound peaks ###
bm <- getBM(attributes=c('entrezgene_id', 'hgnc_symbol'),
            filters = 'chromosome_name',
            values = c(1:22,'X'),
            mart = ensembl)

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
k27.ftVc <- dba.report(contrast2)
## make a table with results from diffbind + gene in the diffbound peak ##
k27.ft.diffbind <- k27.ftVc[k27.ftVc$Fold > 2 & k27.ftVc$FDR < 0.01]
ovp <- findOverlaps(genes(TxDb.Hsapiens.UCSC.hg38.knownGene),
                 k27.ftVc[k27.ftVc$Fold > 2 & k27.ftVc$FDR < 0.01]) 
k27.ft.diffbind <- k27.ft.diffbind[(subjectHits(ovp))]
k27.ft.diffbind$Gene <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)$gene_id[queryHits(ovp)]
aux <- bm[bm$entrezgene_id %in% k27.ft.diffbind$Gene,]
aux <- aux[match(k27.ft.diffbind$Gene, aux$entrezgene_id),]
k27.ft.diffbind$Gene <- aux$hgnc_symbol
k27.ft.diffbind <- as.data.frame(k27.ft.diffbind, row.names = NULL)
k27.ft.diffbind <- k27.ft.diffbind[order(k27.ft.diffbind$Fold, decreasing = T), ]
write.table(k27.ft.diffbind, quote = F, row.names = F, sep = "\t", 
            file = "~/Documents/OCRA_Files/H3k27ac_DiffBind_FTSEC_table.txt")

k27.hgsoc.diffbind <- k27.ftVc[k27.ftVc$Fold < -2 & k27.ftVc$FDR < 0.01]
ovp2 <- findOverlaps(genes(TxDb.Hsapiens.UCSC.hg38.knownGene),
                    k27.ftVc[k27.ftVc$Fold < -2 & k27.ftVc$FDR < 0.01]) 
k27.hgsoc.diffbind <- k27.hgsoc.diffbind[(subjectHits(ovp2))]
k27.hgsoc.diffbind$Gene <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)$gene_id[queryHits(ovp2)]
aux <- bm[bm$entrezgene_id %in% k27.hgsoc.diffbind$Gene,]
aux <- aux[match(k27.hgsoc.diffbind$Gene, aux$entrezgene_id),]
k27.hgsoc.diffbind$Gene <- aux$hgnc_symbol
k27.hgsoc.diffbind <- as.data.frame(k27.hgsoc.diffbind, row.names = NULL)
k27.hgsoc.diffbind <- k27.hgsoc.diffbind[order(k27.hgsoc.diffbind$Fold, decreasing = F), ]
write.table(k27.hgsoc.diffbind, quote = F, row.names = F, sep = "\t", 
            file = "~/Documents/OCRA_Files/H3k27ac_DiffBind_HGSOC_table.txt")

####





##### FTSEC vs All cancer types
contrast3 <- dba.contrast(k27.db.c, group1 = k27.db.c$masks$FT, group2 = k27.db.c$masks$OvCancer)
contrast3 <- dba.analyze(contrast3)
plot(contrast3, contrast=1, colScheme = "Greys")
k27.ftVall <- dba.report(contrast3)
# k27.ftVc <- as.data.frame(k27.ftVc)
# k27.ftVc <- k27.ftVc[k27.ftVc$Fold != 0 & k27.ftVc$FDR < 0.01, ]
# for remap divide by FT and HGSOC specific
# write.table(k27.ftVc[k27.ftVc$Fold > 0, 1:3], sep= "\t", row.names = F, col.names = F, quote = F, 
#             file = "~/Documents/k27ac_FTSEC.bed")
# write.table(k27.ftVc[k27.ftVc$Fold < 0, 1:3], sep= "\t", row.names = F, col.names = F, quote = F, 
#             file = "~/Documents/k27ac_HGSOC.bed")
# 
library(ReMapEnrich)
# Use the function DowloadRemapCatalog
remapCatalog <- bedToGranges("~/Documents/remap2020_nr_macs2_hg38_v1_0.bed.gz")   
### Clean TF names
aux <- as.data.frame(remapCatalog)
aux <- separate(aux, col = "id", into = c("id", "line"), sep = ":")
remapCatalog <- makeGRangesFromDataFrame(aux, keep.extra.columns = T)

enrichment.FTvsAll <- enrichment(query = k27.ftVall[k27.ftVall$Fold > 0 & 
                                k27.ftVall$FDR <= 0.01], remapCatalog, byChrom = T)

enrichment.AllvsFT <- enrichment(query = k27.ftVall[k27.ftVall$Fold < 0 & 
                                k27.ftVall$FDR <= 0.01], remapCatalog, byChrom = T)

enrichmentBarPlot(enrichment.FTvsAll, top = 20)
enrichmentDotPlot(enrichment.FTvsAll, top = 20)

enrichmentBarPlot(enrichment.AllvsFT, top = 20)
enrichmentDotPlot(enrichment.AllvsFT, top = 20)






ggplot(k27.ftVc, aes(Fold, -log10(FDR))) + geom_point() + labs(title = "H3k27ac - FTSEC vs HGSOC")

### Plot heatmap of OCRA gene expression with specific TFs enriched in FTSEC and HGSOC based on REMAP
library(data.table)
ft.tfs <- fread("~/pacoquita_bkp/H3k27ac_FTSEC_REMAP.tab")
ft.tfs$`Transcription Factor` <- toupper(ft.tfs$`Transcription Factor`)

ft.ocra.rna <- ocra.rna[ocra.rna$Genename %in% ft.tfs$`Transcription Factor`[1:50], ]
rownames(ft.ocra.rna) <- ft.ocra.rna$Genename
ft.ocra.rna <- ft.ocra.rna[, -1]
ft.ocra.rna <- t(scale(t(ft.ocra.rna)))
pheatmap(ft.ocra.rna,  color = wes_palette(name = "Zissou1", 100, type = "continuous"),
         main= "FTSEC specific diff bound peaks")

hgsoc.tfs <- fread("~/Documents/H3k27ac_HGSOC_REMAP.tab")
hgsoc.tfs$`Transcription Factor` <- toupper(hgsoc.tfs$`Transcription Factor`)

hgsoc.ocra.rna <- ocra.rna[ocra.rna$Genename %in% hgsoc.tfs$`Transcription Factor`[1:50], ]
rownames(hgsoc.ocra.rna) <- hgsoc.ocra.rna$Genename
hgsoc.ocra.rna <- hgsoc.ocra.rna[, -1]
hgsoc.ocra.rna <- t(scale(t(hgsoc.ocra.rna)))
pheatmap(hgsoc.ocra.rna,  color = wes_palette(name = "Zissou1", 100, type = "continuous"), 
         main = "HGSOC specific diff bound peaks")



## OSEC vs HGSOC 
contrast3 <- dba.contrast(k27.db.c, group1 = k27.db.c$masks$OSEC, group2 = k27.db.c$masks$HGSOC)
contrast3 <- dba.analyze(contrast3)
plot(contrast3, contrast=1, colScheme = "Greys")
k27.oseVc <- dba.report(contrast3)

k27.oseVc <- as.data.frame(k27.oseVc)
k27.oseVc <- k27.oseVc[k27.oseVc$Fold != 0 & k27.oseVc$FDR < 0.01, ]
# for remap divide by ose and HGSOC specific
write.table(k27.oseVc[k27.oseVc$Fold > 0, 1:3], sep= "\t", row.names = F, col.names = F, quote = F, 
            file = "~/Documents/k27ac_OSEC.bed")
write.table(k27.oseVc[k27.oseVc$Fold < 0, 1:3], sep= "\t", row.names = F, col.names = F, quote = F, 
            file = "~/Documents/k27ac_OSEC_HGSOC.bed")


ggplot(k27.oseVc, aes(Fold, -log10(FDR))) + geom_point() + labs(title = "H3k27ac - OSEC vs HGSOC")

### Plot heatmap of OCRA gene expression with specific TFs enriched in OSEC and HGSOC based on REMAP
ose.tfs <- fread("~/Documents/k27ac_OSEC.tab")
ose.tfs$`Transcription Factor` <- toupper(ose.tfs$`Transcription Factor`)

ose.ocra.rna <- ocra.rna[ocra.rna$Genename %in% ose.tfs$`Transcription Factor`[1:50], ]
rownames(ose.ocra.rna) <- ose.ocra.rna$Genename
ose.ocra.rna <- ose.ocra.rna[, -1]
ose.ocra.rna <- t(scale(t(ose.ocra.rna)))
pheatmap(ose.ocra.rna,  color = wes_palette(name = "Zissou1", 100, type = "continuous"),
         main= "OSEC specific diff bound peaks")

hgsoc.tfs <- fread("~/Documents/k27ac_OSEC_HGSOC.tab")
hgsoc.tfs$`Transcription Factor` <- toupper(hgsoc.tfs$`Transcription Factor`)

hgsoc.ocra.rna <- ocra.rna[ocra.rna$Genename %in% hgsoc.tfs$`Transcription Factor`[1:50], ]
rownames(hgsoc.ocra.rna) <- hgsoc.ocra.rna$Genename
hgsoc.ocra.rna <- hgsoc.ocra.rna[, -1]
hgsoc.ocra.rna <- t(scale(t(hgsoc.ocra.rna)))
pheatmap(hgsoc.ocra.rna,  color = wes_palette(name = "Zissou1", 100, type = "continuous"), 
         main = "HGSOC specific diff bound peaks")



