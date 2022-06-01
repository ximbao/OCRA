#### OCRA H3k4me1 peaks ####
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(ReactomePA)
library(biomaRt)
library(DiffBind)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
kelly.colours <- c("gray95", "gray13", "gold2", "plum4", "darkorange1", "lightskyblue2", "firebrick", "burlywood3", "gray51", "springgreen4", "lightpink2", "deepskyblue4", "lightsalmon2", "mediumpurple4", "orange", "maroon", "yellow3", "brown4", "yellow4", "sienna4", "chocolate", "gray19")

ocra.k4me1 <- read.delim("~/Documents/ocra_h3k4me1_peaks.txt", header = F, sep = "\t", stringsAsFactors = F)
names(ocra.k4me1) <- c("SampleID", "mark", "Peaks", "bamReads")
ocra.k4me1$Replicate <- "1"
ocra.k4me1$Tissue <- c(rep("OvCancer", 3), "Normal", rep("OvCancer", 4), rep("Normal", 5), rep("OvCancer", 5))
ocra.k4me1$Factor <- c(rep("CCOC",3), "Endometrioid", "HGSOC", rep("Metastatic HGSOC",2), "LGSOC", rep("FTSEC",3),
                       rep("OSEC",2), rep("Mucinous", 3), rep('HGSOC',2))
ocra.k4me1$PeakCaller <- "narrowPeak"

k4me1.db <- dba(sampleSheet = ocra.k4me1, minOverlap = 1)
k4me1.db.c <- dba.count(k4me1.db, bParallel = T)
dba.plotPCA(k4me1.db.c, label = DBA_ID)
dba.plotHeatmap(k4me1.db.c, colScheme = "Greys")


### Making PCA without diffbind function - PRETTY PCA
aux <- dba.peakset(k4me1.db.c, bRetrieve = T) 
aux <- mcols(aux)
aux <- as.data.frame(aux)
pca1 <- prcomp(t(log2(aux+1)))
pca <- as.data.frame(pca1$x[, 1:3])
pca$anno <- ocra.k4me1$Factor
pca$shape <- ocra.k4me1$Tissue

pdf("~/Documents/OCRA_Files/OCRA_h3k4me1_PCA.pdf", 10,8)
ggplot(pca, aes(PC1, PC2, label = rownames(pca), color = shape)) + 
  geom_point(aes(size= 5,shape=anno)) + geom_text_repel(point.padding = 0.5, size = 5, segment.color = NA) + 
  theme_bw() +
  scale_color_manual(values = c("darkolivegreen4", "deepskyblue4")) +
  scale_shape_manual(values=c(15,18,22,16,3,2,4,2)) +
  coord_fixed() + 
  xlab(paste0("PC1 (", prettyNum(summary(pca1)$importance[2,1]*100, digits =2),"%)")) +
  ylab(paste0("PC2 (", prettyNum(summary(pca1)$importance[2,2]*100, digits =2),"%)")) +
  labs(color = "Histotype", shape = "Group") + guides(size = FALSE) +
  ggtitle("OCRA - H3k4me1 Peaks")
dev.off()



## Normal vs Cancer 
contrast <- dba.contrast(k4me1.db.c, group1 = k4me1.db.c$masks$Normal, group2 = k4me1.db.c$masks$OvCancer)
contrast <- dba.analyze(contrast)
plot(contrast, contrast=1)
k4me1.nVc <- dba.report(contrast)
save(k4me1.db, k4me1.db.c, k4me1.nVc, file = "~/Documents/ocra_h3k4me1_diffbind.Rda")

tagMatrix <- getTagMatrix(k4me1.nVc, windows=promoter)
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




## FTSEC vs HGSOC 
# k4me1.db.c$masks$HGSOC[c("Kuramochi")] <- TRUE
# k4me1.db.c$masks$HGSOC[c("HeyA8")] <- FALSE
contrast2 <- dba.contrast(k4me1.db.c, group1 = k4me1.db.c$masks$FT, group2 = k4me1.db.c$masks$HGSOC)
contrast2 <- dba.analyze(contrast2)
pdf("~/Documents/OCRA_Files/DiffBind_CorHeatmap_H3k4me1_FTSECvsHGSOC.pdf", 10,8)
dba.plotHeatmap(contrast2, contrast=1, colScheme = "Greys",colSideCols = c("darkolivegreen4", "deepskyblue4"), 
                main="H3k4me1 Peaks - FTSEC vs HGSOC")
dev.off()
k4me1.ftVc <- dba.report(contrast2)
k4me1.ftVc <- as.data.frame(k4me1.ftVc)
k4me1.ftVc <- k4me1.ftVc[k4me1.ftVc$Fold != 0 & k4me1.ftVc$FDR < 0.01, ]
# for remap divide by FT and HGSOC specific
write.table(k4me1.ftVc[k4me1.ftVc$Fold > 2, 1:3], sep= "\t", row.names = F, col.names = F, quote = F, 
            file = "~/Documents/OCRA_Files/k4me1_FTSEC.bed")
write.table(k4me1.ftVc[k4me1.ftVc$Fold < -2, 1:3], sep= "\t", row.names = F, col.names = F, quote = F, 
            file = "~/Documents/OCRA_Files/k4me1_HGSOC.bed")


ggplot(k27.ftVc, aes(Fold, -log10(FDR))) + geom_point() + labs(title = "H3k4me1 - FTSEC vs HGSOC")

### Plot remap enrichment ###
enrichment.diffbind.ft <- read.csv("~/Documents/OCRA_Files/Result_enrichment_H3k4me1_FTSEC.csv")
pdf("~/Documents/OCRA_Files/TFBS_H3k4me1_DiffBind_FTSEC.pdf", 8, 10)
enrichmentDotPlot(enrichment.diffbind.ft, top=20, main = "TFBS - H3k4me1 FTSEC Peaks differentially bound FTSEC vs HGSOC")
dev.off()

enrichment.diffbind.hgsoc <- read.csv("~/Documents/OCRA_Files/Result_enrichment_H3k4me1_HGSOC.csv")
pdf("~/Documents/OCRA_Files/TFBS_H3k4me1_DiffBind_HGSOC.pdf", 8, 10)
enrichmentDotPlot(enrichment.diffbind.hgsoc, top=10, main = "TFBS - H3k4me1 HGSOC Peaks differentially bound FTSEC vs HGSOC")
dev.off()

### Get genes/nearest gene on diffbound peaks ###
bm <- getBM(attributes=c('entrezgene_id', 'hgnc_symbol'),
            filters = 'chromosome_name',
            values = c(1:22,'X'),
            mart = ensembl)

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
k4me1.ftVc <- dba.report(contrast2)
## make a table with results from diffbind + gene in the diffbound peak ##
k4me1.ft.diffbind <- k4me1.ftVc[k4me1.ftVc$Fold > 2 & k4me1.ftVc$FDR < 0.01]
ovp <- findOverlaps(genes(TxDb.Hsapiens.UCSC.hg38.knownGene),
                    k4me1.ftVc[k4me1.ftVc$Fold > 2 & k4me1.ftVc$FDR < 0.01]) 
k4me1.ft.diffbind <- k4me1.ft.diffbind[(subjectHits(ovp))]
k4me1.ft.diffbind$Gene <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)$gene_id[queryHits(ovp)]
aux <- bm[bm$entrezgene_id %in% k4me1.ft.diffbind$Gene,]
aux <- aux[match(k4me1.ft.diffbind$Gene, aux$entrezgene_id),]
k4me1.ft.diffbind$Gene <- aux$hgnc_symbol
k4me1.ft.diffbind <- as.data.frame(k4me1.ft.diffbind, row.names = NULL)
k4me1.ft.diffbind <- k4me1.ft.diffbind[order(k4me1.ft.diffbind$Fold, decreasing = T), ]
write.table(k4me1.ft.diffbind, quote = F, row.names = F, sep = "\t", 
            file = "~/Documents/OCRA_Files/H3k4me1_DiffBind_FTSEC_table.txt")

k4me1.hgsoc.diffbind <- k4me1.ftVc[k4me1.ftVc$Fold < -2 & k4me1.ftVc$FDR < 0.01]
ovp2 <- findOverlaps(genes(TxDb.Hsapiens.UCSC.hg38.knownGene),
                     k4me1.ftVc[k4me1.ftVc$Fold < -2 & k4me1.ftVc$FDR < 0.01]) 
k4me1.hgsoc.diffbind <- k4me1.hgsoc.diffbind[(subjectHits(ovp2))]
k4me1.hgsoc.diffbind$Gene <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)$gene_id[queryHits(ovp2)]
aux <- bm[bm$entrezgene_id %in% k4me1.hgsoc.diffbind$Gene,]
aux <- aux[match(k4me1.hgsoc.diffbind$Gene, aux$entrezgene_id),]
k4me1.hgsoc.diffbind$Gene <- aux$hgnc_symbol
k4me1.hgsoc.diffbind <- as.data.frame(k4me1.hgsoc.diffbind, row.names = NULL)
k4me1.hgsoc.diffbind <- k4me1.hgsoc.diffbind[order(k4me1.hgsoc.diffbind$Fold, decreasing = F), ]
write.table(k4me1.hgsoc.diffbind, quote = F, row.names = F, sep = "\t", 
            file = "~/Documents/OCRA_Files/H3k4me1_DiffBind_HGSOC_table.txt")





### Plot heatmap of OCRA gene expression with specific TFs enriched in FTSEC and HGSOC based on REMAP
ft.tfs <- fread("~/Documents/k4me1_FTSEC.tab")
ft.tfs$`Transcription Factor` <- toupper(ft.tfs$`Transcription Factor`)

ft.ocra.rna <- ocra.rna[ocra.rna$Genename %in% ft.tfs$`Transcription Factor`, ]
rownames(ft.ocra.rna) <- ft.ocra.rna$Genename
ft.ocra.rna <- ft.ocra.rna[, -1]
ft.ocra.rna <- t(scale(t(ft.ocra.rna)))
pheatmap(ft.ocra.rna,  color = wes_palette(name = "Zissou1", 100, type = "continuous"),
         main= "FTSEC specific diff bound peaks")

hgsoc.tfs <- fread("~/Documents/k4me1_FTSEC_HGSOC.tab")
hgsoc.tfs$`Transcription Factor` <- toupper(hgsoc.tfs$`Transcription Factor`)

hgsoc.ocra.rna <- ocra.rna[ocra.rna$Genename %in% hgsoc.tfs$`Transcription Factor`, ]
rownames(hgsoc.ocra.rna) <- hgsoc.ocra.rna$Genename
hgsoc.ocra.rna <- hgsoc.ocra.rna[, -1]
hgsoc.ocra.rna <- t(scale(t(hgsoc.ocra.rna)))
pheatmap(hgsoc.ocra.rna,  color = wes_palette(name = "Zissou1", 100, type = "continuous"), 
         main = "HGSOC specific diff bound peaks")






## OSEC vs HGSOC 
contrast3 <- dba.contrast(k4me1.db.c, group1 = k4me1.db.c$masks$OSEC, group2 = k4me1.db.c$masks$HGSOC)
contrast3 <- dba.analyze(contrast3)
plot(contrast3, contrast=1, colScheme = "Greys")
k4me1.oseVc <- dba.report(contrast3)

k4me1.oseVc <- as.data.frame(k4me1.oseVc)
k4me1.oseVc <- k4me1.oseVc[k4me1.oseVc$Fold != 0 & k4me1.oseVc$FDR < 0.01, ]
# for remap divide by ose and HGSOC specific
write.table(k4me1.oseVc[k4me1.oseVc$Fold > 0, 1:3], sep= "\t", row.names = F, col.names = F, quote = F, 
            file = "~/Documents/k4me1ac_OSEC.bed")
write.table(k4me1.oseVc[k4me1.oseVc$Fold < 0, 1:3], sep= "\t", row.names = F, col.names = F, quote = F, 
            file = "~/Documents/k4me1ac_OSEC_HGSOC.bed")


ggplot(k4me1.oseVc, aes(Fold, -log10(FDR))) + geom_point() + labs(title = "H3k4me1 - OSEC vs HGSOC")

### Plot heatmap of OCRA gene expression with specific TFs enriched in FTSEC and HGSOC based on REMAP
ose.tfs <- fread("~/Documents/k4me1_OSEC.tab")
ose.tfs$`Transcription Factor` <- toupper(ose.tfs$`Transcription Factor`)

ose.ocra.rna <- ocra.rna[ocra.rna$Genename %in% ose.tfs$`Transcription Factor`, ]
rownames(ose.ocra.rna) <- ose.ocra.rna$Genename
ose.ocra.rna <- ose.ocra.rna[, -1]
ose.ocra.rna <- t(scale(t(ose.ocra.rna)))
pheatmap(na.omit(ose.ocra.rna),  color = wes_palette(name = "Zissou1", 100, type = "continuous"),
         main= "OSEC specific diff bound peaks")

hgsoc.tfs <- fread("~/Documents/k4me1_OSEC_HGSOC.tab")
hgsoc.tfs$`Transcription Factor` <- toupper(hgsoc.tfs$`Transcription Factor`)

hgsoc.ocra.rna <- ocra.rna[ocra.rna$Genename %in% hgsoc.tfs$`Transcription Factor`, ]
rownames(hgsoc.ocra.rna) <- hgsoc.ocra.rna$Genename
hgsoc.ocra.rna <- hgsoc.ocra.rna[, -1]
hgsoc.ocra.rna <- t(scale(t(hgsoc.ocra.rna)))
pheatmap(hgsoc.ocra.rna,  color = wes_palette(name = "Zissou1", 100, type = "continuous"), 
         main = "HGSOC specific diff bound peaks")





