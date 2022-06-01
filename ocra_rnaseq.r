##### OCRA Rnaseq data ####### - using cpm matrix / downloaded from GENAVi
library(edgeR)
library(biomaRt)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
hmart=useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
bm <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'),
            filters = 'chromosome_name',
            values = c(1:22,'X'),
            mart = hmart)


#ocra.rna <- read.table("~/Documents/ocra_rna_cpm.csv", sep = ",", stringsAsFactors = F, header = T)
#ocra.rna <- ocra.rna[, -c(1,3:8, 26:30)]
#genes <- as.character(ocra.rna$Genename)
#bm <- bm[bm$hgnc_symbol %in% genes,]

aux <- read.delim("~/OCRA/RNASeq/Cell_Lines_RNA-seq_final_table.txt", sep = "\t", header = T, skip = 1)
names(aux)[7:26] <- substr(names(aux)[7:26], 48, 108)
names(aux)[7:26] <- gsub("_.*$", "", names(aux)[7:26])
ocra.rna <- aux[, -c(7,20:21)]
#ocra.rna <- DGEList(counts = aux[, 7:23], genes = aux$Geneid)
#keep <- filterByExpr(ocra.rna, min.count = 1)
#ocra.rna <- ocra.rna[keep,]
#aux <- as.data.frame(ocra.rna$counts)
#aux$Length <- as.character(ocra.rna$genes$genes)
#aux$Length <- substr(aux$Length, 1,15)
#ocra.rna <- aux[aux$Length %in% bm$ensembl_gene_id, ]

## Gotta add FT282 in this data-frame... 
ocra.raw <- read.csv("~/OCRA/RNASeq/OCRA_RNAseq_rawCounts.csv")
ocra.raw$FT282 <- round((ocra.raw$FT282.A + ocra.raw$FT282.B) / 2)
ocra.rna <- merge(ocra.rna, ocra.raw[, c("Geneid" , "FT282")], by = "Geneid", all.x = F) %>% na.omit()

################# PCA ####################
## ocra table for color codes 
ocra.colors <- data.frame("CellLine" = colnames(ocra.rna)[7:24], "Histotype" = NA, "Color" = NA, "Group" = NA)
ocra.colors$Histotype[ocra.colors$CellLine %in% c("CaOV3", "kuramochi", "UWB1")] <- "HGSOC"
#ocra.colors$Histotype[ocra.colors$CellLine %in% c("HEY", "OAW42")] <- "Metastatic HGSOC"
ocra.colors$Histotype[ocra.colors$CellLine %in% c("HEY", "OAW42")] <- "HGSOC"
ocra.colors$Histotype[ocra.colors$CellLine %in% c("FT33", "FT246")] <- "FTSEC"
#ocra.colors$Histotype[ocra.colors$CellLine %in% c("FT282")] <- "FTSEC+tp53"
ocra.colors$Histotype[ocra.colors$CellLine %in% c("FT282")] <- "FTSEC"
ocra.colors$Histotype[ocra.colors$CellLine %in% c("IOSE4", "IOSE11")] <- "OSEC"
ocra.colors$Histotype[ocra.colors$CellLine %in% c("EFO27", "GFTR230", "MCAS")] <- "Mucinous"
ocra.colors$Histotype[ocra.colors$CellLine %in% c("JHOC5", "RMGII", "ES2")] <- "CCOC"
ocra.colors$Histotype[ocra.colors$CellLine %in% c("VOA1056")] <- "LGSOC"
ocra.colors$Histotype[ocra.colors$CellLine %in% c("EEC16")] <- "Endometrioid"
#assign colors to the groups for all the plots to come!!
ocra.colors$Color[ocra.colors$Histotype %in% "HGSOC"] <- "#8EE7E5"
#ocra.colors$Color[ocra.colors$Histotype %in% "Metastatic HGSOC"] <- "#4D95B9"
ocra.colors$Color[ocra.colors$Histotype %in% "FTSEC"] <- "#97E17A"
#ocra.colors$Color[ocra.colors$Histotype %in% "FTSEC+tp53"] <- "#3D8C23"
ocra.colors$Color[ocra.colors$Histotype %in% "OSEC"] <- "#FED966"
ocra.colors$Color[ocra.colors$Histotype %in% "Mucinous"] <- "#9851CB"
ocra.colors$Color[ocra.colors$Histotype %in% "CCOC"] <- "#D30347"
ocra.colors$Color[ocra.colors$Histotype %in% "LGSOC"] <- "#7F7F7F"
ocra.colors$Color[ocra.colors$Histotype %in% "Endometrioid"] <- "#FD9933"
ocra.colors$Group[ocra.colors$Histotype %in% c("HGSOC", "Mucinous", "CCOC", "LGSOC", "Metastatic HGSOC")] <- "OvCancer"
ocra.colors$Group[ocra.colors$Histotype %in% c("FTSEC", "OSEC", "Endometrioid", "FTSEC+tp53")] <- "Normal"
#write.table(ocra.colors, row.names = F, sep = "\t", file = "~/OCRA/OCRA_Colors_metadata.txt")

pca <- prcomp(t(ocra.rna[, 7:24]))
aux <- as.data.frame(pca$x)

shape <- (ocra.colors$Group)
levels(shape) <- c('1', '2')

pdf("~/Documents/OCRA_Files/RNASEQ_PCA_updated3.pdf", 10, 10)
ggplot(aux, aes(x=PC1, y=PC2, color = ocra.colors$Histotype[match(rownames(aux), ocra.colors$CellLine)])) +
  geom_point(aes(shape=shape, size = 5)) + 
  scale_color_manual(values = c("#D4126A", "#FD9933", "#97E17A", "#34CCCC", "#7F7F7F", "#6D69B3",
                                "#FF6699")) +
  xlab(paste0("PC1: ",prettyNum(summary(pca)$importance[2,1]*100,
                                digits = 2, decimal.mark = "."),"% variance")) +
  ylab(paste0("PC2: ",prettyNum(summary(pca)$importance[2,2]*100,
                                digits = 2, decimal.mark = "."),"% variance")) +
  coord_fixed() + ggtitle("PCA - OCRA Cell lines") +
  geom_text_repel(data=aux[, 1:2], label = rownames(aux)) + labs(color = "Histotype", shape = "Group") +
  guides(size = FALSE) + 
  theme_bw() + theme( panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
dev.off()

# rownames(ocra.rna) <- make.names(ocra.rna$Genes, unique = T)
# m <- match(rownames(ocra.rna), bm$ensembl_gene_id)
# rownames(ocra.rna) <- bm$hgnc_symbol[m]
# ocra.rna <- ocra.rna[, -1]
# ocra.rna <- DGEList(ocra.rna)
# y <- calcNormFactors(ocra.rna)
# y <- estimateDisp(y, NULL)
# y <- glmFit(y, NULL)
# head(y$fitted.values)
# pca <- prcomp(t(y$fitted.values))

gtf <- readGFF("~/reference_files/gencode.v37.annotation.gtf")
gtf <- gtf[gtf$type %in% "gene",]
rownames(ocra.rna) <- ocra.rna$Geneid
#### DEG
library(DESeq2)
coldata <- data.frame(row.names = colnames(ocra.rna[7:24]), condition = ocra.colors$Group)
dds <- DESeqDataSetFromMatrix(countData = ocra.rna[, 7:24],
                              colData = coldata,
                              design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
res <- na.omit(res)
res$STATUS <- "Not Sig"
res[res$log2FoldChange < -1 & res$padj < 0.05, ]$STATUS <- "Down"
res[res$log2FoldChange > 1 & res$padj < 0.05, ]$STATUS <- "Up"

aux <- as.data.frame(res)

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
genes.mart <-
  getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
        values = na.omit(ocra.rna$Geneid), mart = mart)
aux$Symbol <- genes.mart$external_gene_name[match(substr(rownames(aux), 1,15), genes.mart$ensembl_gene_id)]

# label top DEGs on volcano
aux$label <- NA
up <- aux[aux$STATUS == "Up", ]
up <- arrange(up, -log2FoldChange, padj)
#up <- arrange(up,  padj)

w <- which(up$Symbol != "")
up$label[w[1:10]] <- up$Symbol[w[1:10]]

dn <- aux[aux$STATUS == "Down", ]
dn <- arrange(dn, log2FoldChange, padj)
#dn <- arrange(dn, padj)

dn$label[1:10] <- dn$Symbol[1:10]

# transfer label to the plotting obj
aux$label[which(aux$Symbol %in% up$label[1:10])] <- up$label[1:10]
aux$label[which(aux$Symbol %in% dn$label[1:10])] <- dn$label[1:10]


pdf("~/Documents/OCRA_Files/Volcano_NormalsvsOvCancer_updated.pdf", 15, 10)
ggplot(aux, aes(log2FoldChange, -log10(padj), col = STATUS, label = label)) + geom_point(size=3) + 
  geom_label_repel(aes(label = aux$label), vjust = "inward", hjust = "inward") + xlab("log2 Fold Change") +
  ylab("-log10 Adjusted Pvalue") + ggtitle("Normals vs OvCancer") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
   scale_color_manual(values = c("darkolivegreen3",  "grey", "deepskyblue3"))
dev.off()

 
# vst <- vst(dds)
# hm <- assay(vst)
# hm <- hm[rownames(res)[!res$STATUS %in% "Not Sig"] , ]

hm <- ocra.rna[rownames(ocra.rna) %in% rownames(res)[res$STATUS != "Not Sig"],]


table.deg.normalXcancer <- res[res$STATUS != "Not Sig", ]
table.deg.normalXcancer$GeneSymbol <- gtf$gene_name[match(rownames(table.deg.normalXcancer), gtf$gene_id)]
write.table(table.deg.ftsecXhgsoc, quote = F, sep = "\t", file = "~/Documents/OCRA_Files/Table_DEG_Normal_vs_Cancer.txt")


# col_anno <- data.frame(names = colnames(hm), sample = c("HGSOC", "Endometriosis", "Mucinous", "CCOC", rep("FTSEC",2),
#                                                             "Mucinous","Metastatic HGSOC", rep("OSEC",2), "CCOC", "Mucinous",
#                                                             "Metastatic HGSOC", "CCOC", "HGSOC", "LGSOC", "HGSOC"),
#                        color = ggplot_build(p)$data[[2]]$fill)
# 
# color_anno <- levels(as.factor(ggplot_build(p)$data[[2]]$fill))
# names(color_anno) <- levels(as.factor(col_anno$sample))
# library(ComplexHeatmap)
# ha = HeatmapAnnotation(
#   condition = coldata$condition, 
#   histotype = col_anno$sample,
#   col = list(condition = c("OvCancer" = "#8362A0", "Normal" = "#76DDFD"),
#              histotype = c("HGSOC" = "#00BE67", "Endometriosis" = "#CD9600", "Mucinous" = "#C77CFF", "CCOC" = "#F8766D", "FTSEC" = "#7CAE00",
#                            "Metastatic HGSOC" = "#00A9FF", "OSEC" = "#FF61CC", "LGSOC" = "#00BFC4")
#   )
# )

#pheatmap(hm, show_rownames = F, annotation_col = col_anno$sample, annotation_colors = col_anno$color)
# pdf("~/Documents/OCRA_RNASeq_HM.pdf", 10,8)
# Heatmap(hm, show_row_names = F, top_annotation = ha, name = "Gene Counts (vst)", col = viridis(50))
# dev.off()

anno <- data.frame("Group" = ocra.colors$Group, "Histotype" = ocra.colors$Histotype, row.names = ocra.colors$CellLine)
Var1        <- unique(ocra.colors$Color)
names(Var1) <- unique(ocra.colors$Histotype)
Var2        <- factor(ocra.colors$Group, labels = c("#658D2E", "#00698E"))
Var2 <- unique(as.character(Var2))
names(Var2) <- unique(ocra.colors$Group)
anno_colors <- list(Group = Var2, Histotype = Var1)
#change color to match pca
anno_colors$Histotype[6] <- "#FF6699"
anno_colors$Histotype[3] <- "#6D69B3"

pdf("~/Documents/OCRA_Files/OCRA_OvCancerNormals_Heatmap_updated.pdf", 10,7)
p1 <- pheatmap(t(scale(t(hm[, 7:24]))), show_rownames = F, col = viridis(50), annotation_col = anno, 
        annotation_colors = anno_colors ,
        main = "RNASeq - OvCancer vs Normal, unsupervised clustering", border_color = NA)
p1
dev.off()

#### DEG                          FTSEC vs HGSOC                                                      ####
library(DESeq2)
coldata.sub1 <- subset(ocra.colors, ocra.colors$Histotype %in% c("HGSOC", "FTSEC", "Metastatic HGSOC", "FTSEC+tp53"))
coldata.sub1 <- data.frame(row.names = coldata.sub1$CellLine, condition = c("HGSOC", "FTSEC", "FTSEC", 
                                                                            rep("HGSOC", 4), "FTSEC"))

# Just FT244, FT33 vs Kura & UWB
aux <- data.frame(row.names = c("FT246", "FT33", "UWB1", "kuramochi"), 
                  condition = c(rep("FTSEC", 2), rep("HGSOC", 2)))
  
dds2 <- DESeqDataSetFromMatrix(countData = ocra.rna[, colnames(ocra.rna) %in% rownames(coldata.sub1)],
                              colData = coldata.sub1,
                              design = ~ condition)

dds2 <- DESeqDataSetFromMatrix(countData = ocra.rna[, colnames(ocra.rna) %in% rownames(aux)],
                               colData = aux,
                               design = ~ condition)

dds2 <- DESeq(dds2)
res2 <- results(dds2)
res2 <- na.omit(res2)
res2$STATUS <- "Not Sig"
res2[res2$log2FoldChange < -1 & res2$padj < 0.05, ]$STATUS <- "Down"
res2[res2$log2FoldChange > 1 & res2$padj < 0.05, ]$STATUS <- "Up"



aux <- as.data.frame(res2)
write.table(aux, quote = F, sep = "\t", file = "~/Documents/OCRA_Files/DEG_HGSOC_FTSEC_no282.txt")

aux$Symbol <- genes.mart$external_gene_name[match(substr(rownames(aux), 1,15), genes.mart$ensembl_gene_id)]
# label top DEGs on volcano
aux$label <- NA
up <- aux[aux$STATUS == "Up", ]
up <- arrange(up, -log2FoldChange, padj)
#up <- arrange(up,  padj)

up$label[1:10] <- up$Symbol[1:10]

dn <- aux[aux$STATUS == "Down", ]
dn <- arrange(dn, log2FoldChange, padj)
#dn <- arrange(dn, padj)

dn$label[1:10] <- dn$Symbol[1:10]

# transfer label to the plotting obj
aux$label[which(aux$Symbol %in% up$label[1:10])] <- up$label[1:10]
aux$label[which(aux$Symbol %in% dn$label[1:10])] <- dn$label[1:10]


pdf("~/Documents/OCRA_Files/Volcano_FTSECvsHGSOC.pdf", 15, 10)
ggplot(aux, aes(log2FoldChange, -log10(padj), col = STATUS, label = label)) + geom_point() + 
  geom_label_repel(aes(label = aux$label), vjust = "inward", hjust = "inward") + xlab("log2 Fold Change") +
  ylab("-log10 Adjusted Pvalue") + ggtitle("FTSEC vs HGSOC") +theme_bw() +
  scale_color_manual(values = c("darkolivegreen4",  "grey", "deepskyblue4"))
dev.off()

# vst <- vst(dds)
# hm <- assay(vst)
# hm <- hm[rownames(res)[!res$STATUS %in% "Not Sig"] , ]
table.deg.ftsecXhgsoc <- res2[res2$STATUS != "Not Sig", ]
table.deg.ftsecXhgsoc$GeneSymbol <- gtf$gene_name[match(rownames(table.deg.ftsecXhgsoc), gtf$gene_id)]
write.table(table.deg.ftsecXhgsoc, quote = F, sep = "\t", file = "~/Documents/OCRA_Files/Table_DEG_FTSEC_vs_HGSOC.txt")

####                                                                                                  ####

#### DEG                          IOSE vs HGSOC                                                      ####
library(DESeq2)
coldata.sub2 <- subset(ocra.colors, ocra.colors$Histotype %in% c("HGSOC", "OSEC", "Metastatic HGSOC"))
coldata.sub2 <- data.frame(row.names = coldata.sub2$CellLine, condition = c("HGSOC", "HGSOC", "OSEC", 
                                                                            "OSEC", rep("HGSOC",3)))


dds3 <- DESeqDataSetFromMatrix(countData = ocra.rna[, colnames(ocra.rna) %in% rownames(coldata.sub2)],
                               colData = coldata.sub2,
                               design = ~ condition)
dds3 <- DESeq(dds3)
res3 <- results(dds3)
res3 <- na.omit(res3)
res3$STATUS <- "Not Sig"
res3[res3$log2FoldChange < -1 & res3$padj < 0.05, ]$STATUS <- "Down"
res3[res3$log2FoldChange > 1 & res3$padj < 0.05, ]$STATUS <- "Up"


aux <- as.data.frame(res3)

aux$Symbol <- genes.mart$external_gene_name[match(substr(rownames(aux), 1,15), genes.mart$ensembl_gene_id)]
aux$Symbol <- make.names(aux$Symbol, unique = T)

# label top DEGs on volcano
aux$label <- NA
up <- aux[aux$STATUS == "Up", ]
up <- arrange(up, -log2FoldChange, padj)
up <- arrange(up,  padj)

up$label[1:10] <- up$Symbol[1:10]

dn <- aux[aux$STATUS == "Down", ]
dn <- arrange(dn, log2FoldChange, padj)
dn <- arrange(dn, padj)

dn$label[1:10] <- dn$Symbol[1:10]

# transfer label to the plotting obj
aux$label[which(aux$Symbol %in% up$label[1:10])] <- up$label[1:10]
aux$label[which(aux$Symbol %in% dn$label[1:10])] <- dn$label[1:10]


pdf("~/Documents/OCRA_Files/Volcano_IOSEvsHGSOC.pdf", 15, 10)
ggplot(aux, aes(log2FoldChange, -log10(padj), col = STATUS, label = label)) + geom_point() + 
  geom_label_repel(aes(label = aux$label), vjust = "inward", hjust = "inward") + xlab("log2 Fold Change") +
  ylab("-log10 Adjusted Pvalue") + ggtitle("IOSE vs HGSOC") +theme_bw() +
  scale_color_manual(values = c("darkolivegreen4",  "grey", "deepskyblue4"))
dev.off()

# vst <- vst(dds)
# hm <- assay(vst)
# hm <- hm[rownames(res)[!res$STATUS %in% "Not Sig"] , ]
table.deg.osecXhgsoc <- res3[res3$STATUS != "Not Sig", ]
table.deg.osecXhgsoc$GeneSymbol <- gtf$gene_name[match(rownames(table.deg.osecXhgsoc), gtf$gene_id)]
write.table(table.deg.osecXhgsoc, quote = F, sep = "\t", file = "~/Documents/OCRA_Files/Table_DEG_OSEC_vs_HGSOC.txt")


####                                                                                                  ####

#### DEG                          IOSE vs FTSEC                                                      ####
coldata.sub3 <- subset(ocra.colors, ocra.colors$Histotype %in% c("OSEC", "FTSEC"))
coldata.sub3 <- data.frame(row.names = coldata.sub3$CellLine, condition = c(rep("FTSEC",2), rep("OSEC",2)))

dds4 <- DESeqDataSetFromMatrix(countData = ocra.rna[, colnames(ocra.rna) %in% rownames(coldata.sub3)],
                               colData = coldata.sub3,
                               design = ~ condition)
dds4 <- DESeq(dds4)
res4 <- results(dds4)
res4 <- na.omit(res4)
res4$STATUS <- "Not Sig"
res4[res4$log2FoldChange < -1 & res4$padj < 0.05, ]$STATUS <- "Down"
res4[res4$log2FoldChange > 1 & res4$padj < 0.05, ]$STATUS <- "Up"

#save(ocra.rna, res, res2, res3, res4, file = "~/OCRA/DEG_Objects_OCRA.Rda")


aux <- as.data.frame(res4)

aux$Symbol <- genes.mart$external_gene_name[match(substr(rownames(aux), 1,15), genes.mart$ensembl_gene_id)]
# label top DEGs on volcano
aux$label <- NA
up <- aux[aux$STATUS == "Up", ]
up <- arrange(up, -log2FoldChange, padj)
up <- arrange(up,  padj)

up$label[1:10] <- up$Symbol[1:10]

dn <- aux[aux$STATUS == "Down", ]
dn <- arrange(dn, log2FoldChange, padj)
dn <- arrange(dn, padj)

dn$label[1:10] <- dn$Symbol[1:10]

# transfer label to the plotting obj
aux$label[which(aux$Symbol %in% up$label[1:10])] <- up$label[1:10]
aux$label[which(aux$Symbol %in% dn$label[1:10])] <- dn$label[1:10]


pdf("~/Documents/OCRA_Files/Volcano_FTSECvsIOSE.pdf", 15, 10)
ggplot(aux, aes(log2FoldChange, -log10(padj), col = STATUS, label = label)) + geom_point() + 
  geom_label_repel(aes(label = aux$label), vjust = "inward", hjust = "inward") + xlab("log2 Fold Change") +
  ylab("-log10 Adjusted Pvalue") + ggtitle("FTSEC vs IOSE") +theme_bw() +
  scale_color_manual(values = c("darkolivegreen4",  "grey", "deepskyblue4"))
dev.off()



# vst <- vst(dds)
# hm <- assay(vst)
# hm <- hm[rownames(res)[!res$STATUS %in% "Not Sig"] , ]
table.deg.osecXftsec <- res4[res4$STATUS != "Not Sig", ]
table.deg.osecXftsec$GeneSymbol <- gtf$gene_name[match(rownames(table.deg.osecXftsec), gtf$gene_id)]
write.table(table.deg.osecXftsec, quote = F, sep = "\t", file = "~/Documents/OCRA_Files/Table_DEG_OSEC_vs_FTSEC.txt")

## Compare results with OSEC vs FTSEC from the origins paper (using tissues)
library(readxl)
origins <- read_excel("~/Documents/OCRA_Files/1-s2.0-S221112471931455X-mmc3.xlsx", skip = 2) 

intersect(table.deg.osecXftsec$GeneSymbol, origins$`Gene ID`)

deseq2.fc=table.deg.osecXftsec$log2FoldChange
#names(deseq2.fc)=rownames(table.deg.osecXftsec)
exp.fc=deseq2.fc
out.suffix="deseq2"
# Next, GAGE for pathway analysis, and Pathview for visualization. We only visualize up-regulated pathways
# here, down-regulated pathways can be done the same way (see also the above native workflow). Notice that
# this step (the same code) is identical for DESeq, edgeR, Limma and Cufflinks workflows, hence skipped in
# other workflows below.

library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
genes.mart <- 
  getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), 
        values = na.omit(table.deg.osecXftsec$GeneSymbol), mart = mart)

rownames(table.deg.osecXftsec) <- substr(rownames(table.deg.osecXftsec), 1, 15)
table.deg.osecXftsec$EntrezID <- genes.mart$entrezgene_id[match(rownames(table.deg.osecXftsec), 
                                                                genes.mart$ensembl_gene_id)]

names(deseq2.fc) <- table.deg.osecXftsec$EntrezID

# require(gage)
data(kegg.gs)
fc.kegg.p <- gage(exp.fc, gsets = kegg.gs, ref = NULL, samp = NULL)
sel <- fc.kegg.p$greater[, "q.val"] < 0.1 & !is.na(fc.kegg.p$greater[, "q.val"])
path.ids <- rownames(fc.kegg.p$greater)[sel]
sel.l <- fc.kegg.p$less[, "q.val"] < 0.1 & !is.na(fc.kegg.p$less[,"q.val"])
path.ids.l <- rownames(fc.kegg.p$less)[sel.l]
path.ids2 <- substr(c(path.ids, path.ids.l), 1, 8)
require(pathview)
#view first 3 pathways as demo
pv.out.list <- sapply(path.ids2[1:3], function(pid) pathview(gene.data = exp.fc, pathway.id = pid, species = "hsa", out.suffix=out.suffix))



##

## prepare data for pathway 
fc.vector <- br_oe$FC
names(fc.vector) <- br_oe$entrez
fc.vector <- na.omit(fc.vector)



keggres = gage(fc.vector, gsets=kegg.sets.hs, same.dir=TRUE)


pathways = data.frame(id=rownames(keggres$greater), keggres$less)
head(pathways)

keggrespathways <- rownames(keggres$greater)[1:10]


######## 
library(org.Hs.eg.db)
ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=substr(rownames(res),1,15), 
                                    columns="SYMBOL",
                                    keytype="ENSEMBL")
ens2symbol <- as_tibble(ens2symbol)
ens2symbol

# OvCancer vs Normal
res2$row <- substr(rownames(res2),1,15)
res2 <- as_tibble(res2)
aux <- inner_join(res2, ens2symbol, by=c("row"="ENSEMBL"))

aux2 <- aux %>% 
  dplyr::select(SYMBOL, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(SYMBOL) %>% 
  summarize(stat=mean(stat))
aux2

library(fgsea)
ranks <- deframe(aux2)
head(ranks, 20)

## Choose the pathway from the downloaded sets on ~/Documents/pathways_dir/
pathways.hallmark <- gmtPathways("~/Documents/pathways_dir/h.all.v7.2.symbols.gmt")
go.bp <- gmtPathways("~/Documents/pathways_dir/c5.go.bp.v7.2.symbols.gmt")
kegg.path <- gmtPathways("~/Documents/pathways_dir/c2.cp.kegg.v7.2.symbols.gmt")

fgseaRes <- fgsea(pathways=kegg.path, stats=ranks)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

pdf("~/Documents/OCRA_Files/Pathway_kegg_FTSECvsHGSOC.pdf", 15, 26)
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Kegg pathway NES - FTSEC vs HGSOC DEGs") + 
  theme_minimal() + scale_fill_viridis_d()
dev.off()

