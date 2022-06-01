#################### PCA Of all epigenetic marks for OCRA #######################
load("~/pacoquita_bkp/env_nov2019.rda")
load("~/Documents/OCRA_Files/OCRA_ALl_Marks_DiffBind.Rda")

# k27.pk <- dba.peakset(k27.db.c, bRetrieve = T) 
# k4me1.pk <- dba.peakset(k4me1.db.c, bRetrieve = T)
# k4me3.pk <- dba.peakset(k4me3.db.c, bRetrieve = T)
# ctcf.pk <- dba.peakset(ctcf.db.c, bRetrieve = T)
# 

names(manifest)[2] <- "mark"
# all.db1 <- rbind(ocra.27ac, ocra.k4me1, ocra.k4me3, manifest[1:18,])
# all.db1$Factor <- all.db1$mark
#                                                                             ALREADY RAN THIS ON PACOQUITA 
# all.db <- dba(sampleSheet = all.db1, minOverlap = 1)
# all.db.c <- dba.count(all.db, bParallel = T)
all.db.c$samples$SampleID[1:18] <- paste0(all.db.c$samples$SampleID[1:18], "_H3k27ac")
all.db.c$samples$SampleID[19:36] <- paste0(all.db.c$samples$SampleID[19:36], "_H3k4me1")
all.db.c$samples$SampleID[37:54] <- paste0(all.db.c$samples$SampleID[37:54], "_H3k4me3")
all.db.c$samples$SampleID[55:72] <- paste0(all.db.c$samples$SampleID[55:72], "_CTCF")


dba.plotPCA(all.db.c, label = DBA_ID, attributes = DBA_FACTOR, vColors = kelly.colours[2:5], dotSize = 1, 
            labelSize = 0.8, score = DBA_SCORE_TMM_READS_FULL_CPM)

color.vec <- c("#4B3A8F", "#24504F", "#CFBE62", "#73889A")
names(color.vec) <- c("CTCF", "H3k27ac", "H3k4me1", "H3k4me3")
color.vec2 <- c("darkolivegreen4", "deepskyblue4")
names(color.vec2) <- c("Normal", "OvCancer")
hmcol <- list("Group" = color.vec2, "Factor" = color.vec) 

pdf("~/Documents/OCRA_Files/All_Marks_CorrHeatmap.pdf", 15,10)
dba.plotHeatmap(all.db.c, attributes = DBA_ID, colScheme = "Greys", score = DBA_SCORE_RPKM,
                 colSideCols = hmcol, main = "Correlation Heatmap - All Marks, unsupervised clustering")
dev.off()
# mark.order.hm <-rownames(aux)
# plot.hm.marks <- aux
# 
# pheatmap(plot.hm.marks, cluster_rows = F, cluster_cols = F)



contrast <- dba.contrast(all.db.c, group1 = all.db.c$masks$Normal, group2 = all.db.c$masks$OvCancer)
contrast <- dba.analyze(contrast)
plot(contrast, contrast=1)
all.nVc <- dba.report(contrast)


#######
all.peakset <- dba.peakset(all.db.c, bRetrieve = T)
aux <- mcols(all.peakset)
aux <- as.data.frame(aux)
colnames(aux)[1:18] <- paste0(colnames(aux)[1:18], "_H3k27ac")
colnames(aux)[19:36] <- paste0(colnames(aux)[19:36], "_H3k4me1")
colnames(aux)[37:54] <- paste0(colnames(aux)[37:54], "_H3k4me3")
colnames(aux)[55:72] <- paste0(colnames(aux)[55:72], "_CTCF")
colnames(aux) <- gsub("\\.[1-3]", "", colnames(aux))

pca1 <- prcomp(t(log2(aux+.1)))
pca <- as.data.frame(pca1$x[, 1:3])
pca$anno <- all.db1$Factor
pca$shape <- all.db1$Tissue

p1 <- ggplot(pca, aes(PC1, PC2, label = rownames(pca), color = anno)) + 
  geom_point(aes(shape=shape, size = 3)) + 
  #geom_text_repel(point.padding = 0.5, size = 3.5, segment.color = NA) + 
  theme_bw() +
 # scale_color_manual(values = kelly.colours[5:dim(pca)[1]]) +
  scale_color_manual(values = c("#4B3A8F", "#24504F", "#CFBE62", "#73889A")) +
  coord_fixed() + 
  xlab(paste0("PC1 (", prettyNum(summary(pca1)$importance[2,1]*100, digits =2),"%)")) +
  ylab(paste0("PC2 (", prettyNum(summary(pca1)$importance[2,2]*100, digits =2),"%)")) +
  labs(color = "Histone modification", shape = "Group") + guides(size = FALSE) +
  ggtitle("OCRA - All Marks")
pdf("~/Documents/OCRA_Files/PCA_All_Marks_NOLABEL.pdf", 10,10)
p1
dev.off()

p2 <- ggplot(pca, aes(PC1, PC3, label = rownames(pca), color = anno)) + 
  geom_point(aes(shape=pca$shape, size = 3)) + geom_text_repel(point.padding = 0.5, size = 2.5, segment.color = NA) + 
  theme_bw() +
  scale_color_manual(values = kelly.colours[5:dim(pca)[1]]) +
  coord_fixed() + 
  xlab(paste0("PC1 (", prettyNum(summary(pca1)$importance[2,1]*100, digits =2),"%)")) +
  ylab(paste0("PC3 (", prettyNum(summary(pca1)$importance[2,3]*100, digits =2),"%)")) +
  labs(color = "Histotype", shape = "Group") + guides(size = FALSE) +
  ggtitle("OCRA - All Marks")
ggarrange(p1,p2)



#### No CTCF
aux2 <- as.data.frame(aux[, 1:54])
pca2 <- prcomp(t(log2(aux2+.1)))
pca <- as.data.frame(pca2$x[, 1:3])
pca$anno <- all.db1$Factor[1:54]
pca$shape <- all.db1$Tissue[1:54]

p3 <- ggplot(pca, aes(PC1, PC2, label = rownames(pca), color = anno)) + 
  geom_point(aes(shape=pca$shape, size = 3)) + geom_text_repel(point.padding = 0.5, size = 2.5, segment.color = NA) + 
  theme_bw() +
  scale_color_manual(values = kelly.colours[5:dim(pca)[1]]) +
  coord_fixed() + 
  xlab(paste0("PC1 (", prettyNum(summary(pca2)$importance[2,1]*100, digits =2),"%)")) +
  ylab(paste0("PC2 (", prettyNum(summary(pca2)$importance[2,2]*100, digits =2),"%)")) +
  labs(color = "Histotype", shape = "Group") + guides(size = FALSE) +
  ggtitle("OCRA - H3k27ac, H3k4me1, H3k4me3")

p4 <- ggplot(pca, aes(PC1, PC3, label = rownames(pca), color = anno)) + 
  geom_point(aes(shape=pca$shape, size = 3)) + geom_text_repel(point.padding = 0.5, size = 2.5, segment.color = NA) + 
  theme_bw() +
  scale_color_manual(values = kelly.colours[5:dim(pca)[1]]) +
  coord_fixed() + 
  xlab(paste0("PC1 (", prettyNum(summary(pca2)$importance[2,1]*100, digits =2),"%)")) +
  ylab(paste0("PC3 (", prettyNum(summary(pca2)$importance[2,3]*100, digits =2),"%)")) +
  labs(color = "Histotype", shape = "Group") + guides(size = FALSE) +
  ggtitle("OCRA - H3k27ac, H3k4me1, H3k4me3")
ggarrange(p3,p4)

















