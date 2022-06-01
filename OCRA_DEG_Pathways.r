###### Pathway analysis of OCRA RNASeq comparisons ######

library(magrittr)
library(clusterProfiler)

# Load DEEG Objects from ocra_rnaseq.R script 
load("~/OCRA/DEG_Objects_OCRA.Rda")

res #                                   OvCancer vs Normal                                                        #
res <- na.omit(res)
norm.can.path <- res[res$padj < 0.05, ]$log2FoldChange
## Name for each gene need to be in entrezid 

library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
genes.mart <- 
  getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), 
        values = substr(rownames(res[res$padj < 0.05, ]),1,15), mart = mart)

res$entrezid <- genes.mart$entrezgene_id[match(substr(rownames(res),1,15), genes.mart$ensembl_gene_id)]
names(norm.can.path) <- res[res$padj < 0.05,]$entrezid


gene.df <- bitr(names(norm.can.path), fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)
gene.df <- gene.df[gene.df$ENSEMBL %in% substr(rownames(res),1,15),]
aux <- res[substr(rownames(res),1,15) %in% gene.df$ENSEMBL, ]
gene.df$FC <- aux$log2FoldChange[match(substr(rownames(aux),1,15), gene.df$ENSEMBL)]

ego <- enrichGO(gene         = gene.df$ENSEMBL,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05)
head(ego)


ego.mf <-  enrichGO(gene         = gene.df$ENSEMBL,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "MF",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05)


ego.cc <-  enrichGO(gene         = gene.df$ENSEMBL,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "CC",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05)

norn.can.kk <- enrichKEGG(gene         = gene.df$ENTREZID,
                       organism     = 'hsa',
                       pvalueCutoff = 0.05, qvalueCutoff = 1)



#geneList = sort(norm.can.path, decreasing = TRUE)
e.go.bp <- setReadable(ego, 'org.Hs.eg.db', 'ENSEMBL')
e.go.mf <- setReadable(ego.mf, 'org.Hs.eg.db', 'ENSEMBL')
e.go.cc <- setReadable(ego.cc, 'org.Hs.eg.db', 'ENSEMBL')


write.table(e.go.bp, quote =F, row.names = F, sep = "\t", 
            file = "~/Documents/OCRA_Files/Pathway_Analysis/NormalvsOvCancer_GO_ORA_BP_table.tsv")
write.table(e.go.mf, quote =F, row.names = F, sep = "\t", 
            file = "~/Documents/OCRA_Files/Pathway_Analysis/NormalvsOvCancer_GO_ORA_MF_table.tsv")
write.table(e.go.cc, quote =F, row.names = F, sep = "\t", 
            file = "~/Documents/OCRA_Files/Pathway_Analysis/NormalvsOvCancer_GO_ORA_CC_table.tsv")


library(enrichplot)
pdf("~/Documents/OCRA_Files/Pathway_DotPlot_NormalsvsCancer_GOBP.pdf", 17, 15)
dotplot(ego, showCategory=50, title = "ORA - Normal vs Ovarian Cancer , GO BP")
dev.off()

pdf("~/Documents/OCRA_Files/Pathway_DotPlot_NormalsvsCancer_GOMF.pdf", 17, 15)
dotplot(ego.mf, showCategory=50,  title = "ORA - Normal vs Ovarian Cancer , GO MF")
dev.off()

pdf("~/Documents/OCRA_Files/Pathway_DotPlot_NormalsvsCancer_GOCC.pdf", 17, 15)
dotplot(ego.cc, showCategory = 50, title = "ORA - Normal vs Ovarian Cancer , GO CC")
dev.off()

#Goplots
pdf("~/Documents/OCRA_Files/Pathway_GOPlot_NormalsvsCancer_GOBP.pdf", 20, 15)
goplot(ego, showCategory=20, title = "ORA - Normals vs ovCancer , GO BP")
dev.off()

pdf("~/Documents/OCRA_Files/Pathway_GOPlot_NormalsvsCancer_GOMF.pdf", 20, 15)
goplot(ego.mf, showCategory=20, title = "ORA - Normals vs ovCancer , GO MF")
dev.off()

pdf("~/Documents/OCRA_Files/Pathway_GOPlot_NormalsvsCancer_GOCC.pdf", 20, 15)
goplot(ego.cc, showCategory=20, title = "ORA - Normals vs ovCancer , GO CC")
dev.off()


#                                           HGSOC vs FTSEC                                                    #
res2

res2 <- na.omit(res2)
hg.ft.path <- res2[res2$padj < 0.05, ]$log2FoldChange
## Name for each gene need to be in entrezid 

#mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
genes.mart <- 
  getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), 
        values = substr(rownames(res2[res2$padj < 0.05, ]),1,15), mart = mart)

res2$entrezid <- genes.mart$entrezgene_id[match(substr(rownames(res2),1,15), genes.mart$ensembl_gene_id)]
names(hg.ft.path) <- res2[res2$padj < 0.05,]$entrezid


gene.df <- bitr(names(hg.ft.path), fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)
# gene.df <- gene.df[gene.df$ENSEMBL %in% substr(rownames(res),1,15),]
# aux <- res[substr(rownames(res),1,15) %in% gene.df$ENSEMBL, ]
# gene.df$FC <- aux$log2FoldChange[match(substr(rownames(aux),1,15), gene.df$ENSEMBL)]

hg.ft.ego.bp <- enrichGO(gene         = gene.df$ENSEMBL,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENSEMBL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05)
head(hg.ft.ego.bp)


hg.ft.ego.mf <-  enrichGO(gene         = gene.df$ENSEMBL,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "MF",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05)


hg.ft.ego.cc <-  enrichGO(gene         = gene.df$ENSEMBL,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "CC",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05)

hg.ft.kk <- enrichKEGG(gene         = gene.df$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05, qvalueCutoff = 1)


#geneList = sort(norm.can.path, decreasing = TRUE)
e.hg.ft.ego.bp <- setReadable(hg.ft.ego.bp, 'org.Hs.eg.db', 'ENSEMBL')
e.hg.ft.ego.mf <- setReadable(hg.ft.ego.mf, 'org.Hs.eg.db', 'ENSEMBL')
e.hg.ft.ego.cc <- setReadable(hg.ft.ego.cc, 'org.Hs.eg.db', 'ENSEMBL')
e.gh.ft.kk <-  setReadable(hg.ft.kk, 'org.Hs.eg.db', 'ENTREZID')
  
write.table(e.hg.ft.ego.bp, quote =F, row.names = F, sep = "\t", 
            file = "~/Documents/OCRA_Files/Pathway_Analysis/HGSOCvsFTSEC_GO_ORA_BP_table.tsv")
write.table(e.hg.ft.ego.mf, quote =F, row.names = F, sep = "\t", 
            file = "~/Documents/OCRA_Files/Pathway_Analysis/HGSOCvsFTSEC_GO_ORA_MF_table.tsv")
write.table(e.hg.ft.ego.cc, quote =F, row.names = F, sep = "\t", 
            file = "~/Documents/OCRA_Files/Pathway_Analysis/HGSOCvsFTSEC_GO_ORA_CC_table.tsv")


pdf("~/Documents/OCRA_Files/Pathway_DotPlot_HGSOCvsFTSEC_GOBP.pdf", 17, 15)
dotplot(hg.ft.ego.bp, showCategory=50, title = "ORA - HGSOC vs FTSEC , GO BP")
dev.off()

pdf("~/Documents/OCRA_Files/Pathway_DotPlot_HGSOCvsFTSEC_GOMF.pdf", 17, 15)
dotplot(hg.ft.ego.mf, showCategory=50, title = "ORA - HGSOC vs FTSEC , GO MF")
dev.off()

pdf("~/Documents/OCRA_Files/Pathway_DotPlot_HGSOCvsFTSEC_GOCC.pdf", 17, 15)
dotplot(hg.ft.ego.cc, showCategory=50, title = "ORA - HGSOC vs FTSEC , GO CC")
dev.off()

#Goplots
pdf("~/Documents/OCRA_Files/Pathway_GOPlot_HGSOCvsFTSEC_GOBP.pdf", 20, 15)
goplot(hg.ft.ego.bp, showCategory=20, title = "ORA - HGSOC vs FTSEC , GO BP")
dev.off()

pdf("~/Documents/OCRA_Files/Pathway_GOPlot_HGSOCvsFTSEC_GOMF.pdf", 20, 15)
goplot(hg.ft.ego.mf, showCategory=20, title = "ORA - HGSOC vs FTSEC , GO MF")
dev.off()

pdf("~/Documents/OCRA_Files/Pathway_GOPlot_HGSOCvsFTSEC_GOCC.pdf", 20, 15)
goplot(hg.ft.ego.cc, showCategory=20, title = "ORA - HGSOC vs FTSEC , GO CC")
dev.off()


#                                           IOSE vs HGSOC                                                    #
res3

res3 <- na.omit(res3)
io.hg.path <- res3[res3$padj < 0.05, ]$log2FoldChange
## Name for each gene need to be in entrezid 

#mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
genes.mart <- 
  getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), 
        values = substr(rownames(res3[res3$padj < 0.05, ]),1,15), mart = mart)

res3$entrezid <- genes.mart$entrezgene_id[match(substr(rownames(res3),1,15), genes.mart$ensembl_gene_id)]
names(io.hg.path) <- res3[res3$padj < 0.05,]$entrezid


gene.df <- bitr(names(io.hg.path), fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)
# gene.df <- gene.df[gene.df$ENSEMBL %in% substr(rownames(res),1,15),]
# aux <- res[substr(rownames(res),1,15) %in% gene.df$ENSEMBL, ]
# gene.df$FC <- aux$log2FoldChange[match(substr(rownames(aux),1,15), gene.df$ENSEMBL)]

io.hg.ego.bp <- enrichGO(gene         = gene.df$ENSEMBL,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)
head(io.hg.ego.bp)


io.hg.ego.mf <-  enrichGO(gene         = gene.df$ENSEMBL,
                          OrgDb         = org.Hs.eg.db,
                          keyType       = 'ENSEMBL',
                          ont           = "MF",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05)


io.hg.ego.cc <-  enrichGO(gene         = gene.df$ENSEMBL,
                          OrgDb         = org.Hs.eg.db,
                          keyType       = 'ENSEMBL',
                          ont           = "CC",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05)

io.hg.kk <- enrichKEGG(gene         = gene.df$ENTREZID,
                       organism     = 'hsa',
                       pvalueCutoff = 0.05, qvalueCutoff = 1)


#geneList = sort(norm.can.path, decreasing = TRUE)
e.io.hg.ego.bp <- setReadable(io.hg.ego.bp, 'org.Hs.eg.db', 'ENSEMBL')
e.io.hg.ego.mf <- setReadable(io.hg.ego.mf, 'org.Hs.eg.db', 'ENSEMBL')
e.io.hg.ego.cc <- setReadable(io.hg.ego.cc, 'org.Hs.eg.db', 'ENSEMBL')
e.io.hg.kk <-  setReadable(io.hg.kk, 'org.Hs.eg.db', 'ENTREZID')

write.table(e.io.hg.ego.bp, quote =F, row.names = F, sep = "\t", 
            file = "~/Documents/OCRA_Files/Pathway_Analysis/IOSEvsHGSOC_GO_ORA_BP_table.tsv")
write.table(e.io.hg.ego.mf, quote =F, row.names = F, sep = "\t", 
            file = "~/Documents/OCRA_Files/Pathway_Analysis/IOSEvsHGSOC_GO_ORA_MF_table.tsv")
write.table(e.io.hg.ego.cc, quote =F, row.names = F, sep = "\t", 
            file = "~/Documents/OCRA_Files/Pathway_Analysis/IOSEvsHGSOC_GO_ORA_CC_table.tsv")


pdf("~/Documents/OCRA_Files/Pathway_DotPlot_IOSEvsHGSOC_GOBP.pdf", 17, 15)
dotplot(io.hg.ego.bp, showCategory=50, title = "ORA - IOSE vs HGSOC , GO BP")
dev.off()

pdf("~/Documents/OCRA_Files/Pathway_DotPlot_IOSECvsHGSOC_GOMF.pdf", 17, 15)
dotplot(io.hg.ego.mf, showCategory=50, title = "ORA - IOSE vs HGSOC , GO MF")
dev.off()

pdf("~/Documents/OCRA_Files/Pathway_DotPlot_IOSEvsHGSOC_GOCC.pdf", 17, 15)
dotplot(io.hg.ego.cc, showCategory=50, title = "ORA - IOSE vs HGSOC , GO CC")
dev.off()

#Goplots
pdf("~/Documents/OCRA_Files/Pathway_GOPlot_IOSEvsHGSOC_GOBP.pdf", 20, 15)
goplot(io.hg.ego.bp, showCategory=20, title = "ORA - IOSE vs FTSEC , GO BP")
dev.off()

pdf("~/Documents/OCRA_Files/Pathway_GOPlot_IOSEvsHGSOC_GOMF.pdf", 20, 15)
goplot(io.hg.ego.mf, showCategory=20, title = "ORA - IOSE vs HGSOC , GO BP")
dev.off()

pdf("~/Documents/OCRA_Files/Pathway_GOPlot_IOSEvsHGSOC_GOCC.pdf", 20, 15)
goplot(io.hg.ego.cc, showCategory=20, title = "ORA - IOSE vs HGSOC , GO BP")
dev.off()




#                                           IOSE vs FTSEC                                                    #
res4

res4 <- na.omit(res4)
io.ft.path <- res4[res4$padj < 0.05, ]$log2FoldChange
## Name for each gene need to be in entrezid 

#mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
genes.mart <- 
  getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), 
        values = substr(rownames(res4[res4$padj < 0.05, ]),1,15), mart = mart)

res4$entrezid <- genes.mart$entrezgene_id[match(substr(rownames(res4),1,15), genes.mart$ensembl_gene_id)]
names(io.ft.path) <- res4[res4$padj < 0.05,]$entrezid


gene.df <- bitr(names(io.ft.path), fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)
# gene.df <- gene.df[gene.df$ENSEMBL %in% substr(rownames(res),1,15),]
# aux <- res[substr(rownames(res),1,15) %in% gene.df$ENSEMBL, ]
# gene.df$FC <- aux$log2FoldChange[match(substr(rownames(aux),1,15), gene.df$ENSEMBL)]

io.ft.ego.bp <- enrichGO(gene         = gene.df$ENSEMBL,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)


io.ft.ego.mf <-  enrichGO(gene         = gene.df$ENSEMBL,
                          OrgDb         = org.Hs.eg.db,
                          keyType       = 'ENSEMBL',
                          ont           = "MF",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05)


io.ft.ego.cc <-  enrichGO(gene         = gene.df$ENSEMBL,
                          OrgDb         = org.Hs.eg.db,
                          keyType       = 'ENSEMBL',
                          ont           = "CC",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05)

io.ft.kk <- enrichKEGG(gene         = gene.df$ENTREZID,
                       organism     = 'hsa',
                       pvalueCutoff = 0.05, qvalueCutoff = 1)


#geneList = sort(norm.can.path, decreasing = TRUE)
e.io.ft.ego.bp <- setReadable(io.ft.ego.bp, 'org.Hs.eg.db', 'ENSEMBL')
e.io.ft.ego.mf <- setReadable(io.ft.ego.mf, 'org.Hs.eg.db', 'ENSEMBL')
e.io.ft.ego.cc <- setReadable(io.ft.ego.cc, 'org.Hs.eg.db', 'ENSEMBL')
e.io.ft.kk<-  setReadable(io.ft.kk, 'org.Hs.eg.db', 'ENTREZID')

write.table(e.io.ft.ego.bp, quote =F, row.names = F, sep = "\t", 
            file = "~/Documents/OCRA_Files/Pathway_Analysis/IOSEvsFTSEC_GO_ORA_BP_table.tsv")
write.table(e.io.ft.ego.mf, quote =F, row.names = F, sep = "\t", 
            file = "~/Documents/OCRA_Files/Pathway_Analysis/IOSEvsFTSEC_GO_ORA_MF_table.tsv")
write.table(e.io.ft.ego.cc, quote =F, row.names = F, sep = "\t", 
            file = "~/Documents/OCRA_Files/Pathway_Analysis/IOSEvsFTSEC_GO_ORA_CC_table.tsv")


pdf("~/Documents/OCRA_Files/Pathway_DotPlot_IOSEvsFTSEC_GOBP.pdf", 17, 15)
dotplot(io.ft.ego.bp, showCategory=50, title = "ORA - IOSE vs FTSEC , GO BP")
dev.off()

pdf("~/Documents/OCRA_Files/Pathway_DotPlot_IOSECvsFTSEC_GOMF.pdf", 17, 15)
dotplot(io.ft.ego.mf, showCategory=50, title = "ORA - IOSE vs FTSEC , GO MF")
dev.off()

pdf("~/Documents/OCRA_Files/Pathway_DotPlot_IOSEvsFTSEC_GOCC.pdf", 17, 15)
dotplot(io.ft.ego.cc, showCategory=50, title = "ORA - IOSE vs FTSEC , GO CC")
dev.off()

#Goplots
pdf("~/Documents/OCRA_Files/Pathway_GOPlot_IOSEvsFTSEC_GOBP.pdf", 20, 15)
goplot(io.ft.ego.bp, showCategory=20, title = "ORA - IOSE vs FTSEC , GO BP")
dev.off()

pdf("~/Documents/OCRA_Files/Pathway_GOPlot_IOSEvsFTSEC_GOMF.pdf", 20, 15)
goplot(io.ft.ego.mf, showCategory=20, title = "ORA - IOSE vs FTSEC , GO MF")
dev.off()

pdf("~/Documents/OCRA_Files/Pathway_GOPlot_IOSEvsFTSEC_GOCC.pdf", 20, 15)
goplot(io.ft.ego.cc, showCategory=20, title = "ORA - IOSE vs FTSEC , GO CC")
dev.off()





###

aux <- read_excel("~/Documents/OCRA_Files/Plummer_OCRA_SciAdv_Suppl_2021.xlsx", sheet = 21, skip = 1)





