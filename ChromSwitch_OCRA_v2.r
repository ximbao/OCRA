###### Chromswitch part 2 - using chromhmm states and different query ######
library(chromswitch)


## prepare states as GRanges list
hmm.states <- dir("~/Documents/OCRA_ChromHMM_hg19/")
hmm.states <- dir("~/Documents/OCRA_Files/ChromHMM_Edited/")

## Consensus states
ft.cons <- fread("~/Documents/OCRA_Files/Consensus_states/FT_consensus_UCSC.bed", skip = 1)
hg.cons <- fread("~/Documents/OCRA_Files/Consensus_states/HGSOC_consensus_UCSC.bed", skip = 1)
iose.cons <- fread("~/Documents/OCRA_Files/Consensus_states/IOSE_consensus_UCSC.bed", skip = 1)

ft.cons <- ft.cons[, 1:4]
colnames(ft.cons) <-  c("chr", "start", "end", "state_n")
ft.cons <- makeGRangesFromDataFrame(na.omit(ft.cons), keep.extra.columns = T)
hg.cons <- hg.cons[, 1:4]
colnames(hg.cons) <-  c("chr", "start", "end", "state_n")
hg.cons <- makeGRangesFromDataFrame(na.omit(hg.cons), keep.extra.columns = T)
iose.cons <- iose.cons[, 1:4]
colnames(iose.cons) <-  c("chr", "start", "end", "state_n")
iose.cons <- makeGRangesFromDataFrame(na.omit(iose.cons), keep.extra.columns = T)


cons.states <- GRangesList(ft.cons, hg.cons, iose.cons)
names(cons.states) <- c("FTSEC_Consensus", "HGSOC_Consensus", "iOSE_Consensus")

#aux <- fread(paste0("~/Documents/OCRA_ChromHMM_hg19/", hmm.states[1]), skip = 1)
aux <- fread(paste0("~/Documents/OCRA_Files/ChromHMM_Edited/", hmm.states[1]), skip = 1) ## cleaner chromhmm calls
aux <- aux[, 1:4]
colnames(aux) <- c("chr", "start", "end", "state_n")
list.states <- GRangesList(makeGRangesFromDataFrame(aux, keep.extra.columns = T))
names(list.states) <- hmm.states[1]
w <- which(list.states[[1]]$state_n %in% c("TRS", "CTCF"))
list.states[[1]] <- list.states[[i]][-c(w)]
for( i in 2:length(hmm.states)) {
 # aux <- fread(paste0("~/Documents/OCRA_ChromHMM_hg19/", hmm.states[i]), skip = 1)
  aux <- fread(paste0("~/Documents/OCRA_Files/ChromHMM_Edited/", hmm.states[i]), skip = 1)
  aux <- aux[, 1:4]
  colnames(aux) <- c("chr", "start", "end", "state_n")
  list.states[[i]] <- makeGRangesFromDataFrame(aux, keep.extra.columns = T)
  names(list.states)[i] <- hmm.states[i]
  w <- which(list.states[[i]]$state_n %in% c("TRS", "CTCF"))
  list.states[[i]] <- list.states[[i]][-c(w)]
}

metadata <- data.frame(Sample = names(list.states), Condition = NA)
metadata$Condition <- c("HGSOC", "Endometrioid", "Muc", "CCOC", rep("FTSEC", 3), "Muc",
                        "HGSOC", rep("IOSE", 2), "CCOC", "HGSOC", "Muc", "HGSOC", "CCOC", "HGSOC", "LGSOC")

cons.metadata <- data.frame(Sample = names(cons.states), Condition = c("FTSEC", "HGSOC", "IOSE"))

## query will be DEG from ocra_rnaseq.R 
deg <- read.table("~/Documents/OCRA_Files/Table_DEG_FTSEC_vs_HGSOC.txt")

library(biomaRt)
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
bm <- getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol'),
      filters=c('hgnc_symbol'),
      values=deg[deg$STATUS %in% "Down", ]$GeneSymbol,
      mart=ensembl)
bm2 <- getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'ensembl_gene_id','hgnc_symbol'),
            filters=c('ensembl_gene_id'),
            values= gsub("\\.[0-9]","", deg[deg$STATUS %in% "Up", ] %>% rownames),
            mart=ensembl)
bm3 <- getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol'),
            filters=c('hgnc_symbol'),
            values=deg$GeneSymbol,
            mart=ensembl)

genes.ft.reg <- bm
colnames(genes.ft.reg) <- c("chr", "start", "end", "Gene")
genes.ft.reg$chr <- paste0("chr", genes.ft.reg$chr)
genes.ft.reg <- makeGRangesFromDataFrame(genes.ft.reg, keep.extra.columns = T)

genes.hg.reg <- bm2
colnames(genes.hg.reg) <- c("chr", "start", "end", "Gene")
genes.hg.reg$chr <- paste0("chr", genes.hg.reg$chr)
genes.hg.reg <- genes.hg.reg[, -5]
genes.hg.reg <- makeGRangesFromDataFrame(genes.hg.reg, keep.extra.columns = T)

genes.deg.all <- bm3
colnames(genes.deg.all) <- c("chr", "start", "end", "Gene")
genes.deg.all$chr <- paste0("chr", genes.deg.all$chr)
genes.deg.all <- makeGRangesFromDataFrame(genes.deg.all, keep.extra.columns = T)


## HGSOC vs FTSEC (no ft282 and only 2 hgsoc lines) ##
deg.new <- read.table("~/Documents/OCRA_Files/DEG_HGSOC_FTSEC_no282.txt")
deg.new.hg <- deg.new[deg.new$STATUS %in% "Up",]
deg.new.ft <- deg.new[deg.new$STATUS %in% "Down",]

# get gene coordinates from hg19 gtf annotation file
gtf.hg19 <- readGFF("~/reference_files/gencode.v34lift37.annotation.gtf.gz")
gtf.hg19 <- gtf.hg19[gtf.hg19$type %in% "gene",]
gtf.hg19$gene_id <- substr(gtf.hg19$gene_id, 1, 15)
aux <- gtf.hg19[gtf.hg19$gene_id %in% substr(rownames(deg.new.hg),1,15),]
aux <- aux[, c(1,4:5,7,11)]
deg.new.hg.gr <- makeGRangesFromDataFrame(aux, keep.extra.columns = T)

aux <- gtf.hg19[gtf.hg19$gene_id %in% substr(rownames(deg.new.ft),1,15),]
aux <- aux[, c(1,4:5,7,11)]
deg.new.ft.gr <- makeGRangesFromDataFrame(aux, keep.extra.columns = T)



# get gene names and coordnates from the rnaseq object ocra.rnaseq
# ocra.rna[ocra.rna$Geneid %in% rownames(deg.new),] -> aux
# ocra.rna[ocra.rna$Geneid %in% rownames(deg.new.hg),] -> aux
# 
# aux <- aux[, c(2:5)]
# aux$gene <- rownames(aux)
# deg.new.gr <- makeGRangesFromDataFrame(aux, keep.extra.columns = T)
# deg.new.hg.gr <- makeGRangesFromDataFrame(aux, keep.extra.columns = T)


### all genes as quuery ###
gtf.hg19 <- readGFF("~/reference_files/gencode.v34lift37.annotation.gtf.gz")
gtf.hg19 <- gtf.hg19[gtf.hg19$type %in% "gene" , ]
gtf.hg19 <- gtf.hg19[gtf.hg19$gene_type %in% c("protein_coding", "lncRNA", "miRNA"),]

query.all.genes <- gtf.hg19[, c(1,4:5,7,11)]
colnames(query.all.genes) <- c("chr","start", "end", "strand")
query.all.genes <- makeGRangesFromDataFrame(query.all.genes, keep.extra.columns = T)

library(BiocParallel)
register(default = T, BPPARAM = DoparParam()) ## to avoid error in parallel processing

# ft genes
out.summ <- callSummary(query = genes.ft.reg,       # Input 1: Query regions
                       # metadata = metadata[c(1,5,7,9,13,15,17),], # Input 2: Metadata 
                        metadata[metadata$Condition %in% c("FTSEC", "HGSOC"),],
                       # peaks = list.states[metadata$Sample[c(1,5,7,9,13,15,17)],],       # Input 3: Peaks
                        peaks = list.states[metadata$Sample[metadata$Condition %in% c("FTSEC", "HGSOC")]],
                        mark = "Chromatin States",  # Arbitrary string describing the data type
                        heatmap = F,
                        normalize = TRUE)   

out.summ.cons <- callSummary(query = genes.deg.all,       # Input 1: Query regions
                        #                                 , # Input 2: Metadata 
                        metadata =  cons.metadata,
                        
                       # metadata =  cons.metadata[cons.metadata$Condition %in% c("FTSEC", "HGSOC"),],
                        #                                    # Input 3: Peaks
                       # peaks = cons.states[cons.metadata$Sample[cons.metadata$Condition %in% c("FTSEC", "HGSOC")]],
                        peaks = cons.states,
                       
                        mark = "Chromatin States",  # Arbitrary string describing the data type
                        heatmap = F,
                        normalize = F)   

# hg genes
out.summ2 <- callSummary(query = genes.hg.reg,       # Input 1: Query regions
                        # metadata = metadata[c(1,5,7,9,13,15,17),], # Input 2: Metadata 
                        metadata[metadata$Condition %in% c("FTSEC", "HGSOC"),],
                        # peaks = list.states[metadata$Sample[c(1,5,7,9,13,15,17)],],       # Input 3: Peaks
                        peaks = list.states[metadata$Sample[metadata$Condition %in% c("FTSEC", "HGSOC")]],
                        mark = "Chromatin States",  # Arbitrary string describing the data type
                        heatmap = F,
                        normalize = TRUE)   

out.summ.all <- callSummary(query = deg.new.gr,       # Input 1: Query regions
                         # metadata = metadata[c(1,5,7,9,13,15,17),], # Input 2: Metadata 
                         metadata[metadata$Condition %in% c("FTSEC", "HGSOC"),][c(2,4,6,8),],
                         # peaks = list.states[metadata$Sample[c(1,5,7,9,13,15,17)],],       # Input 3: Peaks
                         peaks = list.states[metadata$Sample[metadata$Condition %in% c("FTSEC", "HGSOC")]][c(2,4,6,8),],
                         mark = "Chromatin States",  # Arbitrary string describing the data type
                         heatmap = F,
                         normalize = F)  
#break down new deg by hg and ft genes
out.summ.new.hg <- callSummary(query = deg.new.hg.gr,       # Input 1: Query regions
                            # metadata = metadata[c(1,5,7,9,13,15,17),], # Input 2: Metadata 
                            metadata[metadata$Condition %in% c("FTSEC", "HGSOC"),][c(2,4,6,8),],
                            # peaks = list.states[metadata$Sample[c(1,5,7,9,13,15,17)],],       # Input 3: Peaks
                            peaks = list.states[metadata$Sample[metadata$Condition %in% c("FTSEC", "HGSOC")]][c(2,4,6,8),],
                            mark = "Chromatin States",  # Arbitrary string describing the data type
                            heatmap = F,
                            normalize = F)   
write.table(out.summ.new.hg, quote = F, row.names = F, sep = "\t", file = '~/Documents/OCRA_Files/ChromSwitch_DEG_HGSOC_new.txt')

out.summ.new.ft <- callSummary(query = deg.new.ft.gr,       # Input 1: Query regions
                               # metadata = metadata[c(1,5,7,9,13,15,17),], # Input 2: Metadata 
                               metadata[metadata$Condition %in% c("FTSEC", "HGSOC"),][c(2,4,6,8),],
                               # peaks = list.states[metadata$Sample[c(1,5,7,9,13,15,17)],],       # Input 3: Peaks
                               peaks = list.states[metadata$Sample[metadata$Condition %in% c("FTSEC", "HGSOC")]][c(2,4,6,8),],
                               mark = "Chromatin States",  # Arbitrary string describing the data type
                               heatmap = F,
                               normalize = F)   

write.table(out.summ.new.ft, quote = F, row.names = F, sep = "\t", file = '~/Documents/OCRA_Files/ChromSwitch_DEG_FTSEC_new.txt')


out.bin.all <- callBinary(query = deg.new.gr,       # Input 1: Query regions
                            # metadata = metadata[c(1,5,7,9,13,15,17),], # Input 2: Metadata 
                            metadata[metadata$Condition %in% c("FTSEC", "HGSOC"),][c(2,4,6,8),],
                            # peaks = list.states[metadata$Sample[c(1,5,7,9,13,15,17)],],       # Input 3: Peaks
                            peaks = list.states[metadata$Sample[metadata$Condition %in% c("FTSEC", "HGSOC")]][c(2,4,6,8),]
                            )   
write.table(out.bin.all, quote = F, sep = "\t",
            row.names = F, file = "~/Documents/OCRA_Files/ChromSwitch_newDEG_FT_HG_binary.txt")

out.summ.all <- callSummary(query = query.all.genes,       # Input 1: Query regions
                         # metadata = metadata[c(1,5,7,9,13,15,17),], # Input 2: Metadata 
                         metadata[metadata$Condition %in% c("FTSEC", "HGSOC"),][c(2,4,6,8),],
                         # peaks = list.states[metadata$Sample[c(1,5,7,9,13,15,17)],],       # Input 3: Peaks
                         peaks = list.states[metadata$Sample[metadata$Condition %in% c("FTSEC", "HGSOC")]][c(2,4,6,8),],
                         mark = "Chromatin States",  # Arbitrary string describing the data type
                         heatmap = F,
                         normalize = TRUE)   

write.table(out.summ, quote = F, row.names = F, file = '~/Documents/OCRA_Files/ChromSwitch_DEG_FT_Summary.txt')
write.table(out.summ2, quote = F, row.names = F, file = '~/Documents/OCRA_Files/ChromSwitch_DEG_HG_Summary.txt')
write.table(out.summ3, quote = F, row.names = F, file = '~/Documents/OCRA_Files/ChromSwitch_DEG_All_Summary.txt')



out.summ.all <- callSummary(query = query.all.genes,       # Input 1: Query regions
                        metadata = metadata[metadata$Condition %in% c("FTSEC", "HGSOC"),], # Input 2: Metadata dataframe
                        peaks = list.states[metadata$Sample[metadata$Condition %in% c("FTSEC", "HGSOC")]],       # Input 3: Peaks
                        mark = "Chromatin States",  # Arbitrary string describing the data type
                        heatmap = F)

gr.all.switch <- out.summ.all
gr.all.switch <- separate(gr.all.switch, "query", into = c("start", "end"), sep = "-", remove = F)
gr.all.switch <- separate(gr.all.switch, "start", into = c('chr', 'start'), sep = ':', remove = T)
gr.all.switch <- gr.all.switch[, -1]
gr.all.switch <- makeGRangesFromDataFrame(gr.all.switch, keep.extra.columns = T)
###now annotate which switch with a gene near by
grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                 host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
genes_all <- getBM(
  attributes=c("hgnc_symbol","entrezgene_id","chromosome_name","start_position","end_position"),
  mart = mart)
colnames(genes_all)[3:5] <- c("chr", "start", "end")
genes_all$chr <- paste0("chr", genes_all$chr)

genes_all <- makeGRangesFromDataFrame(genes_all, keep.extra.columns = T)
genes_all <- genes_all[seqnames(genes_all) %in% seqnames(gr.all.switch)] 
gr.all.switch <- gr.all.switch[seqnames(gr.all.switch) %in% seqnames(genes_all)]

#ovp <- findOverlaps(gr.all.switch, genes(TxDb.Hsapiens.UCSC.hg19.knownGene))
ovp <- findOverlaps(gr.all.switch, genes_all)
gr.all.switch <- gr.all.switch[queryHits(ovp)]
#gr.all.switch$Gene <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)$gene_id[subjectHits(ovp)]
gr.all.switch$Gene <- genes_all$hgnc_symbol[subjectHits(ovp)]
## filter by Consensus / cluster quality
gr.all.switch.filter <- gr.all.switch[gr.all.switch$Consensus >= 0.75, ]
all.switch.filter <- as.data.frame(gr.all.switch.filter)
all.switch.filter <- all.switch.filter[order(all.switch.filter$Consensus, decreasing = T), ] 
#gr.all.switch.filter <- gr.all.switch.filter[sort(gr.all.switch.filter$Consensus, decreasing = T)]
write.table(all.switch.filter, quote =F, row.names = F, sep = "\t", file = "~/Documents/OCRA_Files/ChromSwitch_All_Genes_Query.txt")

out.summ2 <- callSummary(query = genes.hg.reg ,       # Input 1: Query regions
                        metadata = metadata[c(1,5,7,9,13,15,17),], # Input 2: Metadata dataframe
                        peaks =  list.states[metadata$Sample[c(1,5,7,9,13,15,17)],],       # Input 3: Peaks
                        mark = "Chromatin States",  # Arbitrary string describing the data type
                        heatmap = F)   



a <- callBinary(query = genes.ft.reg,       # Input 1: Query regions
                metadata =metadata[metadata$Condition %in% c("FTSEC", "HGSOC"),],      # Input 2: Metadata dataframe
                peaks = list.states[metadata$Sample[metadata$Condition %in% c("FTSEC", "HGSOC")]],
                heatmap =F)      # Input 3: Peaks


write.table(a, quote = F, row.names = F, file = '~/Documents/OCRA_Files/ChromSwitch_DEG_FT_binary.txt')
## MUC

bm3 <- getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol'),
             filters=c('hgnc_symbol'),
             values=c(genes.muc,genes.muc.hg),
             mart=ensembl)

genes.muc.reg <- bm3
colnames(genes.muc.reg) <- c("chr", "start", "end", "Gene")
genes.muc.reg$chr <- paste0("chr", genes.muc.reg$chr)
genes.muc.reg <- makeGRangesFromDataFrame(genes.muc.reg, keep.extra.columns = T)



out.summ.muc <- callSummary(query = genes.muc.reg,       # Input 1: Query regions
                        metadata[metadata$Condition %in% c("Muc", "HGSOC"),],
                        # peaks = list.states[metadata$Sample[c(1,5,7,9,13,15,17)],],       # Input 3: Peaks
                        peaks = list.states[metadata$Sample[metadata$Condition %in% c("Muc", "HGSOC")]],
                        mark = "Chromatin States",  # Arbitrary string describing the data type
                        heatmap = F,
                        normalize = TRUE)  


##### iOSE
deg <- read.table("~/Documents/OCRA_Files/Table_DEG_OSEC_vs_HGSOC.txt")

bm4 <- getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'ensembl_gene_id','hgnc_symbol'),
             filters=c('ensembl_gene_id'),
             values=gsub("\\.[0-9]","",rownames(deg)),
             mart=ensembl)

genes.iose.reg <- bm4
colnames(genes.iose.reg) <- c("chr", "start", "end", "GeneID", "GeneSymbol")
genes.iose.reg$chr <- paste0("chr", genes.iose.reg$chr)
genes.iose.reg <- makeGRangesFromDataFrame(genes.iose.reg, keep.extra.columns = T)



out.summ.ose <- callSummary(query = genes.iose.reg,       # Input 1: Query regions
                            metadata[metadata$Condition %in% c("IOSE", "HGSOC"),],
                            # peaks = list.states[metadata$Sample[c(1,5,7,9,13,15,17)],],       # Input 3: Peaks
                            peaks = list.states[metadata$Sample[metadata$Condition %in% c("IOSE", "HGSOC")]],
                            mark = "Chromatin States",  # Arbitrary string describing the data type
                            heatmap = F,
                            normalize = TRUE)  
write.table(out.summ.ose, quote = F, row.names = F, file = '~/Documents/OCRA_Files/ChromSwitch_DEG_iOSE_HGSOC_Summary.txt')


##### CCOC
bm5 <- getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol'),
             filters=c('hgnc_symbol'),
             values=genes.ccoc,
             mart=ensembl)

genes.ccoc.reg <- bm5
colnames(genes.ccoc.reg) <- c("chr", "start", "end", "Gene")
genes.ccoc.reg$chr <- paste0("chr", genes.ccoc.reg$chr)
genes.ccoc.reg <- makeGRangesFromDataFrame(genes.ccoc.reg, keep.extra.columns = T)



out.summ.ccoc <- callSummary(query = genes.ccoc.reg,       # Input 1: Query regions
                            metadata[metadata$Condition %in% c("CCOC", "Endometrioid"),],
                            # peaks = list.states[metadata$Sample[c(1,5,7,9,13,15,17)],],       # Input 3: Peaks
                            peaks = list.states[metadata$Sample[metadata$Condition %in% c("CCOC", "Endometrioid")]],
                            mark = "Chromatin States",  # Arbitrary string describing the data type
                            heatmap = F,
                            normalize = TRUE)  



###











