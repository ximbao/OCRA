### Cell lines QC metrics plot ####
## example of how to get the information: grep -FR --include=*.frip.qc . > qc_frip.txt

library(tidyr)

nreads <- read.table("/media/Data01/Felipe/OvData/qc_nreads_all.txt", header= F)
nreads <- separate(nreads, V1, c("file", "n_reads"), sep = ":")
nreads$file <- gsub('^((.*?/){1}.*?)/.*', "\\1", nreads$file)
unames <- unique(nreads$file)
nreads$file <- make.unique(nreads$file)
nreads <- nreads[nreads$file %in% unames, ]
nreads$n_reads <- as.numeric(nreads$n_reads)

encodeqc <- read.table("/media/Data01/Felipe/OvData/qc_encodeMetrics.txt", header = F)
names(encodeqc) <- c("file", "NSC", "RSC", "Qtag")
encodeqc$file <- gsub('^((.*?/){1}.*?)/.*', "\\1", encodeqc$file)
encodeqc$file <- make.unique(encodeqc$file)
encodeqc <- encodeqc[encodeqc$file %in% unames, ]
write.table(encodeqc, row.names = F, sep = "\t", quote = F, file = "~/cellLine_QC.txt")

frip <- read.table("/media/Data01/Felipe/OvData/qc_frip.txt", header = F)
frip <- separate(frip, V1, c("file", "FRiP"), sep = ":")
frip$file <- gsub('^((.*?/){1}.*?)/.*', "\\1", frip$file)
frip$file <- make.unique(frip$file)
frip <- frip[frip$file %in% unames, ]
write.table(frip, row.names = F, sep = "\t", quote = F, file = "~/cellLine_QC_FRiP.txt")


library(ggplot2)

ggplot(nreads, aes(reorder(file, -n_reads), n_reads)) + geom_bar(stat="identity") +
  xlab("Sample") + ylab('Number of Reads') +
  labs(title= "QC - Number of Mapped reads - ChipSeq") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 20000000, linetype='dotted') +
  geom_text(aes(0, 20000000,label='Encode Standard - Acceptable', hjust = -0.1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(data=nreads,aes(reorder(file, -n_reads), n_reads, label = prettyNum(n_reads / 1000000, digits = 4)),vjust=0)

p1 <- ggplot(encodeqc, aes(NA, NSC)) + geom_boxplot(aes(alpha=0)) +
  geom_jitter(aes(size=2)) +
  theme(legend.position="none") +
  theme(axis.title.x =element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  geom_hline(yintercept = 1.05, linetype='dotted') +
  geom_text(aes(0, 1.06,label='Encode Standard - Acceptable', hjust = -0.1)) +
  labs(title = 'Normalized Strand Cross-correlation coefficient (NSC)') +
  theme(plot.title = element_text(hjust = 0.5)) + scale_y_log10()

p2 <- ggplot(encodeqc, aes(NA, RSC)) + geom_boxplot(aes(alpha=0)) +
  geom_jitter(aes(size=2)) +
  theme(legend.position="none") +
  theme(axis.title.x =element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  geom_hline(yintercept = 0.8, linetype='dotted') +
  geom_text(aes(0, 0.81,label='Encode Standard - Acceptable', hjust = -0.1)) +
  labs(title = 'Relative Strand Cross-correlation coefficient (RSC)') +
  theme(plot.title = element_text(hjust = 0.5)) 

p3 <- ggplot(frip, aes(NA, FRiP)) + geom_boxplot(aes(alpha=0)) +
  geom_jitter(aes(size=2)) +
  theme(legend.position="none") +
  theme(axis.title.x =element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  geom_hline(yintercept = 0.2, linetype='dotted') +
  #geom_hline(yintercept = 0.5, linetype='dotted') +
  geom_text(aes(0, 0.21,label='Encode Standard - Acceptable', hjust = -0.1)) +
  labs(title = 'Fraction of Reads in Peaks (FRiP)') +
  theme(plot.title = element_text(hjust = 0.5)) 

grid.arrange(p1, p2, p3, nrow = 1)

#######  tsv parse with all QC data ######
qc.table <- read.delim2("~/chip-seq-pipeline2/utils/qc_jsons_to_tsv/ocra_v1.txt", header = T, skip = 1, stringsAsFactors = F)
qc.table2 <- read.delim2("~/OCRA/QC/OCRA_QC.tsv", header = T, skip = 1, stringsAsFactors = F)
qc.table3 <- read.delim2("~/OCRA/QC/OCRA_redos_QC.tsv", header = T, skip = 1, stringsAsFactors = F)


aux <- qc.table2[c(1,5,9,16,18,20,24,26,30,34,38,42,46,48,50,52,54,58,60,64,66,68,69,71,73), 
                c(6,32,38,36,110,84,181:183,152)]
rownames(aux) <- NULL
aux2 <- read.delim2("~/chip-seq-pipeline2/utils/qc_jsons_to_tsv/ocra_v2.txt", header = T, skip = 1, stringsAsFactors = F)
aux2 <- aux2[c(1,5,9,13), c(6,32,38,36,110,84,181:183,152)]
rownames(aux2) <- NULL

ocra.qc <- rbind(aux, aux2)
ocra.qc$title <- gsub("_hg38", "", ocra.qc$title)
ocra.qc$mapped_pct <- (ocra.qc$mapped/ocra.qc$total)*100
ocra.qc$dupes_pct <- as.numeric(ocra.qc$dupes_pct)
ocra.qc$NSC <- as.numeric(ocra.qc$NSC)
ocra.qc$RSC <- as.numeric(ocra.qc$RSC)
ocra.qc$FRiP <- as.numeric(ocra.qc$FRiP)


p1 <- ggplot(ocra.qc, aes(reorder(title, total), total)) + geom_bar(stat="identity") +
  xlab("Sample") + ylab('Total Number of Reads (M)') +
  labs(title= "OCRA - Encode chipseq pipeline - Number of Reads") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 20000000, linetype='dotted', color= "yellow") +
  #geom_text(aes(1, 20000000,label='Encode Standard - Acceptable', hjust = -0.1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(data=ocra.qc,aes(reorder(title, -total), total, label = prettyNum(total / 1000000, digits = 4)),vjust=0) + coord_flip()

p2 <- ggplot(ocra.qc, aes(reorder(title, total), mapped_pct)) + geom_bar(stat="identity") +
  xlab("Sample") + ylab('% of mapped reads') +
  labs(title= "OCRA - Encode chipseq pipeline - Mapped Reads") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 80, linetype='dotted', color = "yellow") +
 # geom_text(aes(0, 80,label='Encode Standard - Acceptable', hjust = -0.1)) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(data=ocra.qc[order(ocra.qc$total), ]
            ,aes(title, mapped_pct, label = round(mapped_pct,2) ),vjust=0) 
grid.arrange(p1, p2, nrow = 1)

p3 <- ggplot(ocra.qc, aes(reorder(title, total), dupes_pct)) + geom_bar(stat="identity") +
  xlab("Sample") + ylab('% of dup reads') +
  labs(title= "OCRA - Encode chipseq pipeline - % Of Duplicate reads") +
  theme(plot.title = element_text(hjust = 0.5)) +
  #geom_hline(yintercept = 80, linetype='dotted', color = "yellow") +
  # geom_text(aes(0, 80,label='Encode Standard - Acceptable', hjust = -0.1)) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(data=ocra.qc, aes(reorder(title, total), dupes_pct, label = round(dupes_pct,2) ),vjust=0) 


p4 <- ggplot(ocra.qc, aes(reorder(title, total), mapped.2)) + geom_bar(stat="identity") +
  xlab("Sample") + ylab('Filtered mapped reads') +
  labs(title= "OCRA - Encode chipseq pipeline - Filtered mapped reads") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 20000000, linetype='dotted', color= "yellow") +
  #geom_text(aes(1, 20000000,label='Encode Standard - Acceptable', hjust = -0.1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(data=ocra.qc,aes(reorder(title, -total), mapped.2, label = prettyNum(mapped.2 / 1000000, digits = 4)),vjust=0) + coord_flip()
grid.arrange(p3, p4, nrow = 1)

p5 <- ggplot(ocra.qc, aes(reorder(title, total), NSC)) + geom_bar(stat="identity") +
  xlab("Sample") + ylab('Normalized Strand Cross-correlation coefficient') +
  labs(title= "OCRA - Encode chipseq pipeline - NSC ") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 1.05, linetype='dotted', color = "yellow") +
  #geom_text(aes(0, 80,label='Encode Standard - Acceptable', hjust = -0.1)) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(data=ocra.qc, aes(reorder(title, total), NSC, label = round(NSC,2) ),vjust=0) 

p6 <- ggplot(ocra.qc, aes(reorder(title, total), RSC)) + geom_bar(stat="identity") +
  xlab("Sample") + ylab('Relative Strand Cross-correlation coefficient') +
  labs(title= "OCRA - Encode chipseq pipeline - RSC ") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 0.8, linetype='dotted', color = "yellow") +
  #geom_text(aes(0, 80,label='Encode Standard - Acceptable', hjust = -0.1)) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(data=ocra.qc, aes(reorder(title, total), RSC, label = round(RSC,2) ),vjust=0) 
grid.arrange(p5, p6, nrow = 1)

p7 <- ggplot(ocra.qc, aes(reorder(title, total), FRiP)) + geom_bar(stat="identity") +
  xlab("Sample") + ylab('Fraction of reads in peaks') +
  labs(title= "OCRA - Encode chipseq pipeline - FRiP ") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 0.1, linetype='dotted', color = "yellow") +
  #geom_text(aes(0, 80,label='Encode Standard - Acceptable', hjust = -0.1)) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(data=ocra.qc, aes(reorder(title, total), FRiP, label = round(FRiP,2) ),vjust=0) 

p8 <- ggplot(ocra.qc, aes(reorder(title, total), Nt)) + geom_bar(stat="identity") +
  xlab("Sample") + ylab('Number of peaks') +
  labs(title= "OCRA - Encode chipseq pipeline - Number of peaks ") +
  theme(plot.title = element_text(hjust = 0.5)) +
  #geom_hline(yintercept = 0.1, linetype='dotted', color = "yellow") +
  #geom_text(aes(0, 80,label='Encode Standard - Acceptable', hjust = -0.1)) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(data=ocra.qc, aes(reorder(title, total), Nt, label = round(Nt,2) ),vjust=0) 
grid.arrange(p7, p8, nrow = 1)
# 








