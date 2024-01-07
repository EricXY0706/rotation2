# https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

library("DESeq2")
library('RColorBrewer')
library('ggplot2')
library('pheatmap')
library(ggsci)
library("genefilter")
library('data.table')
library('stringr')
library(dplyr)
library(valr)
library('tidyverse')
library(EnhancedVolcano)
library(ggpubr)

##### RepeatMasker RawCount Data ######
setwd(paste0('/analysis2/xuy/RNA-seq/E250005786/5-RawCounts/RepeatMasker'))
FILE<-dir(pattern="*RawCounts*")
FILE<-FILE[grep("RepeatMasker",FILE, value = F)]
RawCountR<-read.table(FILE[1], sep="\t", stringsAsFactors=FALSE, header=TRUE)
for (i in FILE) {RawCountR<-cbind(RawCountR, read.table(i, sep="\t", stringsAsFactors=FALSE, header=TRUE)[,2]);names(RawCountR)[ncol(RawCountR)]<-gsub("\\.RawCounts.*","",i)}
RawCountR<-RawCountR[,-2]
setwd(paste0('/analysis2/xuy/RNA-seq/E250005786/6-Deseq/RepeatMasker'))
FilterCount<-RawCountR[apply(RawCountR[,-1],1,max)>5,]
CountMatrix<-as.matrix(FilterCount[,-1])
rownames(CountMatrix)<-FilterCount[,1]
colnames(CountMatrix)=c("AP2S1_1","AP2S1_2","AP2S1_3", "CMTR1_1", "CMTR1_2", "CMTR1_3", "DLD_1","DLD_2","DLD_3", "ELL_1", "ELL_2", "ELL_3", 
                        "FAM231A_1", "FAM231A_2", "NELFA_1", "NELFA_2", "NELFA_3", "NELFB_1", "NELFB_2", "NELFB_3", "NELFCD_1", "NELFCD_2", "NELFCD_3",
                        "NELFE_1", "NELFE_2", "NELFE_3", "control_1", "control_2", "control_3")
dds <- DESeqDataSetFromMatrix(countData = CountMatrix, colData = data.frame(genotype=c("AP2S1", "AP2S1", "AP2S1", "CMTR1", "CMTR1","CMTR1", "DLD", "DLD", "DLD", "ELL", "ELL", "ELL",
                                                                                       "FAM231A", "FAM231A", "NELFA", "NELFA", "NELFA", "NELFB", "NELFB", "NELFB", "NELFCD", "NELFCD", "NELFCD",
                                                                                       "NELFE", "NELFE", "NELFE", "control", "control", "control")), design = ~ genotype)
dds <- dds[rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)
#V<-as.data.frame(assay(varianceStabilizingTransformation(dds)))
vsd <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
# colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf("sampleDistMatrix.heatmap.pdf",height=11, width=11)
print(pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, cellwidth=20, cellheight=20))
dev.off()

Treat = c("AP2S1", "CMTR1", "DLD", "ELL", "FAM231A", "NELFA", "NELFB", "NELFCD", "NELFE")
for (i in Treat){
  ko=i
  wt="control"
  R <- results(dds, contrast=c("genotype",ko,wt))
  filename=paste0("RepeatMasker.Deseq2.", ko, ".result.txt")
  R2=R
  R2$TE=rownames(R)
  R2=cbind(R2$TE, R)
  colnames(R2)[1]="TE"
  write.table(R2, file=filename, sep="\t", quote=F, row.names=F)
  table (R$padj<.1)
  Rsig=as.data.frame(R[R$padj<.1 &!is.na(R$padj) , ])
  R<-as.data.frame(R)
  filenameSig=paste0("RepeatMasker.Deseq2.", ko, ".sig.padj.0.1.result.txt")
  write.table(Rsig, file=filenameSig, sep="\t", quote=F)
  #################### volcano plot####################
  # shrink 
  R <- results(dds, contrast=c("genotype",ko,wt))
  res1<- lfcShrink(dds,contrast = c("genotype",ko,wt), res=R, type = 'ashr')
  # use padj instead of p value
  res2 <- as.data.frame(res1) %>% rownames_to_column() %>% dplyr::select(-"pvalue")
  # res3 <- left_join(x = l1 ,y = res2,by="rowname")
  EnhancedVolcano(res2,
                  lab =res2$rowname,
                  x = 'log2FoldChange',
                  y = 'padj',
                  title = paste0(i," versus ",wt),
                  pCutoff = 0.05,
                  FCcutoff = 1,
                  pointSize = 3.0,
                  labSize = 3.0,
                  xlim = c(-5, 5),
                  ylim = c(0, 20),
                  legendPosition = "right",
                  ylab = bquote(~-Log[10]~italic(Padj)),
                  subtitle = "",
                  caption =  paste0("Total = ", nrow(res2), " genes")
  )
  ggsave(paste0(ko,"_vs_",wt,".Volcanoplot.pdf"),height = 20,width =30,units = "cm")
}

###### RepeatMasker Unique RawCount Data ######
setwd(paste0('/analysis2/xuy/RNA-seq/E250005786/5-RawCounts/RepeatMasker_Unique'))
TE=fread('/LiuLab/reference/Human/GRCh38/TE/Repeat_Masker_GRCh38_sorted.bed')
TE= TE %>% mutate(id=paste(V5,paste0("chr",V1),(V2+1),V3,V4, sep="_"))
FILE<-dir(pattern="*RawCounts*")
FILE<-FILE[grep ("RepeatMasker",FILE)]
RawCountR<-read.table(FILE[1], sep="\t", stringsAsFactors=FALSE, header=TRUE)
for (i in FILE) {RawCountR<-cbind(RawCountR, read.table(i, sep="\t", stringsAsFactors=FALSE, header=TRUE)[,2]);names(RawCountR)[ncol(RawCountR)]<-gsub("\\.RawCounts.*","",i)}
RawCountR<-RawCountR[,-2]
setwd(paste0('/analysis2/xuy/RNA-seq/E250005786/6-Deseq/RepeatMasker_Unique'))
FilterCount<-RawCountR[apply(RawCountR[,-1],1,max)>0,]
CountMatrix<-as.matrix(FilterCount[,-1])
rownames(CountMatrix)<-FilterCount[,1]
colnames(CountMatrix)=c("AP2S1_1","AP2S1_2","AP2S1_3", "CMTR1_1", "CMTR1_2", "CMTR1_3", "DLD_1","DLD_2","DLD_3", "ELL_1", "ELL_2", "ELL_3", 
                        "FAM231A_1", "FAM231A_2", "NELFA_1", "NELFA_2", "NELFA_3", "NELFB_1", "NELFB_2", "NELFB_3", "NELFCD_1", "NELFCD_2", "NELFCD_3",
                        "NELFE_1", "NELFE_2", "NELFE_3", "control_1", "control_2", "control_3")
dds <- DESeqDataSetFromMatrix(countData = CountMatrix, colData = data.frame(genotype=c("AP2S1", "AP2S1", "AP2S1", "CMTR1", "CMTR1","CMTR1", "DLD", "DLD", "DLD", "ELL", "ELL", "ELL",
                                                                                       "FAM231A", "FAM231A", "NELFA", "NELFA", "NELFA", "NELFB", "NELFB", "NELFB", "NELFCD", "NELFCD", "NELFCD",
                                                                                       "NELFE", "NELFE", "NELFE", "control", "control", "control")), design = ~ genotype)
dds <- dds[rowSums(counts(dds)) > 0, ]
dds <- DESeq(dds)
#V<-as.data.frame(assay(varianceStabilizingTransformation(dds)))
vsd <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
# colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf("sampleDistMatrix.heatmap.pdf",height=11, width=11)
print(pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, cellwidth=20, cellheight=20))
dev.off()

Treat = c("AP2S1", "CMTR1", "DLD", "ELL", "FAM231A", "NELFA", "NELFB", "NELFCD", "NELFE")
for (i in Treat){
  ko=i
  wt="control"
  R <- results(dds, contrast=c("genotype",ko,wt))
  filename=paste0("RepeatMasker.Deseq2.", ko, ".result.txt")
  R2=R
  R2$TE=rownames(R)
  R2=cbind(R2$TE, R)
  colnames(R2)[1]="TE"
  write.table(R2, file=filename, sep="\t", quote=F, row.names=F)
  table (R$padj<.1)
  Rsig=as.data.frame(R[R$padj<.1 &!is.na(R$padj) , ])
  R<-as.data.frame(R)
  filenameSig=paste0("RepeatMasker.Deseq2.", ko, ".sig.padj.0.1.result.txt")
  write.table(Rsig, file=filenameSig, sep="\t", quote=F)
  
  Rsig$id = rownames(Rsig)
  Ranno = Rsig %>% left_join(TE, by="id")
  write.table(Ranno, file=paste0("RepeatMasker.Deseq2.", ko, ".sig.padj.0.1.result.anno.txt"), sep="\t", quote=F,row.names=F)
  write.table(Ranno[which(Ranno$V6=="LTR" & Ranno$log2FoldChange > 0),], file=paste0("RepeatMasker.Deseq2.", ko, ".sig.padj.0.1.result.LTR.increase.anno.txt"), sep="\t", quote=F,row.names=F)
  write.table(Ranno[which(Ranno$V6=="LTR" & Ranno$log2FoldChange < 0),], file=paste0("RepeatMasker.Deseq2.", ko, ".sig.padj.0.1.result.LTR.decrease.anno.txt"), sep="\t", quote=F,row.names=F)
  write.table(Ranno[which(Ranno$V7=="L1" & Ranno$log2FoldChange > 0),], file=paste0("RepeatMasker.Deseq2.", ko, ".sig.padj.0.1.result.L1.increase.anno.txt"), sep="\t", quote=F,row.names=F)
  write.table(Ranno[which(Ranno$V7=="L1" & Ranno$log2FoldChange < 0),], file=paste0("RepeatMasker.Deseq2.", ko, ".sig.padj.0.1.result.L1.decrease.anno.txt"), sep="\t", quote=F,row.names=F)
  write.table(Ranno[which(Ranno$V7=="Alu" & Ranno$log2FoldChange > 0),], file=paste0("RepeatMasker.Deseq2.", ko, ".sig.padj.0.1.result.Alu.increase.anno.txt"), sep="\t", quote=F,row.names=F)
  write.table(Ranno[which(Ranno$V7=="Alu" & Ranno$log2FoldChange < 0),], file=paste0("RepeatMasker.Deseq2.", ko, ".sig.padj.0.1.result.Alu.decrease.anno.txt"), sep="\t", quote=F,row.names=F)
  write.table(Ranno[which(Ranno$V7=="SVA" & Ranno$log2FoldChange > 0),], file=paste0("RepeatMasker.Deseq2.", ko, ".sig.padj.0.1.result.SVA.increase.anno.txt"), sep="\t", quote=F,row.names=F)
  write.table(Ranno[which(Ranno$V7=="SVA" & Ranno$log2FoldChange < 0),], file=paste0("RepeatMasker.Deseq2.", ko, ".sig.padj.0.1.result.SVA.decrease.anno.txt"), sep="\t", quote=F,row.names=F)
  # volcano plot
  # shrink 
  R <- results(dds, contrast=c("genotype",ko,wt))
  res1<- lfcShrink(dds,contrast = c("genotype",ko,wt), res=R, type = 'ashr')
  # use padj instead of p value 
  res2 <- as.data.frame(res1) %>% rownames_to_column() %>% dplyr::select(-"pvalue")
  # res3 <- left_join(x = l1 ,y = res2,by="rowname")
  EnhancedVolcano(res2,
                  lab =res2$rowname,
                  x = 'log2FoldChange',
                  y = 'padj',
                  title = paste0(i," versus ",wt),
                  pCutoff = 0.05,
                  FCcutoff = 1,
                  pointSize = 3.0,
                  labSize = 3.0,
                  xlim = c(-5, 5),
                  ylim = c(0, 20),
                  legendPosition = "right",
                  ylab = bquote(~-Log[10]~italic(Padj)),
                  subtitle = "",
                  caption =  paste0("Total = ", nrow(res2), " genes")
  )
  ggsave(paste0(ko,"_vs_",wt,".Volcanoplot.pdf"),height = 20,width =30,units = "cm")
}

##### KnownGene RawCount Data ######
setwd(paste0('/analysis2/xuy/RNA-seq/E250005786/5-RawCounts/KnownGene'))
FILE<-dir(pattern="*RawCounts*")
FILE<-FILE[grep ("KnownGene",FILE)]
RawCountR<-read.table(FILE[1], sep="\t", stringsAsFactors=FALSE, header=TRUE)
for (i in FILE) {RawCountR<-cbind(RawCountR, read.table(i, sep="\t", stringsAsFactors=FALSE, header=TRUE)[,2]);names(RawCountR)[ncol(RawCountR)]<-gsub("\\.RawCounts.*","",i)}
RawCountR<-RawCountR[,-2]
setwd(paste0('/analysis2/xuy/RNA-seq/E250005786/6-Deseq/KnownGene'))
FilterCount<-RawCountR[apply(RawCountR[,-1],1,max)>20,]
CountMatrix<-as.matrix(FilterCount[,-1])
rownames(CountMatrix)<-FilterCount[,1]
colnames(CountMatrix)=c("AP2S1_1","AP2S1_2","AP2S1_3", "CMTR1_1", "CMTR1_2", "CMTR1_3", "DLD_1","DLD_2","DLD_3", "ELL_1", "ELL_2", "ELL_3", 
                        "FAM231A_1", "FAM231A_2", "NELFA_1", "NELFA_2", "NELFA_3", "NELFB_1", "NELFB_2", "NELFB_3", "NELFCD_1", "NELFCD_2", "NELFCD_3",
                        "NELFE_1", "NELFE_2", "NELFE_3", "control_1", "control_2", "control_3")
dds <- DESeqDataSetFromMatrix(countData = CountMatrix, colData = data.frame(genotype=c("AP2S1", "AP2S1", "AP2S1", "CMTR1", "CMTR1","CMTR1", "DLD", "DLD", "DLD", "ELL", "ELL", "ELL",
                                                                                       "FAM231A", "FAM231A", "NELFA", "NELFA", "NELFA", "NELFB", "NELFB", "NELFB", "NELFCD", "NELFCD", "NELFCD",
                                                                                       "NELFE", "NELFE", "NELFE", "control", "control", "control")), design = ~ genotype)
dds <- dds[rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)
#V<-as.data.frame(assay(varianceStabilizingTransformation(dds)))
vsd <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
#colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf("sampleDistMatrix.heatmap.pdf",height=11, width=11)
print(pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, cellwidth=20, cellheight=20))
dev.off()

Treat = c("AP2S1", "CMTR1", "DLD", "ELL", "FAM231A", "NELFA", "NELFB", "NELFCD", "NELFE")
for (i in Treat){
  ko=i
  wt="control"
  R <- results(dds, contrast=c("genotype",ko,wt))
  filename=paste("KnownGene.Deseq2.", ko, ".result.txt", sep="")
  R2=R
  R2$Gene=rownames(R)
  R2=cbind(R2$Gene, R)
  colnames(R2)[1]="Gene"
  write.table(R2, file=filename, sep="\t", quote=F, row.names=F)
  table (R$padj<.1)
  Rsig=as.data.frame(R[R$padj<.1 &!is.na(R$padj) , ])
  R<-as.data.frame(R)
  filenameSig=paste("KnownGene.Deseq2.", ko, ".sig.result.txt", sep="")
  write.table(Rsig, file=filenameSig, sep="\t", quote=F)
  # volcano plot
  # shrink 
  R <- results(dds, contrast=c("genotype",ko,wt))
  res1<- lfcShrink(dds,contrast = c("genotype",ko,wt), res=R, type = 'ashr')
  # use padj instead of p value 
  res2 <- as.data.frame(res1) %>% rownames_to_column() %>% dplyr::select(-"pvalue")
  # res3 <- left_join(x = l1 ,y = res2,by="rowname")
  EnhancedVolcano(res2,
                  lab =res2$rowname,
                  x = 'log2FoldChange',
                  y = 'padj',
                  title = paste0(i," versus ",wt),
                  pCutoff = 0.05,
                  FCcutoff = 1,
                  pointSize = 3.0,
                  labSize = 3.0,
                  xlim = c(-5, 5),
                  ylim = c(0, 20),
                  legendPosition = "right",
                  ylab = bquote(~-Log[10]~italic(Padj)),
                  subtitle = "",
                  caption =  paste0("Total = ", nrow(res2), " genes")
  )
  ggsave(paste0(ko,"_vs_",wt,".Volcanoplot.pdf"),height = 20,width =30,units = "cm")
}
