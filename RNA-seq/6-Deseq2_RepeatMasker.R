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
setwd(paste0('/analysis2/xuy/RNA-seq/RNA-seq-THG15-panotaz-24-48-72/5-RawCounts/RepeatMasker'))
FILE<-dir(pattern="*RawCounts*")
FILE<-FILE[grep("RepeatMasker",FILE, value = F)]
RawCountR<-read.table(FILE[1], sep="\t", stringsAsFactors=FALSE, header=TRUE)
for (i in FILE) {RawCountR<-cbind(RawCountR, read.table(i, sep="\t", stringsAsFactors=FALSE, header=TRUE)[,2]);names(RawCountR)[ncol(RawCountR)]<-gsub("\\.RawCounts.*","",i)}
RawCountR<-RawCountR[,-2]
setwd(paste0('/analysis2/xuy/RNA-seq/RNA-seq-THG15-panotaz-24-48-72/6-Deseq/RepeatMasker'))
FilterCount<-RawCountR[apply(RawCountR[,-1],1,max)>5,]
CountMatrix<-as.matrix(FilterCount[,-1])
rownames(CountMatrix)<-FilterCount[,1]
colnames(CountMatrix)=c("DMSO_24_1","DMSO_24_2","DMSO_48_1", "DMSO_48_2", "DMSO_72_1", "DMSO_72_2", "PanTaz_24_1","PanTaz_24_2","PanTaz_48_1", "PanTaz_48_2", "PanTaz_72_1", "PanTaz_72_2")
dds <- DESeqDataSetFromMatrix(countData = CountMatrix, colData = data.frame(genotype=c("DMSO_24","DMSO_24","DMSO_48","DMSO_48", "DMSO_72","DMSO_72", "PanTaz_24","PanTaz_24","PanTaz_48","PanTaz_48", "PanTaz_72","PanTaz_72")), design = ~ genotype)
dds <- dds[rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)
V<-as.data.frame(assay(varianceStabilizingTransformation(dds)))
vsd <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
# colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf("sampleDistMatrix.heatmap.pdf",height=6, width=6)
	print(pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, cellwidth=20, cellheight=20))
dev.off()
#################################### stop here 
Treat = c("PANO", "PANOTAZE", "TAZE")
for (i in Treat){
		ko=i
		wt="DMSO"
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
}

Treat = c("PANOTAZE")
for (i in Treat){
		ko=i
		wt="PANO"
		R <- results(dds, contrast=c("genotype",ko,wt))
		filename=paste0("RepeatMasker.Deseq2.", ko, "_", wt, ".result.txt")
		R2=R
		R2$TE=rownames(R)
		R2=cbind(R2$TE, R)
		colnames(R2)[1]="TE"
		write.table(R2, file=filename, sep="\t", quote=F, row.names=F)
		table (R$padj<.1)
		Rsig=as.data.frame(R[R$padj<.1 &!is.na(R$padj) , ])
		R<-as.data.frame(R)
		filenameSig=paste0("RepeatMasker.Deseq2.", ko, "_", wt,".sig.padj.0.1.result.txt")
		write.table(Rsig, file=filenameSig, sep="\t", quote=F)
}

Treat = c("PANOTAZE", "PANO")
for (i in Treat){
		ko=i
		wt="TAZE"
		R <- results(dds, contrast=c("genotype",ko,wt))
		filename=paste0("RepeatMasker.Deseq2.", ko, "_", wt, ".result.txt")
		R2=R
		R2$TE=rownames(R)
		R2=cbind(R2$TE, R)
		colnames(R2)[1]="TE"
		write.table(R2, file=filename, sep="\t", quote=F, row.names=F)
		table (R$padj<.1)
		Rsig=as.data.frame(R[R$padj<.1 &!is.na(R$padj) , ])
		R<-as.data.frame(R)
		filenameSig=paste0("RepeatMasker.Deseq2.", ko, "_", wt,".sig.padj.0.1.result.txt")
		write.table(Rsig, file=filenameSig, sep="\t", quote=F)
}


###### RepeatMasker Unique RawCount Data ######
setwd(paste0('/analysis2/xuy/RNA-seq/RNA-seq-THG15-panotaz-24-48-72/5-RawCounts/RepeatMasker_Unique'))
TE=fread('/LiuLab/reference/Human/GRCh38/TE/Repeat_Masker_GRCh38_sorted.bed')
TE= TE %>% mutate(id=paste(V5,paste0("chr",V1),(V2+1),V3,V4, sep="_"))
FILE<-dir(pattern="*RawCounts*")
FILE<-FILE[grep ("RepeatMasker",FILE)]
RawCountR<-read.table(FILE[1], sep="\t", stringsAsFactors=FALSE, header=TRUE)
for (i in FILE) {RawCountR<-cbind(RawCountR, read.table(i, sep="\t", stringsAsFactors=FALSE, header=TRUE)[,2]);names(RawCountR)[ncol(RawCountR)]<-gsub("\\.RawCounts.*","",i)}
RawCountR<-RawCountR[,-2]
setwd(paste0('/analysis2/xuy/RNA-seq/RNA-seq-THG15-panotaz-24-48-72/6-Deseq/RepeatMasker_Unique'))
FilterCount<-RawCountR[apply(RawCountR[,-1],1,max)>0,]
CountMatrix<-as.matrix(FilterCount[,-1])
rownames(CountMatrix)<-FilterCount[,1]
colnames(CountMatrix)=c("DMSO_1","DMSO_2","PANO_1", "PANO_2", "PANOTAZE_1", "PANOTAZE_2", "TAZE_1", "TAZE_2")
dds <- DESeqDataSetFromMatrix(countData = CountMatrix, colData = data.frame(genotype=c("DMSO","DMSO","PANO", "PANO", "PANOTAZE", "PANOTAZE", "TAZE", "TAZE")), design = ~ genotype)
dds <- dds[rowSums(counts(dds)) > 0, ]
dds <- DESeq(dds)
V<-as.data.frame(assay(varianceStabilizingTransformation(dds)))
vsd <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
# colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf("sampleDistMatrix.heatmap.pdf",height=7, width=7)
	print(pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, cellwidth=20, cellheight=20))
dev.off()

Treat = c("PANO", "PANOTAZE", "TAZE")
for (i in Treat){
	ko=i
	wt="DMSO"
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
}



Treat = c("PANOTAZE")
for (i in Treat){
	ko=i
	wt="PANO"
	R <- results(dds, contrast=c("genotype",ko,wt))
	filename=paste0("RepeatMasker.Deseq2.", ko, "_", wt, ".result.txt")
	R2=R
	R2$TE=rownames(R)
	R2=cbind(R2$TE, R)
	colnames(R2)[1]="TE"
	write.table(R2, file=filename, sep="\t", quote=F, row.names=F)
	table (R$padj<.1)
	Rsig=as.data.frame(R[R$padj<.1 &!is.na(R$padj) , ])
	R<-as.data.frame(R)
	filenameSig=paste0("RepeatMasker.Deseq2.", ko, "_", wt, ".sig.padj.0.1.result.txt")
	write.table(Rsig, file=filenameSig, sep="\t", quote=F)

	Rsig$id = rownames(Rsig)
	Ranno = Rsig %>% left_join(TE, by="id")
	write.table(Ranno, file=paste0("RepeatMasker.Deseq2.", ko, "_", wt, ".sig.padj.0.1.result.anno.txt"), sep="\t", quote=F,row.names=F)
	write.table(Ranno[which(Ranno$V6=="LTR" & Ranno$log2FoldChange > 0),], file=paste0("RepeatMasker.Deseq2.", ko, "_", wt, ".sig.padj.0.1.result.LTR.increase.anno.txt"), sep="\t", quote=F,row.names=F)
	write.table(Ranno[which(Ranno$V6=="LTR" & Ranno$log2FoldChange < 0),], file=paste0("RepeatMasker.Deseq2.", ko, "_", wt, ".sig.padj.0.1.result.LTR.decrease.anno.txt"), sep="\t", quote=F,row.names=F)
	write.table(Ranno[which(Ranno$V7=="L1" & Ranno$log2FoldChange > 0),], file=paste0("RepeatMasker.Deseq2.", ko, "_", wt, ".sig.padj.0.1.result.L1.increase.anno.txt"), sep="\t", quote=F,row.names=F)
	write.table(Ranno[which(Ranno$V7=="L1" & Ranno$log2FoldChange < 0),], file=paste0("RepeatMasker.Deseq2.", ko, "_", wt, ".sig.padj.0.1.result.L1.decrease.anno.txt"), sep="\t", quote=F,row.names=F)
	write.table(Ranno[which(Ranno$V7=="Alu" & Ranno$log2FoldChange > 0),], file=paste0("RepeatMasker.Deseq2.", ko, "_", wt, ".sig.padj.0.1.result.Alu.increase.anno.txt"), sep="\t", quote=F,row.names=F)
	write.table(Ranno[which(Ranno$V7=="Alu" & Ranno$log2FoldChange < 0),], file=paste0("RepeatMasker.Deseq2.", ko, "_", wt, ".sig.padj.0.1.result.Alu.decrease.anno.txt"), sep="\t", quote=F,row.names=F)
	write.table(Ranno[which(Ranno$V7=="SVA" & Ranno$log2FoldChange > 0),], file=paste0("RepeatMasker.Deseq2.", ko, "_", wt, ".sig.padj.0.1.result.SVA.increase.anno.txt"), sep="\t", quote=F,row.names=F)
	write.table(Ranno[which(Ranno$V7=="SVA" & Ranno$log2FoldChange < 0),], file=paste0("RepeatMasker.Deseq2.", ko, "_", wt, ".sig.padj.0.1.result.SVA.decrease.anno.txt"), sep="\t", quote=F,row.names=F)
}


Treat = c("PANOTAZE", "PANO")
for (i in Treat){
	ko=i
	wt="TAZE"
	R <- results(dds, contrast=c("genotype",ko,wt))
	filename=paste0("RepeatMasker.Deseq2.", ko, "_", wt, ".result.txt")
	R2=R
	R2$TE=rownames(R)
	R2=cbind(R2$TE, R)
	colnames(R2)[1]="TE"
	write.table(R2, file=filename, sep="\t", quote=F, row.names=F)
	table (R$padj<.1)
	Rsig=as.data.frame(R[R$padj<.1 &!is.na(R$padj) , ])
	R<-as.data.frame(R)
	filenameSig=paste0("RepeatMasker.Deseq2.", ko, "_", wt, ".sig.padj.0.1.result.txt")
	write.table(Rsig, file=filenameSig, sep="\t", quote=F)

	Rsig$id = rownames(Rsig)
	Ranno = Rsig %>% left_join(TE, by="id")
	write.table(Ranno, file=paste0("RepeatMasker.Deseq2.", ko, "_", wt, ".sig.padj.0.1.result.anno.txt"), sep="\t", quote=F,row.names=F)
	write.table(Ranno[which(Ranno$V6=="LTR" & Ranno$log2FoldChange > 0),], file=paste0("RepeatMasker.Deseq2.", ko, "_", wt, ".sig.padj.0.1.result.LTR.increase.anno.txt"), sep="\t", quote=F,row.names=F)
	write.table(Ranno[which(Ranno$V6=="LTR" & Ranno$log2FoldChange < 0),], file=paste0("RepeatMasker.Deseq2.", ko, "_", wt, ".sig.padj.0.1.result.LTR.decrease.anno.txt"), sep="\t", quote=F,row.names=F)
	write.table(Ranno[which(Ranno$V7=="L1" & Ranno$log2FoldChange > 0),], file=paste0("RepeatMasker.Deseq2.", ko, "_", wt, ".sig.padj.0.1.result.L1.increase.anno.txt"), sep="\t", quote=F,row.names=F)
	write.table(Ranno[which(Ranno$V7=="L1" & Ranno$log2FoldChange < 0),], file=paste0("RepeatMasker.Deseq2.", ko, "_", wt, ".sig.padj.0.1.result.L1.decrease.anno.txt"), sep="\t", quote=F,row.names=F)
	write.table(Ranno[which(Ranno$V7=="Alu" & Ranno$log2FoldChange > 0),], file=paste0("RepeatMasker.Deseq2.", ko, "_", wt, ".sig.padj.0.1.result.Alu.increase.anno.txt"), sep="\t", quote=F,row.names=F)
	write.table(Ranno[which(Ranno$V7=="Alu" & Ranno$log2FoldChange < 0),], file=paste0("RepeatMasker.Deseq2.", ko, "_", wt, ".sig.padj.0.1.result.Alu.decrease.anno.txt"), sep="\t", quote=F,row.names=F)
	write.table(Ranno[which(Ranno$V7=="SVA" & Ranno$log2FoldChange > 0),], file=paste0("RepeatMasker.Deseq2.", ko, "_", wt, ".sig.padj.0.1.result.SVA.increase.anno.txt"), sep="\t", quote=F,row.names=F)
	write.table(Ranno[which(Ranno$V7=="SVA" & Ranno$log2FoldChange < 0),], file=paste0("RepeatMasker.Deseq2.", ko, "_", wt, ".sig.padj.0.1.result.SVA.decrease.anno.txt"), sep="\t", quote=F,row.names=F)
}


##### KnownGene RawCount Data ######
setwd(paste0('/analysis2/xuy/RNA-seq/RNA-seq-THG15-panotaz-24-48-72/5-RawCounts/KnownGene'))
FILE<-dir(pattern="*RawCounts*")
FILE<-FILE[grep ("KnownGene",FILE)]
RawCountR<-read.table(FILE[1], sep="\t", stringsAsFactors=FALSE, header=TRUE)
for (i in FILE) {RawCountR<-cbind(RawCountR, read.table(i, sep="\t", stringsAsFactors=FALSE, header=TRUE)[,2]);names(RawCountR)[ncol(RawCountR)]<-gsub("\\.RawCounts.*","",i)}
RawCountR<-RawCountR[,-2]
setwd(paste0('/analysis2/xuy/RNA-seq/RNA-seq-THG15-panotaz-24-48-72/6-Deseq/KnownGene'))
FilterCount<-RawCountR[apply(RawCountR[,-1],1,max)>20,]
CountMatrix<-as.matrix(FilterCount[,-1])
rownames(CountMatrix)<-FilterCount[,1]
colnames(CountMatrix)=c("DMSO_1","DMSO_2","PANO_1", "PANO_2", "PANOTAZE_1", "PANOTAZE_2", "TAZE_1", "TAZE_2")
dds <- DESeqDataSetFromMatrix(countData = CountMatrix, colData = data.frame(genotype=c("DMSO","DMSO","PANO", "PANO", "PANOTAZE", "PANOTAZE", "TAZE", "TAZE")), design = ~ genotype)
dds <- dds[rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)
V<-as.data.frame(assay(varianceStabilizingTransformation(dds)))
vsd <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf("sampleDistMatrix.heatmap.pdf",height=6, width=7)
	print(pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, cellwidth=20, cellheight=20))
dev.off()

Treat = c("PANO", "PANOTAZE", "TAZE")
for (i in Treat){
	ko=i
	wt="DMSO"
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
   
   	#  note , here pCutoff actually means padj Cutoff
    EnhancedVolcano(res2,
                        lab =res2$rowname,
                        x = 'log2FoldChange',
                        y = 'padj',
                        title = paste0(i," versus WT"),
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


Treat = c("PANOTAZE")
for (i in Treat){
	ko=i
	wt="PANO"
	R <- results(dds, contrast=c("genotype",ko,wt))
	filename=paste("KnownGene.Deseq2.", ko, "_", wt, ".result.txt", sep="")
	R2=R
	R2$Gene=rownames(R)
	R2=cbind(R2$Gene, R)
	colnames(R2)[1]="Gene"
	write.table(R2, file=filename, sep="\t", quote=F, row.names=F)
	table (R$padj<.1)
	Rsig=as.data.frame(R[R$padj<.1 &!is.na(R$padj) , ])
	R<-as.data.frame(R)
	filenameSig=paste("KnownGene.Deseq2.", ko, "_", wt, ".sig.result.txt", sep="")
	write.table(Rsig, file=filenameSig, sep="\t", quote=F)
	

	# volcano plot
   	# shrink 
	R <- results(dds, contrast=c("genotype",ko,wt))
   	res1<- lfcShrink(dds,contrast = c("genotype",ko,wt), res=R, type = 'ashr')
  	# use padj instead of p value 
   	res2 <- as.data.frame(res1) %>% rownames_to_column() %>% dplyr::select(-"pvalue")

   	# res3 <- left_join(x = l1 ,y = res2,by="rowname")
   
   	#  note , here pCutoff actually means padj Cutoff
    EnhancedVolcano(res2,
                        lab =res2$rowname,
                        x = 'log2FoldChange',
                        y = 'padj',
                        title = paste0(i," versus WT"),
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

Treat = c("PANOTAZE", "PANO")
for (i in Treat){
	ko=i
	wt="TAZE"
	R <- results(dds, contrast=c("genotype",ko,wt))
	filename=paste("KnownGene.Deseq2.", ko, "_", wt, ".result.txt", sep="")
	R2=R
	R2$Gene=rownames(R)
	R2=cbind(R2$Gene, R)
	colnames(R2)[1]="Gene"
	write.table(R2, file=filename, sep="\t", quote=F, row.names=F)
	table (R$padj<.1)
	Rsig=as.data.frame(R[R$padj<.1 &!is.na(R$padj) , ])
	R<-as.data.frame(R)
	filenameSig=paste("KnownGene.Deseq2.", ko, "_", wt, ".sig.result.txt", sep="")
	write.table(Rsig, file=filenameSig, sep="\t", quote=F)
	# volcano plot
   	# shrink 
	R <- results(dds, contrast=c("genotype",ko,wt))
   	res1<- lfcShrink(dds,contrast = c("genotype",ko,wt), res=R, type = 'ashr')
  	# use padj instead of p value 
   	res2 <- as.data.frame(res1) %>% rownames_to_column() %>% dplyr::select(-"pvalue")

   	# res3 <- left_join(x = l1 ,y = res2,by="rowname")
   
   	#  note , here pCutoff actually means padj Cutoff
    EnhancedVolcano(res2,
                        lab =res2$rowname,
                        x = 'log2FoldChange',
                        y = 'padj',
                        title = paste0(i," versus WT"),
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

