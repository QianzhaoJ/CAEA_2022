#!/pfs1/liuguanghui/liuzunpeng/04_Softwares/bin/Rscript

#----------------------------------------------------------------------------------#
#  Title: "Annotation ChIP/ATAC-seq peak pipeline"                                 #
#  Author: "Zunpeng Liu"                                                           #
#  Email: zunpengliu\@163.com                                                      #
#  Date: "1/6/2018"                                                                #
#  Usages:                                                                         #
#  Rscript ChIPseeker_peak_anno.r [input file] [out_direction] [sample name]       #
#  Example:                                                                        #
#  /pfs1/liuguanghui/liuzunpeng/04_Softwares/bin/Rscript ChIPseeker_peak_anno.r \  #
#  [input]                                                                         #
#  [out_prefix]                                                                    #
#                                                                                  #
#  Only used for Liu Lab !                                                         #
#----------------------------------------------------------------------------------#

argv <- commandArgs(T)

if(length(argv)!= 3){stop("
    Annotation ATAC-seq peak by ChIPseeker
    Usages:
    Rscript ChIPseeker_peak_anno.r [input] [out] [name]
    input:       profile
    out_prefix:  /pfs1/liuguanghui/liuzunpeng/05_Result/sample
    ")
}

dat<-argv[1]
out<-argv[2]
sample<-argv[3]
prefix<- paste(out,"/",sample,sep="")

setwd(out)
#------------------------------------------------------#
## loading all packages needed

library(GenomicFeatures)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(ChIPseeker)
library(ggplot2)
library(ggthemes)

#------------------------------------------------------#

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
chr<-c(paste("chr",seq(22),sep = "") ,"chrX")

bedPeaksFile = file.path(dat)
peak<-readPeakFile(bedPeaksFile)
tagMatrix <- getTagMatrix(peak, windows=promoter)

peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnno_df <- as.data.frame(peakAnno)
write.csv(peakAnno_df,paste(sample,'_peakAnno_df.csv',sep=""))

pdf('p1_ChIP Peaks over Chromosomes.pdf',width=8,height=12)
covplot(peak, chrs = c(paste("chr",seq(22),sep = ""),"chrX"),weightCol = "V5")
dev.off()

png('p2_TagHeatmap of peaks binding to TSS region.tiff',width=500,height=1200)
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red",xlab="Enrichment at TSS region")
dev.off()

pdf('p2_TagHeatmap of peaks binding to TSS region.pdf',width=5,height=12)
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red",xlab="Enrichment at TSS region")
dev.off()

pdf('p3_Average Profile of peaks binding to TSS region.pdf',width=10,height=10)
plotAvgProf(tagMatrix, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.off()

pdf('p4_Pie-summarize the distribution of peaks.pdf',width=8,height=8)
plotAnnoPie(peakAnno)
dev.off()

pdf('p5_Bar-summarize the distribution of peaks.pdf',width=10,height=3.6)
plotAnnoBar(peakAnno)
dev.off()

pdf('p6_vennpie-summarize the distribution of peaks.pdf',width=10,height=7)
vennpie(peakAnno)
dev.off()

pdf('p7_upsetplot-summarize the distribution of peaks.pdf',width=12,height=8)
upsetplot(peakAnno)
dev.off()

pdf('p8_upsetplot_vennpie-summarize the distribution of peaks.pdf',width=12,height=8)
upsetplot(peakAnno, vennpie=TRUE)
dev.off()

pdf('p9_Distribution of transcription factor-binding loci relative to TSS.pdf',width=8,height=2.6)
plotDistToTSS(peakAnno, title="Distribution of peaks relative to TSS")
dev.off()








