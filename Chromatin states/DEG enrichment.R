library(pheatmap)
library(ggplot2)
library(ggthemes)
library(pheatmap)
library(pheatmap)

setwd("/Volumes/progeria2/00_intergration/00_chromHMM_2021/S10/00_reorder/29_DEGs/")


hic<-read.csv("./WT_P3.txt",sep = "\t")
hic2<-hic[,-1]


colnames(hic2)



a<-c("RS_up.bed",
     "HGPS_Hetero_up.bed",
     "HGPS_Homo_up.bed",
     "WS_up.bed",
     
     "RS_down.bed",
     "HGPS_Hetero_down.bed",
     "HGPS_Homo_down.bed",
     "WS_down.bed")

hic3<-hic2[-11,a]

hic3<-hic3[c(1,2,4,3,5,6,7,8,9,10),]


res2<-hic3



colorl = c(colorRampPalette(colors = c("white","white"))(20),
           colorRampPalette(colors = c("white","brown2"))(120),
           colorRampPalette(colors = c("brown2","darkred"))(160))

pheatmap(res2,
         scale = "row",
         show_rownames = T,
         cluster_rows  = F,
         cluster_cols  = F,
         color = colorl,
         border_color = "white" ,
         cellwidth = 20,cellheight = 20 ,  gaps_col = 4,
         height=12, width=20 ,
         breaks = c(seq(0,2,0.01)),
         legend_breaks = seq(0,2,1),
         legend_labels = seq(0,2,1)   ,
         filename =  "DEG enrichment.pdf")


