setwd("/Volumes/progeria2/00_intergration/00_chromHMM_2021/S10/00_reorder/06_LAD/02_around")



hic<-read.csv("./WT_P3.txt",sep = "\t")
hic2<-hic[,-1]


colnames(hic2)


WT_P3_B.compartment.w25.bed

a<-c(paste0("WT_P3_LAD.compartment.w",seq(26,50,1),".bed"),
     paste0("WT_P3_A.compartment.w",seq(1,50,1),".bed"),
     paste0("WT_P3_B.compartment.w",seq(1,25,1),".bed"))



a<-c(paste0("WT_P3_iLAD.w",seq(26,50,1),".bed"),
     paste0("WT_P3_LAD.w",seq(1,50,1),".bed"),
     paste0("WT_P3_iLAD.w",seq(1,25,1),".bed"))

hic3<-hic2[-11,a]
hic3<-hic3[c(1,2,4,3,5,6,7,8,9,10),]



res3<-hic3

library(pheatmap)

colorl = c(  colorRampPalette(colors = c("white","white"))(80),
  colorRampPalette(colors = c("white","#3067D6"))(120),
  colorRampPalette(colors = c("#3067D6","#003399"))(100))
pheatmap(res3,
         scale = "none",
         show_rownames = T,
         cluster_rows  = F,
         cluster_cols  = F,
         color = colorl,
         border_color = NA ,
         cellwidth = 1,cellheight = 20 ,
         height=12, width=6  ,
         breaks = c(seq(0,3,0.01)),
         legend_breaks = seq(0,3,1),
         legend_labels = seq(0,3,1)   ,
         filename =  "LAD region final.2.pdf")

dev.off()
