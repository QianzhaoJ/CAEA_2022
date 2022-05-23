setwd("/Volumes/progeria2/00_intergration/00_chromHMM_2021/S10/00_reorder")


library(pheatmap)
library(ggplot2)
library(ggthemes)

clus<-read.csv("./emissions_10.txt",sep = "\t")
clus2<-clus[,-1]


clus2$cluster<-paste0("cluster_",seq(1,10,1))
colnames(clus2)
clus3<-clus2[,c("H3K9me3","H4K20me3","H3K27me3","ATAC","H3K27ac","H3K4me1","H3K4me3","H3K36me3","cluster")]

clus3$cluster2<-factor(clus3$cluster,levels =  c(paste0("cluster_",seq(1,10,1))))


clus3<-clus3[order(clus3$cluster2,decreasing = F),]
clus4<-clus3[,-c(9:10)]
rownames(clus4)<-clus3$cluster2

anno<-data.frame(old=c( paste0("cluster_",seq(1,10,1)),
                        new=factor(paste0("cluster_",seq(1,10,1)),
                                   levels = paste0("cluster_",seq(1,10,1)))))

anno<-data.frame(anno[,-1])
rownames(anno)<-clus3$cluster2


colorl = c(colorRampPalette(colors = c("white","white"))(80),
           colorRampPalette(colors = c("white","dodgerblue2"))(110),
           colorRampPalette(colors = c("dodgerblue2","blue"))(100),
           colorRampPalette(colors = c("blue","navy"))(280))


pheatmap(clus4,
         scale = "none",
         show_rownames = T,
         cluster_rows = F,cluster_cols = F,
         color = colorl,
         border_color = "grey80" ,
         cellwidth = 20,cellheight = 20,
         #annotation_row = anno,
         #labels_row= paste0("cluster_",seq(1,22,1)),
         height=8, width=6  ,
         filename =  "all.cluster_new2.pdf")

dev.off()




