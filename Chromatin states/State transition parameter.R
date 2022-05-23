setwd("/Volumes/progeria2/00_intergration/00_chromHMM_2021/S10/00_reorder")


library(pheatmap)
library(ggplot2)
library(ggthemes)

clus<-read.csv("./transitions_10.txt",sep = "\t")
clus2<-clus[,-1]
clus2<-clus2[c(1,2,4,3,5,6,7,8,9,10),c(1,2,4,3,5,6,7,8,9,10)]

colorl = c(colorRampPalette(colors = c("white","#1F5CE5"))(10),
           colorRampPalette(colors = c("#1F5CE5","#0033A5"))(50))


pheatmap(clus2,
         scale = "none",
         show_rownames = T,
         cluster_rows = F,
         cluster_cols = F ,
         color = colorl,
         border_color = "grey90" ,
         cellwidth = 10,
         cellheight = 10  ,
         #annotation_row = anno,
         #labels_row= paste0("cluster_",seq(1,22,1)),
         height=8, width=6 ,
         filename =  "transitions_10.finalize.pdf")

dev.off()




