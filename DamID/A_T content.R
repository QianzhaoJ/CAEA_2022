library(ggpubr)

setwd("/Volumes/progeria/00_intergration/20_LAD_switch/nuc")

cont="H3"

treat<-c("R3","M3","W3","H9","R9","M9","W9")



for (i in treat) {
  tmp1<-read.csv(paste0(i,"_vs_",cont,"_stable_LAD.nuc.txt"),sep = "\t")
  tmp1$X4_usercol<-"stable_LAD"
  tmp2<-read.csv(paste0(i,"_vs_",cont,"_LAD_2_iLAD.nuc.txt"),sep = "\t")
  tmp2$X4_usercol<-"LAD_2_iLAD"
  tmp3<-read.csv(paste0(i,"_vs_",cont,"_iLAD_2_LAD.nuc.txt"),sep = "\t")
  tmp3$X4_usercol<-"iLAD_2_LAD"
  tmp4<-read.csv(paste0(i,"_vs_",cont,"_stable_iLAD.nuc.txt"),sep = "\t")
  tmp4$X4_usercol<-"stable_iLAD"
  
  tmp<-rbind(tmp1,tmp2,tmp3,tmp4)
  
  tmp$cluster<-factor(tmp$X4_usercol,levels = c("stable_LAD","LAD_2_iLAD","iLAD_2_LAD","stable_iLAD"))
  
  my_comparisons <- list(c("stable_LAD", "LAD_2_iLAD"), 
                         c("LAD_2_iLAD","iLAD_2_LAD"),
                         c("iLAD_2_LAD", "stable_iLAD"),
                         c("stable_LAD","iLAD_2_LAD"),
                         c("LAD_2_iLAD","stable_iLAD"),
                         c("stable_LAD", "stable_iLAD"))
  pdf(paste0(i,"_vs_",cont,".AT.percentage.pdf"),height =7 ,width = 5)
  p<-ggviolin(tmp, x="cluster", y="X7_pct_at",
              fill = "cluster", 
              color=NA, palette = c("#060608", "#a2a9af", "#dd7777","#dd0a35"),
              size=0,alpha = 0.9,remove = "point",add = "median_iqr",add.params = list(color="white",size=0.6,fill= "grey20"),
              legend = "right")+ 
    stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, label = "p.signif") +
    scale_y_continuous(limits = c(0,1.2),labels =seq(0,1.2,0.2),breaks =seq(0,1.2,0.2) )+
    labs(y="AT percentage")
  print(p)
  dev.off()
  
}

cont="H9"

treat<-c("R9","M9","W9")

for (i in treat) {
  tmp1<-read.csv(paste0(i,"_vs_",cont,"_stable_LAD.nuc.txt"),sep = "\t")
  tmp1$X4_usercol<-"stable_LAD"
  tmp2<-read.csv(paste0(i,"_vs_",cont,"_LAD_2_iLAD.nuc.txt"),sep = "\t")
  tmp2$X4_usercol<-"LAD_2_iLAD"
  tmp3<-read.csv(paste0(i,"_vs_",cont,"_iLAD_2_LAD.nuc.txt"),sep = "\t")
  tmp3$X4_usercol<-"iLAD_2_LAD"
  tmp4<-read.csv(paste0(i,"_vs_",cont,"_stable_iLAD.nuc.txt"),sep = "\t")
  tmp4$X4_usercol<-"stable_iLAD"
  
  tmp<-rbind(tmp1,tmp2,tmp3,tmp4)
  
  tmp$cluster<-factor(tmp$X4_usercol,levels = c("stable_LAD","LAD_2_iLAD","iLAD_2_LAD","stable_iLAD"))
  
  my_comparisons <- list(c("stable_LAD", "LAD_2_iLAD"), 
                         c("LAD_2_iLAD","iLAD_2_LAD"),
                         c("iLAD_2_LAD", "stable_iLAD"),
                         c("stable_LAD","iLAD_2_LAD"),
                         c("LAD_2_iLAD","stable_iLAD"),
                         c("stable_LAD", "stable_iLAD"))
  pdf(paste0(i,"_vs_",cont,".AT.percentage.pdf"),height =7 ,width = 5)
  p<-ggviolin(tmp, x="cluster", y="X7_pct_at",
              fill = "cluster", 
              color=NA, palette = c("#060608", "#a2a9af", "#dd7777","#dd0a35"),
              size=0,alpha = 0.9,remove = "point",add = "median_iqr",add.params = list(color="white",size=0.6,fill= "grey20"),
              legend = "right")+ 
    stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, label = "p.signif") +
    scale_y_continuous(limits = c(0,1.2),labels =seq(0,1.2,0.2),breaks =seq(0,1.2,0.2) )+
    labs(y="AT percentage")
  print(p)
  dev.off()
  
}
