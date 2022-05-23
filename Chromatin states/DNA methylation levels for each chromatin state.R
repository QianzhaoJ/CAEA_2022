setwd("/Volumes/progeria2/00_intergration/00_chromHMM_2021/S10/00_reorder/03_DNA_methylation/")
library(data.table)
library(ggplot2)
library(ggthemes)


H3<-fread("./WT_P3.segments.CpG_methylation.bed",sep = "\t",header = F)
R3<-fread("./Hetero_P3.segments.CpG_methylation.bed",sep = "\t",header = F)
M3<-fread("./Homo_P3.segments.CpG_methylation.bed",sep = "\t",header = F)
W3<-fread("./WS_P3.segments.CpG_methylation.bed",sep = "\t",header = F)

H9<-fread("./WT_P9.segments.CpG_methylation.bed",sep = "\t",header = F)
R9<-fread("./Hetero_P9.segments.CpG_methylation.bed",sep = "\t",header = F)
M9<-fread("./Homo_P9.segments.CpG_methylation.bed",sep = "\t",header = F)
W9<-fread("./WS_P9.segments.CpG_methylation.bed",sep = "\t",header = F)


dna_me<-rbind(H3,R3,M3,W3,
              H9,R9,M9,W9)


dna_me$samp<-c(rep("H3",nrow(H3)),
               rep("R3",nrow(R3)),
               rep("M3",nrow(M3)),
               rep("W3",nrow(W3)),
               rep("H9",nrow(H9)),
               rep("R9",nrow(R9)),
               rep("M9",nrow(M9)),
               rep("W9",nrow(W9))) 


dna_me$samp2<-factor(dna_me$samp,levels = c("H3","R3","M3","W3","H9","R9","M9","W9"))

dna_me$stat2<-factor(dna_me$V4,levels=c(paste0("U",seq(1,10))))



pdf("DNA_methy_for_all_sample_in_each_Cluster.pdf",width = 15,height = 6)
ggplot(dna_me,aes(samp2, V5,fill=samp2)) +
  geom_boxplot(outlier.shape = NA,width=0.6)+
  #geom_text(stat = "stratum", label.strata = TRUE) +
  scale_fill_manual(values=c("navyblue","dodgerblue2","deepskyblue2",
                             "lightblue3","magenta4","red2","hotpink","darksalmon"))+
  theme_classic() +
  facet_grid(.~stat2)+
  ylim(c(0,1))
dev.off()




dna_me2<-subset(dna_me,dna_me$samp2=="H3")

pdf("DNA_methy_for_all_sample_in_each_Cluster in WT.pdf",width = 9,height = 2.5)
ggplot(dna_me2,aes(stat2, V5,fill=stat2)) +
  geom_boxplot(outlier.shape = NA,width=0.6)+
  #geom_text(stat = "stratum", label.strata = TRUE) +
  scale_fill_manual(values=c("#003B99",
                             "#0046CB",
                             "#66CCFF",
                             "#0099FF",
                             "#00AABE",
                             
                             "#D02090",
                             "#F97DC1",
                             "#0C8140",
                             "#EE2025",
                             "#CCCCCC"))+
  theme_few() +
  ylim(c(0,1))
dev.off()




