#!/pfs1/liuguanghui/liuzunpeng/04_Softwares/bin/Rscript

#Author: "Zunpeng Liu"
#Used in Liu Lab
#Email: zunpengliu@163.com
#date: "1/2/2018"
#Usages:
#Rscript Fragment_analysis.r [input] [out_prefix]

#------------------------------------------------------#
argv <- commandArgs(T)

if(length(argv) == 0){stop("
    Usages:
    Rscript Fragment_analysis.r [input] [out_prefix]
    input:       profile.
    out_prefix:  /pfs1/liuguanghui/liuzunpeng/05_Result/sample
    ")
}

dat<-read.csv(argv[1], header = TRUE)
out_prefix<- as.vector(argv[2])
#------------------------------------------------------#

library(ggplot2)
library(ggthemes)

pdf(paste(out_prefix,"_fragment_size_analysis.pdf",sep=""),width=8, height=6)

ggplot(dat,aes(dat$fragment_size))+
  geom_density(color="red")+
  xlim(c(0,1000))+
  labs(x="Fragment length(bp)",y="Normalized read density",title = "ATAC-seq fragment size analysis")+
  theme_classic()

dev.off()





