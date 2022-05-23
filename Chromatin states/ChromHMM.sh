#!/bin/sh
java -mx1000000M -jar \
/work2/liugh/liuzunpeng/04_Softwares/ChromHMM/ChromHMM/ChromHMM.jar \
BinarizeBam \
/work2/liugh/liuzunpeng/04_Softwares/ChromHMM/ChromHMM/CHROMSIZES/hg19.sim.txt \
/work2/liugh/liuzunpeng/05_Results/72_ChromHMM_progeria_new/01_bam \
/work2/liugh/liuzunpeng/05_Results/72_ChromHMM_progeria_new/00_src/sam.tab \
/work2/liugh/liuzunpeng/05_Results/72_ChromHMM_progeria_new/02_binarizebam


java -mx1000000M -jar /work2/liugh/liuzunpeng/04_Softwares/ChromHMM/ChromHMM/ChromHMM.jar LearnModel -p 24 \
/work2/liugh/liuzunpeng/05_Results/72_ChromHMM_progeria_new/02_binarizebam \
/work2/liugh/liuzunpeng/05_Results/72_ChromHMM_progeria_new/03_LearnModel_200bp/S10 10 hg19


java -mx4000M -jar /work2/liugh/liuzunpeng/04_Softwares/ChromHMM/ChromHMM/ChromHMM.jar Reorder \
-o /work2/liugh/liuzunpeng/05_Results/72_ChromHMM_progeria_new/03_LearnModel_200bp/S10/reorder/stateorderingfile.txt \
/work2/liugh/liuzunpeng/05_Results/72_ChromHMM_progeria_new/03_LearnModel_200bp/S10/model_10.txt  \
/work2/liugh/liuzunpeng/05_Results/72_ChromHMM_progeria_new/03_LearnModel_200bp/S10/reorder


#!/bin/sh

out="/work2/liugh/liuzunpeng/05_Results/72_ChromHMM_progeria_new/03_LearnModel_200bp/S10/reorder/01_Segmentation/split"

for i in Hetero_P3 Homo_P3 WS_P3 WT_P3 Hetero_P9 Homo_P9 WS_P9 WT_P9

do

file="/work2/liugh/liuzunpeng/05_Results/72_ChromHMM_progeria_new/03_LearnModel_200bp/S10/reorder/01_Segmentation/"$i"_10_segments.bed"

grep U1 $file |awk 'OFS="\t" {print $1,$2,$3,$4,".\t+"}' > ${file%%.}"_U1.bed"
grep U2 $file |awk 'OFS="\t" {print $1,$2,$3,$4,".\t+"}' > ${file%%.}"_U2.bed"
grep U3 $file |awk 'OFS="\t" {print $1,$2,$3,$4,".\t+"}' > ${file%%.}"_U3.bed"
grep U4 $file |awk 'OFS="\t" {print $1,$2,$3,$4,".\t+"}' > ${file%%.}"_U4.bed"
grep U5 $file |awk 'OFS="\t" {print $1,$2,$3,$4,".\t+"}' > ${file%%.}"_U5.bed"
grep U6 $file |awk 'OFS="\t" {print $1,$2,$3,$4,".\t+"}' > ${file%%.}"_U6.bed"
grep U7 $file |awk 'OFS="\t" {print $1,$2,$3,$4,".\t+"}' > ${file%%.}"_U7.bed"
grep U8 $file |awk 'OFS="\t" {print $1,$2,$3,$4,".\t+"}' > ${file%%.}"_U8.bed"
grep U9 $file |awk 'OFS="\t" {print $1,$2,$3,$4,".\t+"}' > ${file%%.}"_U9.bed"
grep U10 $file |awk 'OFS="\t" {print $1,$2,$3,$4,".\t+"}' > ${file%%.}"_U10.bed"



done

