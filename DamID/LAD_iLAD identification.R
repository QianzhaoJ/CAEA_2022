#!/bin/sh
java -mx100000M -jar /work2/liugh/liuzunpeng/04_Softwares/ChromHMM/ChromHMM/ChromHMM.jar BinarizeBam -b 100000 /work2/liugh/liuzunpeng/04_Softwares/ChromHMM/ChromHMM/CHROMSIZES/hg19.sim.txt /work2/liugh/liuzunpeng/05_Results/59_Lamin_KI_DamID_part2/02_Bam /work2/liugh/liuzunpeng/05_Results/59_Lamin_KI_DamID_part2/03_binarizebam/dam.tab /work2/liugh/liuzunpeng/05_Results/59_Lamin_KI_DamID_part2/03_binarizebam

java -mx100000M -jar /work2/liugh/liuzunpeng/04_Softwares/ChromHMM/ChromHMM/ChromHMM.jar BinarizeBam -b 100000 /work2/liugh/liuzunpeng/04_Softwares/ChromHMM/ChromHMM/CHROMSIZES/hg19.sim.txt /work2/liugh/liuzunpeng/05_Results/59_Lamin_KI_DamID_part2/02_Bam /work2/liugh/liuzunpeng/05_Results/59_Lamin_KI_DamID_part2/03_binarizebam/dam.tab /work2/liugh/liuzunpeng/05_Results/59_Lamin_KI_DamID_part2/03_binarizebam

java -mx40000M -jar /work2/liugh/liuzunpeng/04_Softwares/ChromHMM/ChromHMM/ChromHMM.jar LearnModel -b 100000 -p 24 /work2/liugh/liuzunpeng/05_Results/59_Lamin_KI_DamID_part2/03_binarizebam /work2/liugh/liuzunpeng/05_Results/59_Lamin_KI_DamID_part2/04_LearnModel 2 hg19
