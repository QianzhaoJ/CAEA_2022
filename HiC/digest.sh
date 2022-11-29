#!/bin/bash
#BSUB -q c_liugh2
#BSUB -o digest.out
#BSUB -o digest.err
#BSUB -J digest
#BSUB -n 1
#BSUB -R "span[ptile=1]"

command=XXXXXX/HiC-Pro_2.11.1/bin/utils/digest_genome.py
out=XXXXXXXXX/MboI_resfrag_hg19.bed
fa=XXXXXX/hg19.fa

$command -r mboi -o $out $fa
