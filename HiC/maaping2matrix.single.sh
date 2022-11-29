#!/bin/bash
#BSUB -J HiCpro_s
#BSUB -e HiCpro_s2_.e
#BSUB -o HiCpro_s2_.o
#BSUB -q c_liugh2
#BSUB -n 24

command=XXXXX/HiCpro/HiC-Pro_2.11.1/bin/HiC-Pro
input=XXXXXXXXX
output=XXXXXX
config=./config-hicpro.txt

$command -i $input -o $output -c $config -s merge_persample -s build_contact_maps -s ice_norm

