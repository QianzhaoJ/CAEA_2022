#!/bin/bash
command=XXXXX/HiCpro/HiC-Pro_2.11.1/bin/HiC-Pro
input=XXXXXXXXX
output=XXXXXX
config=./config-hicpro.txt

$command -i $input -o $output -c $config -p
