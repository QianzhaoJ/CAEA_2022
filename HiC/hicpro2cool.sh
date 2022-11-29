#!/bin/bash

raw=XXXXXXXXXXXXXX
tar=XXXXXXXXXXXXXX

mkdir -p $tar/scripts

for sample in `ls $raw`
do

for res in {5000,10000,25000,100000}
do

echo -e '#!/bin/bash' > $tar/scripts/${sample}_${res}_cool.sh
echo '#SBATCH -N 1' >> $tar/scripts/${sample}_${res}_cool.sh
echo "#SBATCH -n 1" >> $tar/scripts/${sample}_${res}_cool.sh
echo "#SBATCH -c 6" >> $tar/scripts/${sample}_${res}_cool.sh
echo "#SBATCH -p compute" >> $tar/scripts/${sample}_${res}_cool.sh
echo "#SBATCH --job-name=${sample}_${res}" >> $tar/scripts/${sample}_${res}_cool.sh
echo '#SBATCH --export=ALL' >> $tar/scripts/${sample}_${res}_cool.sh
echo -e sample=$sample\\nraw=$raw\\ntar=$tar\\nres=$res\\n'
echo Convert: $sample is Starting

out=$tar/$sample
mkdir -p $out

input=$raw/$sample/raw/$res/${sample}_${res}.matrix
chrom=XXXXXXXXXXXX/annotation/chrom_hg19.rm.sizes

echo #### 1.create  bin file
cooler makebins $chrom $res > $out/bin${res}.bed &&
echo "#### 2.raw data to cool files" &&
cooler load -f coo --one-based $out/bin${res}.bed $input $out/${sample}.${res}.cool &&
echo "#### 3.cooler normalization" &&
cooler balance $out/${sample}.${res}.cool &&
echo Convert: $sample is Done !

' >>  $tar/scripts/${sample}_${res}_cool.sh
sbatch $tar/scripts/${sample}_${res}_cool.sh
done
done
