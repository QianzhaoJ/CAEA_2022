#!/bin/bash
tar=xxxxxxxxxxx
raw=xxxxxxxxxxxx

mkdir -p $tar/scripts
mkdir -p $tar/logs

for sample in `ls $raw`
do

echo -e '#!/bin/bash' > $tar/scripts/${sample}_j.sh
echo "#SBATCH -N 1" >> $tar/scripts/${sample}_j.sh
echo "#SBATCH -n 1" >> $tar/scripts/${sample}_j.sh
echo "#SBATCH -c 48" >> $tar/scripts/${sample}_j.sh
echo "#SBATCH -p compute" >> $tar/scripts/${sample}_j.sh
echo "#SBATCH --job-name Juicer_$sample" >> $tar/scripts/${sample}_j.sh
echo "#SBATCH --export=ALL" >> $tar/scripts/${sample}_j.sh
echo -e sample=${sample}\\ntar=$tar\\nraw=$raw'

echo Convert: $sample is Starting
### options
Vp=$raw/$sample/${sample}.allValidPairs
gsize=xxxxxxxxxxxxx/annotation/chrom_hg19.rm.sizes
juicer=xxxxxxxxxxxxx/juicer_tools_1.22.01.jar
resfrag=xxxxxxxxxxxxx/MboI_resfrag_hg19.bed
temp=$tar/$sample/tmp
out=$tar/$sample
mkdir -p $out/$sample
###
sh=xxxxxxx/bin/utils/hicpro2juicebox.sh
###
$sh -i $Vp -g $gsize -j $juicer -r $resfrag -t $temp -o $out 2>$tar/logs/${sample}.log

java -Xmx160g -jar $juicer pre --threads 40 -f $temp/*_resfrag.juicebox $temp/*_allValidPairs.pre_juicebox_sorted $out/${sample}.hic $gsize 2>$tar/logs/${sample}.log

echo Convert: $sample is Done' >> $tar/scripts/${sample}_j.sh
done

echo '#!/bin/bash
scripts=./scripts
cd $scripts
for i in `ls *.sh`

do 
   sbatch $i &

done' >$tar/run_j.sh

