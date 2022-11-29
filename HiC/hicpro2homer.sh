#!/bin/bash
path=XXXXXXXXXXX
out=XXXXXXXXXXXXXX  
mkdir -p $out/convert/logs
mkdir -p $out/convert/scripts
for sample in `ls $path`
do

echo -e '#!/bin/bash' > $out/convert/scripts/${sample}_c.sh
echo '#SBATCH -N 1' >> $out/convert/scripts/${sample}_c.sh
echo "#SBATCH -n 1" >> $out/convert/scripts/${sample}_c.sh
echo "#SBATCH -c 24" >> $out/convert/scripts/${sample}_c.sh
echo "#SBATCH -p compute" >> $out/convert/scripts/${sample}_c.sh
echo "#SBATCH --job-name=${sample}_convert" >> $out/convert/scripts/${sample}_c.sh
echo '#SBATCH --export=ALL' >> $out/convert/scripts/${sample}_c.sh
echo -e sample=${sample}\\nout=$out\\npath=$path'

ori=$path/$sample/${sample}.allValidPairs
data=$path/$sample/${sample}.summary

cut -f 1-7 $ori > $data

mkdir -p $out/convert/$sample

makeTagDirectory $out/convert/$sample -format HiCsummary $data && rm $data

'>> $out/convert/scripts/${sample}_c.sh 

done

echo '#!/bin/bash
scripts=./scripts
cd $scripts
for i in `ls *_c.sh`

do 
   sbatch $i &

done' >$out/convert/run_convert.sh
