#!/bin/bash
tar=XXXXXXXXXXX
raw=XXXXXXXXX
fp=$tar/01_fastp

for sample in `ls $raw`
do

echo -e '#!/bin/bash' > $fp/scripts/${sample}_fp.sh
echo '#BSUB -q XX' >> $fp/scripts/${sample}_fp.sh
echo "#BSUB -o $sample.out" >> $fp/scripts/${sample}_fp.sh
echo "#BSUB -e $sample.err" >> $fp/scripts/${sample}_fp.sh
echo "#BSUB -J fastp_$sample" >> $fp/scripts/${sample}_fp.sh
echo '#BSUB -n 12' >> $fp/scripts/${sample}_fp.sh
echo '#BSUB -R "span[ptile=12]"' >> $fp/scripts/${sample}_fp.sh

echo -e sample=${sample}\\ntar=$tar\\nraw=$raw'

echo Fastp: $sample is Starting

command=/work2/liugh/wuqiong/software/fastp/fastp

fp=$tar/01_fastp

fq1=$raw/$sample/${sample}_R1.fq.gz
fq2=$raw/$sample/${sample}_R2.fq.gz

result=$fp/$sample
mkdir $result

out1=$fp/$sample/${sample}_R1.clean.fq.gz
out2=$fp/$sample/${sample}_R2.clean.fq.gz
md5=$fp/$sample/${sample}.md5

html=$fp/$sample/${sample}.html
json=$fp/$sample/${sample}.json
log=$fp/logs/${sample}.log

md5sum $fq1 > $md5
md5sum $fq2 >> $md5

$command -i $fq1 -I $fq2 -o $out1 -O $out2 --html $html --json $json --thread 12 -R $sample 2>$log  &&

echo Fastp: $sample is Done' >> $fp/scripts/${sample}_fp.sh
done

echo '#!/bin/bash
scripts=./scripts
cd $scripts
for i in `ls *.sh`

do 
   bsub < $i &

done' >$fp/run_fp.sh

