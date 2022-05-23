#!/usr/bin/env python3
# -*- coding:utf-8 -*-

from optparse import OptionParser
import os
import sys
import re
from collections import Counter

'''
    Title:        WGBS-seq pipeline
    Description:  This pipeline is used for WGBS-seq analysis and dowanstream analysis.
    
    Author:       Zunpeng Liu
    Email:        zunpengliuAT163.com
    Date:         23/7/2018
    Version:      1.0
    Python:       based on python 3.6
    Citation:     please
    '''

MY_USAGE = '''
    python3  WGBS.py  [options]
    
    Example: python3 WGBS.py -i /pfs1/liuguanghui/liuzunpeng/02_Data/03_WGBS-seq/04_Lamin-KI/raw_data -p M91_HHVK7CCXY_L1 \
    -o /pfs1/liuguanghui/liuzunpeng/05_Results/20_Lamin_KI
    -p /pfs1/liuguanghui/liuzunpeng/06_Qsub/20_Lamin_KI
    -c /pfs1/liuguanghui/liuzunpeng/01_Script/03_WGBS-seq/07_Lamin_KI/config.txt
    -r /pfs1/liuguanghui/liuzunpeng/03_Database/00_genome/01_hg19/03_Ensembl/01_fa_gtf_bed_chr/hg19.fa
    -g /pfs1/liuguanghui/liuzunpeng/03_Database/00_genome/01_hg19/03_Ensembl/01_fa_gtf_bed_chr/Homo_sapiens.GRCh37.87_chr.chr.gtf
    -x /pfs1/liuguanghui/liuzunpeng/03_Database/00_genome/01_hg19/03_Ensembl/01_fa_gtf_bed_chr/02_bowtie2_index/hg19
    -b /pfs1/liuguanghui/liuzunpeng/03_Database/00_genome/01_hg19/04_bed/hg19_RefSeq.bed
    
    Author:       Zunpeng Liu
    Email:        zunpengliu@163.com
    Date:         16/5/2018
    
    Blessings:    If you have any question or query of this pipeline, please feel free to contact with me.
    '''

parser = OptionParser(MY_USAGE, version='Version 1.1')
parser.add_option('-i', '--in_dir', dest='in_dir', type='string',
                  help='You should sepecify the input file path, wich stored all raw data.')
parser.add_option('-p', '--prefix', dest='prefix', type='string',
                  help='You should specifty the prefix of samples, which will be used to label all processed results. For example, if the prefix was defined as "M91_HHVK7CCXY_L1", the prefix of file names will be precessed as "M91".')
parser.add_option('-o', '--out_dir', dest='out_dir', type='string',
                  help='You should specifty the output directory, which will store all processed results.')
parser.add_option('-q', '--qsub_dir', dest='qsub_dir', type='string',
                  help='You should specifty the qub directory, which will store all shell script derived from this pipeline.')
parser.add_option('-c', '--configure_file', dest='configure_file', type='string',
                  default='/BIGDATA/grid029/01_Script/01_WGBS_LZP/config.txt',
                  help='The configure lib file should include the absolute path of softwares used in this pipeline')

(options, args) = parser.parse_args()

''' Predefine variables '''
in_dir = str(options.in_dir)
prefix = str(options.prefix)
sample = prefix
# sample  = str(re.findall(r'(\w+?)_\w+$',prefix)[0])
out_dir = str(options.out_dir)
qsu_dir = str(options.qsub_dir)

''' ensure all directory has been created, if not exist,then creat the new ones. '''
def ensure_dir(directory):
    # directory = os.path.dirname(f)
    if not os.path.exists(directory):
        os.makedirs(directory)


''' make new directions '''
ensure_dir(out_dir)
ensure_dir(qsu_dir)

''' Configure all softwares and library '''
# configure all softwares and sources from configure.txt
sft = {}
config = options.configure_file
with open(config, 'r') as config:
    for line in config:
        #        if len(line)!=0 and not line.startswith('#'):
        if len(re.findall('=', line)) == 1:
            sft[re.findall(r'(\w+?)=.*', line)[0]] = re.findall(r'\w+=(.*)', line)[0]
globals().update(sft)

''' Raw data files '''
fq1 = in_dir + '/' + prefix + '_1.clean.fq.gz'
fq2 = in_dir + '/' + prefix + '_2.clean.fq.gz'

''' 1. FastQC rawdata '''
ensure_dir(str(out_dir + '/' + sample + '/01_fastQC_raw'))
ensure_dir(str(qsu_dir + '/01_fastQC_raw'))

fqv1 = out_dir + '/' + sample + '/01_fastQC_raw/' + sample + '.R1.basic.xls'
fqv2 = out_dir + '/' + sample + '/01_fastQC_raw/' + sample + '.R2.basic.xls'

FastQC_str = '#!/bin/sh \n\
echo "Start fastQC and fqvalue raw data of ' + sample + '" && \n' \
+ fastqc + ' --outdir ' + out_dir + '/' + sample + '/01_fastQC_raw -f fastq ' + fq1 + ' ' + fq2 + ' && \n' \
+ fqvalue + ' -q 33 ' + fq1 + ' > ' + fqv1 + ' && \n' \
+ fqvalue + ' -q 33 ' + fq2 + ' > ' + fqv2 + ' && \n' \
+ 'echo "FastQC and fqvalue raw data of ' + sample + ' has been done! " '

FastQC_shell_file = open(qsu_dir + '/01_fastQC_raw/fastQC_raw.' + sample + '.sh', 'w')
FastQC_shell_file.write(FastQC_str)
FastQC_shell_file.close()

''' 2. fastp clean '''
ensure_dir(out_dir + '/' + sample + '/02_fastp')
ensure_dir(qsu_dir + '/02_fastp')

clean_fq1 = out_dir + '/' + sample + '/02_fastp/' + prefix + '_clean_1.fq.gz'
clean_fq2 = out_dir + '/' + sample + '/02_fastp/' + prefix + '_clean_2.fq.gz'

fastp_log = out_dir + '/' + sample + '/02_fastp/' + sample + '.fastp.log'

html_out = out_dir + '/' + sample + '/02_fastp/' + prefix + '.html'
json_out = out_dir + '/' + sample + '/02_fastp/' + prefix + '.json'

fastp_str = '#!/bin/sh \n\
echo "Start trim raw data of ' + sample + '" && \n' \
+ fastp + ' -i ' + fq1 + ' -I ' + fq2 + ' --thread 8 -o ' + clean_fq1 + ' -O ' + clean_fq2 + ' --html ' + html_out + ' --json ' + json_out + ' -R ' + sample + ' 2>' + fastp_log + ' && \n' \
+ 'echo "Clean and FastQC of ' + sample + ' done! All done! \"'
# thread suggest 1-8
fastp_shell_file = open(qsu_dir + '/02_fastp/fastp.' + sample + '.sh', 'w')
fastp_shell_file.write(fastp_str)
fastp_shell_file.close()


''' 3. bsmap '''
ensure_dir(out_dir + '/' + sample + '/03_bsmap')
ensure_dir(qsu_dir + '/03_bsmap')

sam_file = out_dir + '/' + sample + '/03_bsmap/' + sample + '.sam'

bsmap_log = out_dir + '/' + sample + '/03_bsmap/' + sample + '.bsmap.log'

bsmap_str = '#!/bin/bash\n' \
+ 'echo "Start mapping raw data of ' + sample + '" && \n' \
+ bsmap + ' -a ' + clean_fq1 + ' -b ' + clean_fq2 + ' -d ' + ref + ' -v 0.1 -g 1 -p 7 -R -u -o ' + sam_file + ' 2>' + bsmap_log + ' && \n' \
+ 'echo "Alignment to reference of ' + sample + ' by bsmap has been done! " '

bsmap_shell_file = open(qsu_dir + '/03_bsmap/bsmap.' + sample + '.sh', 'w')
bsmap_shell_file.write(bsmap_str)
bsmap_shell_file.close()
os.system('chmod 755 ' + qsu_dir + '/03_bsmap/bsmap.' + sample + '.sh')


''' 5. Reads quality control, and ChrY reads,then finally get unique mapped reads '''
ensure_dir(out_dir+'/'+sample+'/05_unique_mapped_reads')
ensure_dir(qsu_dir+'/05_unique_mapped_reads')

bam=out_dir+'/'+sample+'/05_unique_mapped_reads/'+sample+'_no_chrY.bam'
unique_bam=out_dir+'/'+sample+'/05_unique_mapped_reads/'+sample+'_no_chrY_unique.bam'
srt_bam=out_dir+'/'+sample+'/05_unique_mapped_reads/'+sample+'_no_chrY.sort.bam'
srt_rmdup_bam=out_dir+'/'+sample+'/05_unique_mapped_reads/'+sample+'_no_chrY.sort.rmdup.bam'
METRICS_FILE=out_dir+'/'+sample+'/05_unique_mapped_reads/'+sample+'.picard.matrix'

depth_txt    = out_dir+'/'+sample+'/05_unique_mapped_reads/'+sample+'.sort.rmdup.depth.txt'
idxstats_txt = out_dir+'/'+sample+'/05_unique_mapped_reads/'+sample+'.sort.rmdup.idxstats.txt'
flagstat_txt = out_dir+'/'+sample+'/05_unique_mapped_reads/'+sample+'.sort.rmdup.flagstat.txt'
all_stats_txt= out_dir+'/'+sample+'/05_unique_mapped_reads/'+sample+'.sort.rmdup.all_stats.txt'

srt_bam=out_dir+'/'+sample+'/03_bsmap/'+sample+'.srt.bam'

unique_reads_str='\
#!/bin/sh\n\
echo "Start transfering the sam file of '+sample+' to bam file and remove chrY DNA!" && \n'\
+'awk \'$3!="chrY" \' '+sam_file+' | '+samtools+' view -@ 24 -S -b > '+bam+' && \n'\
+'echo "Start sorting no_chrY_unique_rmdup.bam of '+sample+' on position by samtools!" && \n'\
+samtools+' sort -@ 24 -l 9 '+bam+' -o '+srt_bam+' && \n'\
+'echo "Start removing duplicates of '+sample+' by picard MarkDuplicates tools!" && \n'\
+java+' -jar '+picard+' VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true INPUT='+srt_bam+' METRICS_FILE='+METRICS_FILE+' OUTPUT='+srt_rmdup_bam+' && \n'\
+samtools+' index -@ 24 '+srt_rmdup_bam+' && \n'\
+samtools+' flagstat -@ 24  '+srt_rmdup_bam+' > '+flagstat_txt+' && \n'\
+'echo "Samtools and Picard of '+sample+' done!" '

unique_reads_shell_file=open(qsu_dir+'/05_unique_mapped_reads/unique_mapped.'+sample+'.sh','w')
unique_reads_shell_file.write(unique_reads_str)
unique_reads_shell_file.close()
os.system('chmod 755 ' + qsu_dir+'/05_unique_mapped_reads/unique_mapped.'+sample+'.sh')


''' 6. methratio_combine '''
ensure_dir(out_dir + '/' + sample + '/06_methratio_combine')
ensure_dir(qsu_dir + '/06_methratio_combine')
srt_rmdup_bam=out_dir+'/'+sample+'/05_unique_mapped_reads/'+sample+'_no_chrY.sort.rmdup.bam'


methyratio_shell_file = open(qsu_dir + '/06_methratio_combine' + '/methratio_combine.' + sample + '.sh', 'w')

chr_list = ['chr' + str(x) for x in range(1, 23)] + ['chrX', 'chrM']

header = '#!/bin/sh\n\necho "Counting the methylated ratio of ' + sample + '" && \n'
methyratio_shell_file.write(header)

chr_dir = '/work2/liugh/liuzunpeng/03_database/00_genome/01_hg19/02_UCSC/chr/'

for chr in chr_list:
    Cout = out_dir + '/' + sample + '/06_methratio_combine/' + sample + '.' + chr + '.combine.cout'
    methyratio_log = out_dir + '/' + sample + '/06_methratio_combine/' + sample + '.' + chr + '.methratio.combine.log'
    
    chr_file = chr_dir + chr + '.fa'
    
    if chr in ['chr8', 'chr16','chrM']:
        methyratio_str = python + ' ' + meth_ratio + ' --ref ' + chr_file + ' --out ' + Cout + ' --sam-path ' + samtool_dir + ' --unique -t 2 --zero-meth --chr ' + chr + ' --combine-CpG --context CG ' + srt_rmdup_bam + ' 2>' + methyratio_log + '\nwait \n'
        methyratio_shell_file.write(methyratio_str)
    else:
        methyratio_str = python + ' ' + meth_ratio + ' --ref ' + chr_file + ' --out ' + Cout + ' --sam-path ' + samtool_dir + ' --unique -t 2 --zero-meth --chr ' + chr + ' --combine-CpG --context CG ' + srt_rmdup_bam + ' 2>' + methyratio_log + ' & \n'
        methyratio_shell_file.write(methyratio_str)

tail = 'echo "Count the methylated ratio of ' + sample + ' done! "'
methyratio_shell_file.write(tail)

methyratio_shell_file.close()
os.system('chmod 755 ' + qsu_dir + '/06_methratio_combine' + '/methratio_combine.' + sample + '.sh')


''' 7. methratio_merge_chr '''
ensure_dir(out_dir + '/' + sample + '/07_methratio_merge_chr')
ensure_dir(qsu_dir + '/07_methratio_merge_chr')

all_bed = out_dir + '/' + sample + '/07_methratio_merge_chr/' + sample + '.all_chr.combine_all.bed'
all_bw = out_dir + '/' + sample + '/07_methratio_merge_chr/' + sample + '.all_chr.combine_all.bw'

merge_chr_str = '#!/bin/sh\n' \
+ 'echo "Start merging cout of ' + sample + '" && \n' \
+ python3 + ' ' + merge_chr_methratio + ' -i ' + out_dir + '/' + sample + '/06_methratio_combine -o ' + out_dir + '/' + sample + '/07_methratio_merge_chr -p ' + sample + ' && \n' \
+ bedGraphToBigWig + ' ' + all_bed + ' ' + chromosomesize + ' ' + all_bw + ' && \n' \
+ 'echo "Merge cout of ' + sample + ' done! " '

merge_chr_shell_file = open(qsu_dir + '/07_methratio_merge_chr/merge_chr_' + sample + '.sh', 'w')
merge_chr_shell_file.write(merge_chr_str)
merge_chr_shell_file.close()
os.system('chmod 755 ' + qsu_dir + '/07_methratio_merge_chr/merge_chr_' + sample + '.sh')


''' 08_CHG_meth '''
ensure_dir(out_dir + '/' + sample + '/08_CHG_meth')
ensure_dir(qsu_dir + '/08_CHG_meth')

srt_rmdup_bam=out_dir+'/'+sample+'/05_unique_mapped_reads/'+sample+'_no_chrY.sort.rmdup.bam'


methyratio_shell_file = open(qsu_dir + '/08_CHG_meth' + '/CHG_meth.' + sample + '.sh', 'w')

chr_list = ['chr' + str(x) for x in range(1, 23)] + ['chrX', 'chrM']

header = '#!/bin/sh\n\necho "Counting the methylated ratio of ' + sample + '" && \n'
methyratio_shell_file.write(header)

chr_dir = '/work2/liugh/liuzunpeng/03_database/00_genome/01_hg19/02_UCSC/chr/'

for chr in chr_list:
    Cout = out_dir + '/' + sample + '/08_CHG_meth/' + sample + '.' + chr + '.CHG.cout'
    methyratio_log = out_dir + '/' + sample + '/08_CHG_meth/' + sample + '.' + chr + '.methratio.combine.log'
    
    chr_file = chr_dir + chr + '.fa'
    
    if chr in ['chr8', 'chr16','chrM']:
        methyratio_str = python + ' ' + meth_ratio + ' --ref ' + chr_file + ' --out ' + Cout + ' --sam-path ' + samtool_dir + ' --unique -t 2 --zero-meth --chr ' + chr + ' --context CHG ' + srt_rmdup_bam + ' 2>' + methyratio_log + '\nwait \n'
        methyratio_shell_file.write(methyratio_str)
    else:
        methyratio_str = python + ' ' + meth_ratio + ' --ref ' + chr_file + ' --out ' + Cout + ' --sam-path ' + samtool_dir + ' --unique -t 2 --zero-meth --chr ' + chr + ' --context CHG ' + srt_rmdup_bam + ' 2>' + methyratio_log + ' & \n'
        methyratio_shell_file.write(methyratio_str)

tail = 'echo "Count the methylated ratio of ' + sample + ' done! "'
methyratio_shell_file.write(tail)

methyratio_shell_file.close()
os.system('chmod 755 ' + qsu_dir + '/08_CHG_meth' + '/CHG_meth.' + sample + '.sh')


''' 09_CHG_merge_chr '''
ensure_dir(out_dir + '/' + sample + '/09_CHG_merge_chr')
ensure_dir(qsu_dir + '/09_CHG_merge_chr')

all_bed= out_dir + '/' + sample + '/09_CHG_merge_chr/' + sample + '.CHG.all_chr.combine_all.bed'
fw_bed = out_dir + '/' + sample + '/09_CHG_merge_chr/' + sample + '.CHG.all_chr.combine_fw.bed'
rev_bed = out_dir + '/' + sample + '/09_CHG_merge_chr/' + sample + '.CHG.all_chr.combine_rev.bed'
all_bw = out_dir + '/' + sample + '/09_CHG_merge_chr/' + sample + '.CHG.all_chr.combine_all.bw'
fw_bw = out_dir + '/' + sample + '/09_CHG_merge_chr/' + sample + '.CHG.all_chr.combine_fw.bw'
rev_bw = out_dir + '/' + sample + '/09_CHG_merge_chr/' + sample + '.CHG.all_chr.combine_rev.bw'


merge_chr_str = '#!/bin/sh\n' \
+ 'echo "Start merging cout of ' + sample + '" && \n' \
+ python3 + ' ' + merge_chr_methratio_CHH_CHG + ' -i ' + out_dir + '/' + sample + '/08_CHG_meth -o ' + out_dir + '/' + sample + '/09_CHG_merge_chr -p ' + sample + '.CHG && \n' \
+ bedGraphToBigWig + ' ' + all_bed + ' ' + chromosomesize + ' ' + all_bw + ' && \n' \
+ bedGraphToBigWig + ' ' + fw_bed + ' ' + chromosomesize + ' ' + fw_bw + ' && \n' \
+ bedGraphToBigWig + ' ' + rev_bed + ' ' + chromosomesize + ' ' + rev_bw + ' && \n' \
+ 'echo "Merge cout of ' + sample + ' done! " '

merge_chr_shell_file = open(qsu_dir + '/09_CHG_merge_chr/merge_chr_' + sample + '.sh', 'w')
merge_chr_shell_file.write(merge_chr_str)
merge_chr_shell_file.close()
os.system('chmod 755 ' + qsu_dir + '/09_CHG_merge_chr/merge_chr_' + sample + '.sh')


''' 10_CHH_meth '''
ensure_dir(out_dir + '/' + sample + '/10_CHH_meth')
ensure_dir(qsu_dir + '/10_CHH_meth')

methyratio_shell_file = open(qsu_dir + '/10_CHH_meth' + '/CHH_meth.' + sample + '.sh', 'w')

chr_list = ['chr' + str(x) for x in range(1, 23)] + ['chrX', 'chrM']

header = '#!/bin/sh\n\necho "Counting the methylated ratio of ' + sample + '" && \n'
methyratio_shell_file.write(header)

chr_dir = '/work2/liugh/liuzunpeng/03_database/00_genome/01_hg19/02_UCSC/chr/'

for chr in chr_list:
    Cout = out_dir + '/' + sample + '/10_CHH_meth/' + sample + '.' + chr + '.CHH.cout'
    methyratio_log = out_dir + '/' + sample + '/10_CHH_meth/' + sample + '.' + chr + '.methratio.combine.log'
    
    chr_file = chr_dir + chr + '.fa'
    
    if chr in ['chr8', 'chr16','chrM']:
        methyratio_str = python + ' ' + meth_ratio + ' --ref ' + chr_file + ' --out ' + Cout + ' --sam-path ' + samtool_dir + ' --unique -t 2 --zero-meth --chr ' + chr + ' --context CHH ' + srt_rmdup_bam + ' 2>' + methyratio_log + '\nwait \n'
        methyratio_shell_file.write(methyratio_str)
    else:
        methyratio_str = python + ' ' + meth_ratio + ' --ref ' + chr_file + ' --out ' + Cout + ' --sam-path ' + samtool_dir + ' --unique -t 2 --zero-meth --chr ' + chr + ' --context CHH ' + srt_rmdup_bam + ' 2>' + methyratio_log + ' & \n'
        methyratio_shell_file.write(methyratio_str)

tail = 'echo "Count the methylated ratio of ' + sample + ' done! "'
methyratio_shell_file.write(tail)

methyratio_shell_file.close()
os.system('chmod 755 ' + qsu_dir + '/10_CHH_meth' + '/CHH_meth.' + sample + '.sh')


''' 11_CHH_merge_chr '''
ensure_dir(out_dir + '/' + sample + '/11_CHH_merge_chr')
ensure_dir(qsu_dir + '/11_CHH_merge_chr')

all_bed= out_dir + '/' + sample + '/11_CHH_merge_chr/' + sample + '.CHH.all_chr.combine_all.bed'
fw_bed = out_dir + '/' + sample + '/11_CHH_merge_chr/' + sample + '.CHH.all_chr.combine_fw.bed'
rev_bed = out_dir + '/' + sample + '/11_CHH_merge_chr/' + sample + '.CHH.all_chr.combine_rev.bed'
all_bw = out_dir + '/' + sample + '/11_CHH_merge_chr/' + sample + '.CHH.all_chr.combine_all.bw'
fw_bw = out_dir + '/' + sample + '/11_CHH_merge_chr/' + sample + '.CHH.all_chr.combine_fw.bw'
rev_bw = out_dir + '/' + sample + '/11_CHH_merge_chr/' + sample + '.CHH.all_chr.combine_rev.bw'


merge_chr_str = '#!/bin/sh\n' \
+ 'echo "Start merging cout of ' + sample + '" && \n' \
+ python3 + ' ' + merge_chr_methratio_CHH_CHG + ' -i ' + out_dir + '/' + sample + '/10_CHH_meth -o ' + out_dir + '/' + sample + '/11_CHH_merge_chr -p ' + sample + '.CHH && \n' \
+ bedGraphToBigWig + ' ' + all_bed + ' ' + chromosomesize + ' ' + all_bw + ' && \n' \
+ bedGraphToBigWig + ' ' + fw_bed + ' ' + chromosomesize + ' ' + fw_bw + ' && \n' \
+ bedGraphToBigWig + ' ' + rev_bed + ' ' + chromosomesize + ' ' + rev_bw + ' && \n' \
+ 'echo "Merge cout of ' + sample + ' done! " '

merge_chr_shell_file = open(qsu_dir + '/11_CHH_merge_chr/merge_chr_' + sample + '.sh', 'w')
merge_chr_shell_file.write(merge_chr_str)
merge_chr_shell_file.close()
os.system('chmod 755 ' + qsu_dir + '/11_CHH_merge_chr/merge_chr_' + sample + '.sh')






''' 12_PMD_HMD '''

ensure_dir(qsu_dir + '/12_PMD_HMD')

ensure_dir(out_dir + '/' + sample + '/12_PMD_HMD')
    
k10="/work2/liugh/liuzunpeng/04_Softwares/WGBS_src/hg19.10k.bed"
d5_bed = out_dir + '/' + sample + '/07_methratio_merge_chr/' + sample + '.all_chr.combine_depth_5.bed'

k10_methylation_intersect=out_dir + '/' + sample + '/12_PMD_HMD/' + sample + '.10kb.CpG_meth.intersect.bed'
k10_methylation_groupby=out_dir + '/' + sample + '/12_PMD_HMD/' + sample + '.10kb.CpG_meth.groupby.bed'
k10_methylation_groupby_subset=out_dir + '/' + sample + '/12_PMD_HMD/' + sample + '.10kb.CpG_meth.groupby_sub.bed'
    
PMD= out_dir + '/' + sample + '/12_PMD_HMD/' + sample + '.PMD.bed'
PMD_state=out_dir + '/' + sample + '/12_PMD_HMD/' + sample + '.PMD.state.txt'
HMD = out_dir + '/' + sample + '/12_PMD_HMD/' + sample + '.HMD.bed'
PMD_anno_enrich=out_dir + '/' + sample + '/12_PMD_HMD/' + sample + '.PMD.anno.txt'
PMD_anno_txt=out_dir + '/' + sample + '/12_PMD_HMD/' + sample + '.PMD.enrichment.txt'
PMD_anno_enrich_p1=out_dir + '/' + sample + '/12_PMD_HMD/' + sample + '.PMD.anno.enrich_p1.txt'
PMD_anno_enrich_p2=out_dir + '/' + sample + '/12_PMD_HMD/' + sample + '.PMD.anno.enrich_p2.txt'
    
genome="/work2/liugh/liuzunpeng/04_Softwares/bedtools2/genomes/human.hg19.genome.txt"
PMD_HMD_str = '#!/bin/sh\n' \
+ 'echo "Start merging cout of ' + sample + '" && \n' \
+ bedtools + '/bedtools intersect -a  ' + k10 + ' -b  ' + d5_bed + ' -wa  -wb > ' + k10_methylation_intersect + ' && \n' \
+ bedtools + '/bedtools groupby -i  ' + k10_methylation_intersect + ' -g 1,2,3 -opCols 7,7,7,7 -ops min,max,mean,count > ' + k10_methylation_groupby + ' && \n' \
+ 'awk \'$7>=10 && $6<0.7\' ' + k10_methylation_groupby + ' > ' + k10_methylation_groupby_subset + ' && \n' \
+ bedtools + '/bedtools merge -i ' + k10_methylation_groupby_subset +' -c 1 -o count -d 100000 | awk \'$4>2\' > ' +PMD + ' && \n' \
+ bedtools + '/bedtools complement -i ' + PMD +' -g ' + genome + ' > ' + HMD + ' && \n' \
+ perl + ' ' + macs_stat + ' ' +PMD + ' > ' + PMD_state + ' && \n' \
+ annotatePeaks + ' ' + PMD + ' hg19 -annStats ' + PMD_anno_enrich  + ' > ' + PMD_anno_txt + ' && \n' \
+ 'head -n 12 ' + PMD_anno_enrich + ' |grep -E "3UTR|TTS|Exon|Intron|Intergenic|Promoter|5UTR|Ann" |sed \'1c Annotation\\tPeak_num\\tTotal_size\\tLog2_Enrichment\' >' +PMD_anno_enrich_p1 + ' && \n' \
+ 'tail -n 34 ' + PMD_anno_enrich + ' |grep -E "SINE|LINE|LTR|Satellite|rRNA|RC|Low_complexity|Simple_repeat|Other|DNA|RNA|snoRNA|AnnncRNA|srpRNA|tRNA|snRNA|scRNA|CpG-Island|Unknown" |grep -v ? |sed \'1c Annotation\\tPeak_num\\tTotal_size\\tLog2_Enrichment\' >' +PMD_anno_enrich_p2 + ' && \n' \
+ 'echo "Merge cout of ' + sample + ' done! " '
    
PMD_HMD_shell_file = open(qsu_dir + '/12_PMD_HMD/PMD_HMD_' + sample + '.sh', 'w')
PMD_HMD_shell_file.write(PMD_HMD_str)
PMD_HMD_shell_file.close()
os.system('chmod 755 ' + qsu_dir + '/12_PMD_HMD/PMD_HMD_' + sample + '.sh')



######################################
''' Src '''
def getbed(bed):
    global beds
    with open(bed,'r') as b:
        for lines in b:
            lines=lines.strip()
            if len(lines) != 0 and not lines.startswith('#'):
                if len(re.findall('=', lines)) == 1:
                    beds[re.findall(r'(\w+?)=.*', lines)[0]] = re.findall(r'\w+=(.*)', lines)[0]


beds = {}
getbed('/work2/liugh/liuzunpeng/01_Script/01_pipeline/02_ChIP-seq/00_ChIP_seq_pipeline/bed.lamin.txt')

''' 08_PMD_HMD '''

ensure_dir(qsu_dir + '/13_CpG_meth_of_bed')

d5_bed = out_dir + '/' + sample + '/07_methratio_merge_chr/' + sample + '.all_chr.combine_depth_5.bed'
ensure_dir(out_dir + '/' + sample +'/13_CpG_meth_of_bed')
ensure_dir(qsu_dir+'/13_CpG_meth_of_bed')
    
for bed_name, bed in beds.items():
    k2_methylation_intersect=out_dir + '/' + sample +'/13_CpG_meth_of_bed/'+sample+'.'+bed_name+'.CpG.intersect.bed'
    k2_methylation_groupby=out_dir + '/' + sample +'/13_CpG_meth_of_bed/'+sample+'.'+bed_name+'.CpG.bed'
        
    CpG_meth_of_bed_str = '#!/bin/sh\n' \
    + 'echo "Get CpG mthylation of bed ' + sample + '" && \n' \
    + bedtools + '/bedtools  intersect -a  ' + bed + ' -b  ' + d5_bed + ' -wa  -wb > ' + k2_methylation_intersect + ' && \n' \
    + bedtools + '/bedtools  groupby -i  ' + k2_methylation_intersect + ' -g 1,2,3 -opCols 10 -ops mean > ' + k2_methylation_groupby + ' && \n' \
    + 'rm '+k2_methylation_intersect+ ' && \n' \
    + 'echo "' + sample + ' done! " '
        
    CpG_meth_of_bed_shell_file = open(qsu_dir+'/13_CpG_meth_of_bed/' +bed_name+'.'+ sample  + '.sh','w')
    CpG_meth_of_bed_shell_file.write(CpG_meth_of_bed_str)
    CpG_meth_of_bed_shell_file.close()
    os.system('chmod 755 ' + qsu_dir+'/13_CpG_meth_of_bed/' +bed_name+'.'+ sample  + '.sh')


