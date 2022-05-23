#!/usr/bin/env python3
# coding: utf-8

from optparse import OptionParser
import os
import sys
import re

'''
    Title:        ChIP-seq pipeline
    Description:  This pipeline is used for WGBS analysis and downstream analysis.
    
    Author:       Zunpeng Liu
    Email:        zunpengliuAT163.com
    Date:         7/5/2019
    Version:      3.0
    Python:       based on python 3.6
    Citation:     please
    '''

MY_USAGE = '''
    python3  WGBS.part_2.py [options]
    
    Example:
    
    python3 WGBS.part_2.py \
    -i /work1/liugh/liuzunpeng/05_Results/33_WGBS_progerin\
    -o /work1/liugh/liuzunpeng/05_Results/33_WGBS_progerin_part2 \
    -q /work1/liugh/liuzunpeng/06_Qsub/33_WGBS_progerin_part2 \
    -c /work2/liugh/liuzunpeng/01_Script/01_pipeline/09_WGBS/02_progerin/config.txt \
    -g H9-MSC-P3,H9-MSC-P8,R644-MSC-P3,R-MSC-P8,homo-MSC-P3,Homo-P8,WRN-MSC-P3,WS-MSC-P8 \
    -p 05_Lamin_KI_H3K9me3 \
    -b /work2/liugh/liuzunpeng/01_Script/01_pipeline/02_ChIP-seq/00_ChIP_seq_pipeline/bed.txt
    

    Author:       Zunpeng Liu
    Email:        zunpengliu@163.com
    Date:         7/5/2019
    
    Blessings:    If you have any question or query of this pipeline, please feel free to contact with me.
    '''

parser = OptionParser(MY_USAGE, version='Version 1.1')
parser.add_option('-i', '--in_dir', dest='in_dir', type='string',
                  help='You should sepecify the input file path, wich stored all raw data.')
parser.add_option('-o', '--out_dir', dest='out_dir', type='string',
                  help='You should sepecify the output file path, wich stored all raw data.')
parser.add_option('-g', '--group', dest='group', type='string',
                  help='You should specifty the prefix of samples, which will be used to label all processed results. For example, if the prefix was defined as "M91_HHVK7CCXY_L1", the prefix of file names will be precessed as "M91".')
parser.add_option('-q', '--qsub_dir', dest='qsub_dir', type='string',
                  help='You should specifty the qub directory, which will store all shell script derived from this pipeline.')
parser.add_option('-c', '--configure_file', dest='configure_file', type='string',
                  default='/pfs1/liuguanghui/liuzunpeng/01_Script/03_ATAC-seq/07_ATAC_seq_pipeline/config.txt',
                  help='The configure lib file should include the absolute path of softwares used in this pipeline')
parser.add_option('-b', '--bed_file', dest='bed_file', type='string',
                  default='/pfs1/liuguanghui/liuzunpeng/01_Script/03_ATAC-seq/07_ATAC_seq_pipeline/config.txt',
                  help='The configure lib file should include the absolute path of softwares used in this pipeline')
parser.add_option('-C', '--control', dest='control', type='string',
                  help='You should sepecify the name of control samples.')
parser.add_option('-T', '--treat', dest='treat', type='string',
                  help='You should sepecify the name of treat samples.')
parser.add_option('-p', '--prefix', dest='prefix', type='string',
                  help='You should sepecify the name of treat samples.')

(options, args) = parser.parse_args()

''' Predefine variables '''
in_dir  = options.in_dir
group   = options.group
out_dir = options.out_dir
qsu_dir = options.qsub_dir
#control = options.control
#treat   = options.treat
prefix  = options.prefix
bed_file= options.bed_file

''' ensure all directory has been created, if not exist,then creat the new ones. '''

def ensure_dir(directory):
    # directory = os.path.dirname(f)
    if not os.path.exists(directory):
        os.makedirs(directory)

lsd = os.listdir(in_dir)
lsd.sort()

groups=group.split(',')
samples={}
#{g1:[g1_a,g1_b,g1_c],g2:{g2_e,g2_f,g2_g}}

for g in groups :
    for sam in lsd :
        print(sam)
        print(g)
        if sam.startswith(g):
            if g in samples.keys():
                samples[g].append(sam)
            else:
                samples[g]=[sam]
for k,v in samples.items():
    print(k,v)

''' make new directions '''
ensure_dir(out_dir)
ensure_dir(qsu_dir)

''' Src '''

''' Configure all softwares and library '''

sft = {}
config = options.configure_file
with open(config, 'r') as config:
    for line in config:
        #        if len(line)!=0 and not line.startswith('#'):
        if len(re.findall('=', line)) == 1:
            sft[re.findall(r'(\w+?)=.*', line)[0]] = re.findall(r'\w+=(.*)', line)[0]
globals().update(sft)


''' 1. Merge bams from for each group '''
ensure_dir(qsu_dir + '/01_merge_bam')

for sample in samples.keys():
    ensure_dir(out_dir + '/' + sample + '/01_merge_bam')
    merge_bam_log = out_dir + '/' + sample + '/01_merge_bam/' + sample + '.bmerge.bam.log.txt'
    merge_bam = out_dir + '/' + sample + '/01_merge_bam/' + sample + '.merge.bam'
    srt_bam   = out_dir + '/' + sample + '/01_merge_bam/' + sample + '.merge.srt.bam'
    flagstat_txt = out_dir+'/'+sample+'/01_merge_bam/'+sample+'.merge.srt.flagstat.txt'

    bams=[in_dir + '/' + i + '/05_unique_mapped_reads/' + i + '_no_chrY.sort.rmdup.bam' for i in samples[sample]]
    flagstat_txt = out_dir+'/'+sample+'/05_unique_mapped_reads/'+sample+'.sort.rmdup.flagstat.txt'

    merge_bam_str = '#/bin/sh\necho "Start merging bam of ' + sample + '!" && \n' \
    + samtools + ' merge -@ 20 ' + merge_bam + ' ' + ' '.join(bams) + ' 2>' + merge_bam_log + '  && \n' \
    +samtools+' sort -@ 20 -l 9 '+merge_bam+' -o '+srt_bam+' && \n'\
    +samtools+' index -@ 20 '+srt_bam+' && \n'\
    +samtools+' flagstat -@ 20  '+srt_bam+' > '+flagstat_txt+' && \n'\
    + 'echo "Merge bam ' + sample + ' done! \"'
    
    merge_bam_shell_file = open(qsu_dir + '/01_merge_bam/merge_bam.' + sample + '.sh', 'w')
    merge_bam_shell_file.write(merge_bam_str)
    merge_bam_shell_file.close()
    os.system('chmod 755 ' + qsu_dir + '/01_merge_bam/merge_bam.' + sample + '.sh')

''' 02_methratio_CpG_combine '''
ensure_dir(qsu_dir + '/02_methratio_CpG_combine')

chr_list = ['chr' + str(x) for x in range(1, 22)] + ['chrX', 'chrM']

for sample in samples.keys():
    srt_bam   = out_dir + '/' + sample + '/01_merge_bam/' + sample + '.merge.srt.bam'
    ensure_dir(out_dir + '/' + sample + '/02_methratio_CpG_combine')
    
    methyratio_shell_file = open(qsu_dir + '/02_methratio_CpG_combine' + '/methratio_combine.' + sample + '.sh', 'w')
    header = '#!/bin/sh\n\necho "Counting the methylated ratio of ' + sample + '" && \n'
    methyratio_shell_file.write(header)
    chr_dir = '/work2/liugh/liuzunpeng/03_database/00_genome/01_hg19/02_UCSC/chr/'

    for chr in chr_list:
        Cout = out_dir + '/' + sample + '/02_methratio_CpG_combine/' + sample + '.' + chr + '.combine.cout'
        methyratio_log = out_dir + '/' + sample + '/02_methratio_CpG_combine/' + sample + '.' + chr + '.methratio.combine.log'
    
        chr_file = chr_dir + chr + '.fa'
    
        if chr in ['chr8', 'chr16','chrM']:
            methyratio_str = python + ' ' + meth_ratio + ' --ref ' + chr_file + ' --out ' + Cout + ' --sam-path ' + samtool_dir + ' --unique -t 2 --zero-meth --chr ' + chr + ' --combine-CpG --context CG ' + srt_bam + ' 2>' + methyratio_log + '\nwait \n'
            methyratio_shell_file.write(methyratio_str)
        else:
            methyratio_str = python + ' ' + meth_ratio + ' --ref ' + chr_file + ' --out ' + Cout + ' --sam-path ' + samtool_dir + ' --unique -t 2 --zero-meth --chr ' + chr + ' --combine-CpG --context CG ' + srt_bam + ' 2>' + methyratio_log + ' & \n'
            methyratio_shell_file.write(methyratio_str)

    tail = 'echo "Count the methylated ratio of ' + sample + ' done! "'
    methyratio_shell_file.write(tail)

    methyratio_shell_file.close()
    os.system('chmod 755 ' + qsu_dir + '/02_methratio_CpG_combine' + '/methratio_combine.' + sample + '.sh')


''' 03_CpG_methratio_merge_chr '''
ensure_dir(qsu_dir + '/03_CpG_methratio_merge_chr')

for sample in samples.keys():
    ensure_dir(out_dir + '/' + sample + '/03_CpG_methratio_merge_chr')

    all_bed = out_dir + '/' + sample + '/03_CpG_methratio_merge_chr/' + sample + '.all_chr.combine_all.bed'
    all_bw = out_dir + '/' + sample + '/03_CpG_methratio_merge_chr/' + sample + '.all_chr.combine_all.bw'

    merge_chr_str = '#!/bin/sh\n' \
    + 'echo "Start merging cout of ' + sample + '" && \n' \
    + python3 + ' ' + merge_chr_methratio + ' -i ' + out_dir + '/' + sample + '/02_methratio_CpG_combine -o ' + out_dir + '/' + sample + '/03_CpG_methratio_merge_chr -p ' + sample + ' && \n' \
    + bedGraphToBigWig + ' ' + all_bed + ' ' + chromosomesize + ' ' + all_bw + ' && \n' \
    + 'echo "Merge cout of ' + sample + ' done! " '

    merge_chr_shell_file = open(qsu_dir + '/03_CpG_methratio_merge_chr/merge_chr_' + sample + '.sh', 'w')
    merge_chr_shell_file.write(merge_chr_str)
    merge_chr_shell_file.close()
    os.system('chmod 755 ' + qsu_dir + '/03_CpG_methratio_merge_chr/merge_chr_' + sample + '.sh')


''' 04_methratio_CHG '''
ensure_dir(qsu_dir + '/04_methratio_CHG')

chr_list = ['chr' + str(x) for x in range(1, 23)] + ['chrX', 'chrM']

for sample in samples.keys():
    ensure_dir(out_dir + '/' + sample + '/04_methratio_CHG')
    
    srt_bam   = out_dir + '/' + sample + '/01_merge_bam/' + sample + '.merge.srt.bam'
    methyratio_shell_file = open(qsu_dir + '/04_methratio_CHG' + '/04_methratio_CHG.' + sample + '.sh', 'w')
    header = '#!/bin/sh\n\necho "Counting the methylated ratio of ' + sample + '" && \n'
    methyratio_shell_file.write(header)
    chr_dir = '/work2/liugh/liuzunpeng/03_database/00_genome/01_hg19/02_UCSC/chr/'
    
    for chr in chr_list:
        Cout = out_dir + '/' + sample + '/04_methratio_CHG/' + sample + '.' + chr + '.CHG.cout'
        methyratio_log = out_dir + '/' + sample + '/04_methratio_CHG/' + sample + '.' + chr + '.methratio.CHG.log'
        
        chr_file = chr_dir + chr + '.fa'
        
        if chr in ['chr8', 'chr16','chrM']:
            methyratio_str = python + ' ' + meth_ratio + ' --ref ' + chr_file + ' --out ' + Cout + ' --sam-path ' + samtool_dir + ' --unique -t 2 --zero-meth --chr ' + chr + ' --context CHG ' + srt_bam + ' 2>' + methyratio_log + '\nwait \n'
            methyratio_shell_file.write(methyratio_str)
        else:
            methyratio_str = python + ' ' + meth_ratio + ' --ref ' + chr_file + ' --out ' + Cout + ' --sam-path ' + samtool_dir + ' --unique -t 2 --zero-meth --chr ' + chr + ' --context CHG ' + srt_bam + ' 2>' + methyratio_log + ' & \n'
            methyratio_shell_file.write(methyratio_str)

    tail = 'echo "Count the methylated ratio of ' + sample + ' done! "'
    methyratio_shell_file.write(tail)
    
    methyratio_shell_file.close()
    os.system('chmod 755 ' + qsu_dir + '/04_methratio_CHG' + '/04_methratio_CHG.' + sample + '.sh')


''' 05_CHG_merge_chr '''
ensure_dir(qsu_dir + '/05_CHG_merge_chr')

for sample in samples.keys():
    ensure_dir(out_dir + '/' + sample + '/05_CHG_merge_chr')
    
    all_bed= out_dir + '/' + sample + '/05_CHG_merge_chr/' + sample + '.CHG.all_chr.combine_all.bed'
    fw_bed = out_dir + '/' + sample + '/05_CHG_merge_chr/' + sample + '.CHG.all_chr.combine_fw.bed'
    rev_bed = out_dir + '/' + sample + '/05_CHG_merge_chr/' + sample + '.CHG.all_chr.combine_rev.bed'
    all_bw = out_dir + '/' + sample + '/05_CHG_merge_chr/' + sample + '.CHG.all_chr.combine_all.bw'
    fw_bw = out_dir + '/' + sample + '/05_CHG_merge_chr/' + sample + '.CHG.all_chr.combine_fw.bw'
    rev_bw = out_dir + '/' + sample + '/05_CHG_merge_chr/' + sample + '.CHG.all_chr.combine_rev.bw'

    merge_chr_str = '#!/bin/sh\n' \
    + 'echo "Start merging cout of ' + sample + '" && \n' \
    + python3 + ' ' + merge_chr_methratio_CHH_CHG + ' -i ' + out_dir + '/' + sample + '/04_methratio_CHG -o ' + out_dir + '/' + sample + '/05_CHG_merge_chr -p ' + sample + '.CHG && \n' \
    + bedGraphToBigWig + ' ' + all_bed + ' ' + chromosomesize + ' ' + all_bw + ' && \n' \
    + bedGraphToBigWig + ' ' + fw_bed + ' ' + chromosomesize + ' ' + fw_bw + ' && \n' \
    + bedGraphToBigWig + ' ' + rev_bed + ' ' + chromosomesize + ' ' + rev_bw + ' && \n' \
    + 'echo "Merge cout of ' + sample + ' done! " '

    merge_chr_shell_file = open(qsu_dir + '/05_CHG_merge_chr/merge_chr_' + sample + '.sh', 'w')
    merge_chr_shell_file.write(merge_chr_str)
    merge_chr_shell_file.close()
    os.system('chmod 755 ' + qsu_dir + '/05_CHG_merge_chr/merge_chr_' + sample + '.sh')



''' 06_methratio_CHH '''
ensure_dir(qsu_dir + '/06_methratio_CHH')

for sample in samples.keys():
    ensure_dir(out_dir + '/' + sample + '/06_methratio_CHH')
    srt_bam   = out_dir + '/' + sample + '/01_merge_bam/' + sample + '.merge.srt.bam'
    
    methyratio_shell_file = open(qsu_dir + '/06_methratio_CHH' + '/methratio_CHH.' + sample + '.sh', 'w')
    header = '#!/bin/sh\n\necho "Counting the methylated ratio of ' + sample + '" && \n'
    methyratio_shell_file.write(header)
    chr_dir = '/work2/liugh/liuzunpeng/03_database/00_genome/01_hg19/02_UCSC/chr/'
    
    for chr in chr_list:
        Cout = out_dir + '/' + sample + '/06_methratio_CHH/' + sample + '.' + chr + '.CHH.cout'
        methyratio_log = out_dir + '/' + sample + '/06_methratio_CHH/' + sample + '.' + chr + '.methratio.CHH.log'
        
        chr_file = chr_dir + chr + '.fa'
        
        if chr in ['chr8', 'chr16','chrM']:
            methyratio_str = python + ' ' + meth_ratio + ' --ref ' + chr_file + ' --out ' + Cout + ' --sam-path ' + samtool_dir + ' --unique -t 2 --zero-meth --chr ' + chr + ' --context CHH ' + srt_bam + ' 2>' + methyratio_log + '\nwait \n'
            methyratio_shell_file.write(methyratio_str)
        else:
            methyratio_str = python + ' ' + meth_ratio + ' --ref ' + chr_file + ' --out ' + Cout + ' --sam-path ' + samtool_dir + ' --unique -t 2 --zero-meth --chr ' + chr + ' --context CHH ' + srt_bam + ' 2>' + methyratio_log + ' & \n'
            methyratio_shell_file.write(methyratio_str)

    tail = 'echo "Count the methylated ratio of ' + sample + ' done! "'
    methyratio_shell_file.write(tail)
    
    methyratio_shell_file.close()
    os.system('chmod 755 ' + qsu_dir + '/06_methratio_CHH' + '/methratio_CHH.' + sample + '.sh')



''' 07_CHH_merge_chr '''

ensure_dir(qsu_dir + '/07_CHH_merge_chr')

for sample in samples.keys():
    ensure_dir(out_dir + '/' + sample + '/07_CHH_merge_chr')
    
    all_bed= out_dir + '/' + sample + '/07_CHH_merge_chr/' + sample + '.CHH.all_chr.combine_all.bed'
    fw_bed = out_dir + '/' + sample + '/07_CHH_merge_chr/' + sample + '.CHH.all_chr.combine_fw.bed'
    rev_bed = out_dir + '/' + sample + '/07_CHH_merge_chr/' + sample + '.CHH.all_chr.combine_rev.bed'
    all_bw = out_dir + '/' + sample + '/07_CHH_merge_chr/' + sample + '.CHH.all_chr.combine_all.bw'
    fw_bw = out_dir + '/' + sample + '/07_CHH_merge_chr/' + sample + '.CHH.all_chr.combine_fw.bw'
    rev_bw = out_dir + '/' + sample + '/07_CHH_merge_chr/' + sample + '.CHH.all_chr.combine_rev.bw'


    merge_chr_str = '#!/bin/sh\n' \
    + 'echo "Start merging cout of ' + sample + '" && \n' \
    + python3 + ' ' + merge_chr_methratio_CHH_CHG + ' -i ' + out_dir + '/' + sample + '/06_methratio_CHH -o ' + out_dir + '/' + sample + '/07_CHH_merge_chr -p ' + sample + '.CHH && \n' \
    + bedGraphToBigWig + ' ' + all_bed + ' ' + chromosomesize + ' ' + all_bw + ' && \n' \
    + bedGraphToBigWig + ' ' + fw_bed + ' ' + chromosomesize + ' ' + fw_bw + ' && \n' \
    + bedGraphToBigWig + ' ' + rev_bed + ' ' + chromosomesize + ' ' + rev_bw + ' && \n' \
    + 'echo "Merge cout of ' + sample + ' done! " '

    merge_chr_shell_file = open(qsu_dir + '/07_CHH_merge_chr/merge_chr_' + sample + '.sh', 'w')
    merge_chr_shell_file.write(merge_chr_str)
    merge_chr_shell_file.close()
    os.system('chmod 755 ' + qsu_dir + '/07_CHH_merge_chr/merge_chr_' + sample + '.sh')





''' 08_PMD_HMD '''

ensure_dir(qsu_dir + '/08_PMD_HMD')

for sample in samples.keys():
    ensure_dir(out_dir + '/' + sample + '/08_PMD_HMD')
    
    k10="/work2/liugh/liuzunpeng/04_Softwares/WGBS_src/hg19.10k.bed"
    d5_bed = out_dir + '/' + sample + '/03_CpG_methratio_merge_chr/' + sample + '.all_chr.combine_depth_5.bed'
    
    k10_methylation_intersect=out_dir + '/' + sample + '/08_PMD_HMD/' + sample + '.10kb.CpG_meth.intersect.bed'
    k10_methylation_groupby=out_dir + '/' + sample + '/08_PMD_HMD/' + sample + '.10kb.CpG_meth.groupby.bed'
    k10_methylation_groupby_subset=out_dir + '/' + sample + '/08_PMD_HMD/' + sample + '.10kb.CpG_meth.groupby_sub.bed'
    
    PMD= out_dir + '/' + sample + '/08_PMD_HMD/' + sample + '.PMD.bed'
    PMD_state=out_dir + '/' + sample + '/08_PMD_HMD/' + sample + '.PMD.state.txt'
    HMD = out_dir + '/' + sample + '/08_PMD_HMD/' + sample + '.HMD.bed'
    PMD_anno_enrich=out_dir + '/' + sample + '/08_PMD_HMD/' + sample + '.PMD.anno.txt'
    PMD_anno_txt=out_dir + '/' + sample + '/08_PMD_HMD/' + sample + '.PMD.enrichment.txt'
    PMD_anno_enrich_p1=out_dir + '/' + sample + '/08_PMD_HMD/' + sample + '.PMD.anno.enrich_p1.txt'
    PMD_anno_enrich_p2=out_dir + '/' + sample + '/08_PMD_HMD/' + sample + '.PMD.anno.enrich_p2.txt'
    
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
    
    PMD_HMD_shell_file = open(qsu_dir + '/08_PMD_HMD/PMD_HMD_' + sample + '.sh', 'w')
    PMD_HMD_shell_file.write(PMD_HMD_str)
    PMD_HMD_shell_file.close()
    os.system('chmod 755 ' + qsu_dir + '/08_PMD_HMD/PMD_HMD_' + sample + '.sh')



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
getbed('/work2/liugh/liuzunpeng/01_Script/01_pipeline/02_ChIP-seq/00_ChIP_seq_pipeline/bed.txt')

''' 08_PMD_HMD '''

ensure_dir(qsu_dir + '/09_CpG_meth_of_bed')

for sample in samples.keys():
    d5_bed = out_dir + '/' + sample + '/03_CpG_methratio_merge_chr/' + sample + '.all_chr.combine_depth_5.bed'
    ensure_dir(out_dir + '/' + sample +'/09_CpG_meth_of_bed')
    ensure_dir(qsu_dir+'/09_CpG_meth_of_bed')
    
    for bed_name, bed in beds.items():
        k2_methylation_intersect=out_dir + '/' + sample +'/09_CpG_meth_of_bed/'+sample+'.'+bed_name+'.CpG.intersect.bed'
        k2_methylation_groupby=out_dir + '/' + sample +'/09_CpG_meth_of_bed/'+sample+'.'+bed_name+'.CpG.bed'
        
        CpG_meth_of_bed_str = '#!/bin/sh\n' \
        + 'echo "Get CpG mthylation of bed ' + sample + '" && \n' \
        + bedtools + '/bedtools  intersect -a  ' + bed + ' -b  ' + d5_bed + ' -wa  -wb > ' + k2_methylation_intersect + ' && \n' \
        + bedtools + '/bedtools  groupby -i  ' + k2_methylation_intersect + ' -g 1,2,3 -opCols 10 -ops mean > ' + k2_methylation_groupby + ' && \n' \
        + 'rm '+k2_methylation_intersect+ ' && \n' \
        + 'echo "' + sample + ' done! " '
        
        CpG_meth_of_bed_shell_file = open(qsu_dir+'/09_CpG_meth_of_bed/' +bed_name+'.'+ sample  + '.sh','w')
        CpG_meth_of_bed_shell_file.write(CpG_meth_of_bed_str)
        CpG_meth_of_bed_shell_file.close()
        os.system('chmod 755 ' + qsu_dir+'/09_CpG_meth_of_bed/' +bed_name+'.'+ sample  + '.sh')



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
getbed('/work2/liugh/liuzunpeng/01_Script/01_pipeline/02_ChIP-seq/00_ChIP_seq_pipeline/bed.progerin.txt')

getbed(bed_file)

''' 8. computeMatrix plotHeatmap plotProfile '''
#out_dir + '/' + sample + '/07_Deeptools_bamCoverage/02_bin_100/' + sample + '.no_chrY_unique.sort.rmdup.bin_100.bw'

for sample in samples.keys():
    all_bw = out_dir + '/' + sample + '/03_CpG_methratio_merge_chr/' + sample + '.all_chr.combine_all.bw'
    for bed_name, bed in beds.items():
        ensure_dir(out_dir+'/001_computeMatrix_plotHeatmap/'+bed_name)
        ensure_dir(qsu_dir+'/001_computeMatrix_plotHeatmap')
        matrix=out_dir+'/001_computeMatrix_plotHeatmap/'+bed_name+'/'+sample+'.matrix.mat.gz'
        FileSortedRegions=out_dir+'/001_computeMatrix_plotHeatmap/'+bed_name+'/'+sample+'.'+bed_name+'.Sorted.Regions.bed'
        computeMatrix_log=out_dir+'/001_computeMatrix_plotHeatmap/'+bed_name+'/'+sample+'.'+bed_name+'.computeMatrix.log'
        profile_data =out_dir+'/001_computeMatrix_plotHeatmap/'+bed_name+'/'+sample+'.'+bed_name+'.profile.tab'
        heatmap_pdf =out_dir+'/001_computeMatrix_plotHeatmap/'+bed_name+'/'+sample+'.'+bed_name+'.heatmap.pdf'
        profile_pdf =out_dir+'/001_computeMatrix_plotHeatmap/'+bed_name+'/'+sample+'.'+bed_name+'.profile.pdf'
        
        computeMatrix_str='#!/bin/sh\necho "Start deeptools computeMatrix!" && \n'\
        +computeMatrix+' scale-regions --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 -bs 50 -R ' + bed + ' -S ' + all_bw  +' --skipZeros -o ' + matrix + ' --outFileSortedRegions ' + FileSortedRegions+' 2> '+computeMatrix_log+' && \n'\
        +'echo "Finish computeMatrix of '+bed_name+' at bin size 100 bp!"  && \n'\
        +plotHeatmap + ' -m ' + matrix + ' --colorMap seismic -out ' + heatmap_pdf + ' --plotTitle "' + sample + '" --startLabel "' + bed_name + ' start" --endLabel "' + bed_name +' end" && \n'\
        +plotProfile + ' -m ' + matrix + ' --outFileNameData ' + profile_data + ' --colors blue --perGroup -out ' + profile_pdf +' --plotTitle "'+sample+ '" --startLabel "' + bed_name + ' start" --endLabel "' + bed_name +' end" && \n'\
        +'echo "Finish deeptools plotHeatmap and plotProfile!" '
        
        computeMatrix_shell_file=open(qsu_dir+'/001_computeMatrix_plotHeatmap/'+bed_name+'.'+sample+'.sh','w')
        computeMatrix_shell_file.write(computeMatrix_str)
        computeMatrix_shell_file.close()
        os.system('chmod 755 ' + qsu_dir+'/001_computeMatrix_plotHeatmap/'+bed_name+'.'+sample+'.sh')


