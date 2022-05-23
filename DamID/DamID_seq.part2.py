#!/usr/bin/env python3
# coding: utf-8

from optparse import OptionParser
import os
import sys
import re

'''
    Title:        ChIP-seq pipeline
    Description:  This pipeline is used for ChIP-seq analysis and downstream analysis.
    
    Author:       Zunpeng Liu
    Email:        zunpengliuAT163.com
    Date:         22/3/2019
    Version:      2.0
    Python:       based on python 3.6
    Citation:     please
    '''

MY_USAGE = '''
    python3  ChIP_seq_v1.py  [options]
    
    Example:
    
    python3 DamID_seq.part2.py \
    -i /work2/liugh/liuzunpeng/05_Results/59_Lamin_KI_DamID \
    -o /work2/liugh/liuzunpeng/05_Results/59_Lamin_KI_DamID_part2 \
    -q /work2/liugh/liuzunpeng/06_Qsub/59_Lamin_KI_DamID_part2 \
    -c /work2/liugh/liuzunpeng/01_Script/01_pipeline/02_ChIP-seq/00_ChIP_seq_pipeline/config.txt \
    -g Dam-P3-H9,Dam-P3-Homo,Dam-P3-R,Dam-P3-WS,Dam-P9-H9,Dam-P9-Homo,Dam-P9-R,Dam-P9-WS,EMD-P3-H9,EMD-P3-Homo,EMD-P3-R,EMD-P3-WS,EMD-P9-H9,EMD-P9-Homo,EMD-P9-R,EMD-P9-WS \
    -p Lamin_KI_DamID \
    -b /work2/liugh/liuzunpeng/01_Script/01_pipeline/02_ChIP-seq/00_ChIP_seq_pipeline/bed.progerin.txt
    

    python3 DamID_seq.part2.py \
    -i /work2/liugh/liuzunpeng/05_Results/59_Lamin_KI_DamID \
    -o /work2/liugh/liuzunpeng/05_Results/59_Lamin_KI_DamID_part2 \
    -q /work2/liugh/liuzunpeng/06_Qsub/59_Lamin_KI_DamID_part2 \
    -c /work2/liugh/liuzunpeng/01_Script/01_pipeline/02_ChIP-seq/00_ChIP_seq_pipeline/config.txt \
    -g EMD-P3-H9,EMD-P3-Homo,EMD-P3-R,EMD-P3-WS,EMD-P9-H9,EMD-P9-Homo,EMD-P9-R,EMD-P9-WS \
    -p Lamin_KI_DamID \
    -b /work2/liugh/liuzunpeng/01_Script/01_pipeline/02_ChIP-seq/00_ChIP_seq_pipeline/bed.progerin.txt
    
    
    
    Author:       Zunpeng Liu
    Email:        zunpengliu@163.com
    Date:         22/3/2019
    
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
    bams=[in_dir + '/' + i + '/05_unique_mapped_reads/' + i + '_no_chrY.sort.rmdup.bam' for i in samples[sample]]
    merge_bam_bf_select = out_dir + '/' + sample + '/01_merge_bam/' + sample + '.merge.bf_select.bam'
    
    merge_bam_str = '#/bin/sh\necho "Start merging bam of ' + sample + '!" && \n' \
    + samtools + ' merge -@ 24 ' + merge_bam_bf_select + ' ' + ' '.join(bams) + ' 2>' + merge_bam_log + '  && \n' \
    + bam_selection + ' -i ' +merge_bam_bf_select + ' -p ' + sample + ' -o ' + out_dir + '/' + sample + '/01_merge_bam  -n 100000000 2>> ' + merge_bam_log + ' && \n' \
    + 'echo "Merge bam ' + sample + ' done! \"'
    
    merge_bam_shell_file = open(qsu_dir + '/01_merge_bam/merge_bam.' + sample + '.sh', 'w')
    merge_bam_shell_file.write(merge_bam_str)
    merge_bam_shell_file.close()
    os.system('chmod 755 ' + qsu_dir + '/01_merge_bam/merge_bam.' + sample + '.sh')


''' 2. Bedtools_genomeCoverageBed bedGraph, bedGraphToBigWig bigwig '''

for sample in samples.keys():
    merge_bam = out_dir + '/' + sample + '/01_merge_bam/' + sample + '.merge.bam'
    
    ensure_dir(out_dir + '/' + sample + '/02_Bedtools_bedGraph_bigwig')
    ensure_dir(qsu_dir + '/02_Bedtools_bedGraph_bigwig')
    
    Bedtools_bedGraph_file = out_dir + '/' + sample + '/02_Bedtools_bedGraph_bigwig/' + sample + '.no_chrY.sort.rmdup.bedGraph'
    bedGraphToBigWig_file = out_dir + '/' + sample + '/02_Bedtools_bedGraph_bigwig/' + sample + '.no_chrY.sort.rmdup.bw'
    
    bedGraph_bigwig_str = '#/bin/sh\necho "Start converting the bam file of ' + sample + ' to bedGraph file by Bedtools_genomeCoverageBed !" && \n' \
    + genomeCoverageBed + ' -bg -split -ibam ' + merge_bam + ' -g ' + chromosomesize + ' > ' + Bedtools_bedGraph_file + ' && \n' \
    + 'echo "Make bedGraph of ' + sample + 'has been done!" && \n' \
    + 'echo "bedGraphToBigWig convert bedGraph of ' + sample + '!" && \n' \
    + bedGraphToBigWig + ' ' + Bedtools_bedGraph_file + ' ' + chromosomesize + ' ' + bedGraphToBigWig_file + ' && \n' \
    + 'echo "bedGraphToBigWig make bedGraph of ' + sample + ' has been done!" '
    
    bedGraph_bigwig_shell_file = open(qsu_dir + '/02_Bedtools_bedGraph_bigwig/bedGraph_bigwig.' + sample + '.sh', 'w')
    bedGraph_bigwig_shell_file.write(bedGraph_bigwig_str)
    bedGraph_bigwig_shell_file.close()
    os.system('chmod 755 ' + qsu_dir + '/02_Bedtools_bedGraph_bigwig/bedGraph_bigwig.' + sample + '.sh')



''' 05_damidseq_pipeline '''
for sample in samples.keys():
    merge_bam = out_dir + '/' + sample + '/01_merge_bam/' + sample + '.merge.bam'
    
    ensure_dir(out_dir+'/'+sample+'/05_damidseq_pipeline')
    ensure_dir(qsu_dir+'/05_damidseq_pipeline')

    merge_bam = out_dir + '/' + sample + '/01_merge_bam/' + sample + '.merge.bam'
    control_n="Dam"+re.findall(r"(-.*)",sample)[0]
    control_bam=out_dir + '/' + control_n + '/01_merge_bam/' + control_n + '.merge.bam'

    damidseq_pipeline='/work2/liugh/liuzunpeng/01_Script/01_pipeline/07_DamID_seq/DamID_pipeline-master/damidseq_pipeline'

    damidseq_pipeline_str='#!/bin/sh\necho "Start damidseq_pipeline of '+sample+'!" && \n'\
    +damidseq_pipeline+' --load_defaults=hg19 --threads=24 --bamfiles ' +merge_bam+ ' --dam='+control_bam+' \n'\
    +'echo "Finish damidseq_pipeline of '+sample+' !" '

    damidseq_pipeline_shell_file=open(qsu_dir+'/05_damidseq_pipeline/damidseq.'+sample+'.sh','w')
    damidseq_pipeline_shell_file.write(damidseq_pipeline_str)
    damidseq_pipeline_shell_file.close()
    os.system('chmod 755 ' + qsu_dir+'/05_damidseq_pipeline/damidseq.'+sample+'.sh')



''' 07.  Deeptools bamCoverage and plotCoverage  '''
for sample in samples.keys():
    merge_bam = out_dir + '/' + sample + '/01_merge_bam/' + sample + '.merge.bam'
    
    ensure_dir(out_dir + '/' + sample + '/06_Deeptools_bamCoverage')
    ensure_dir(qsu_dir + '/06_Deeptools_bamCoverage')
    
    ''' 7.1  deeptools_bin_10  '''
    ensure_dir(out_dir + '/' + sample + '/06_Deeptools_bamCoverage/01_bin_10')
    ensure_dir(qsu_dir + '/06_Deeptools_bamCoverage/01_bin_10')
    
    deeptools_bin_10_bedGraph = out_dir + '/' + sample + '/06_Deeptools_bamCoverage/01_bin_10/' + sample + '.no_chrY_unique.sort.rmdup.bin_10.bedGraph'
    deeptools_bin_10_BigWig = out_dir + '/' + sample + '/06_Deeptools_bamCoverage/01_bin_10/' + sample + '.no_chrY_unique.sort.rmdup.bin_10.bw'
    
    bin_10_log = out_dir + '/' + sample + '/06_Deeptools_bamCoverage/01_bin_10/' + sample + '.bin_10.log'
    
    bedGraph_bigwig_bin_10_str = '#/bin/sh\necho "Start deeptools bamCoverage of ' + sample + ' at bin size 10 bp !" && \n' \
    + bamCoverage + ' -p 24  -b ' + merge_bam + ' -o ' + deeptools_bin_10_BigWig + ' --normalizeUsing RPKM --binSize 10 --ignoreForNormalization chrX chrM chrY 2> ' + bin_10_log + ' &&\n' \
    + 'echo "Finish deeptools bamCoverage of ' + sample + ' at bin size 10 bp!" '
    
    bedGraph_bigwig_bin_10_shell_file = open(qsu_dir + '/06_Deeptools_bamCoverage/01_bin_10/bedGraph_bigwig_bin_10.' + sample + '.sh', 'w')
    
    bedGraph_bigwig_bin_10_shell_file.write(bedGraph_bigwig_bin_10_str)
    bedGraph_bigwig_bin_10_shell_file.close()
    os.system('chmod 755 ' + qsu_dir + '/06_Deeptools_bamCoverage/01_bin_10/bedGraph_bigwig_bin_10.' + sample + '.sh')
    
    ''' 7.2  deeptools_bin_100 '''
    ensure_dir(out_dir + '/' + sample + '/06_Deeptools_bamCoverage/02_bin_100')
    ensure_dir(qsu_dir + '/06_Deeptools_bamCoverage/02_bin_100')

    deeptools_bin_100_bedGraph = out_dir + '/' + sample + '/06_Deeptools_bamCoverage/02_bin_100/' + sample + '.no_chrY_unique.sort.rmdup.bin_100.bedGraph'
    deeptools_bin_100_bedGraph_chr = out_dir + '/' + sample + '/06_Deeptools_bamCoverage/02_bin_100/' + sample + '.no_chrY_unique.sort.rmdup.bin_100_chr.bedGraph'
    deeptools_bin_100_bedGraph_chr_zcore = out_dir + '/' + sample + '/06_Deeptools_bamCoverage/02_bin_100/' + sample + '.no_chrY_unique.sort.rmdup.bin_100_chr_Zscore.bedGraph'
    
    deeptools_bin_100_BigWig = out_dir + '/' + sample + '/06_Deeptools_bamCoverage/02_bin_100/' + sample + '.no_chrY_unique.sort.rmdup.bin_100.bw'
    deeptools_bin_100_BigWig_Zscore = out_dir + '/' + sample + '/06_Deeptools_bamCoverage/02_bin_100/' + sample + '.no_chrY_unique.sort.rmdup.bin_100_Zscore.bw'
    
    promoter_5K = out_dir + '/' + sample + '/06_Deeptools_bamCoverage/02_bin_100/' + sample + '.promoter_5K.txt'
    promoter_3K = out_dir + '/' + sample + '/06_Deeptools_bamCoverage/02_bin_100/' + sample + '.promoter_3K.txt'
    promoter_5K_Zscore = out_dir + '/' + sample + '/06_Deeptools_bamCoverage/02_bin_100/' + sample + '.promoter_5K_Zscore.txt'
    promoter_3K_Zscore = out_dir + '/' + sample + '/06_Deeptools_bamCoverage/02_bin_100/' + sample + '.promoter_3K_Zscore.txt'
    
    upstream_25K= out_dir + '/' + sample + '/06_Deeptools_bamCoverage/02_bin_100/' + sample + '.upstream_25K.txt'
    upstream_25K_Zscore= out_dir + '/' + sample + '/06_Deeptools_bamCoverage/02_bin_100/' + sample + '.upstream_25K_Zscore.txt'
    downstream_25K= out_dir + '/' + sample + '/06_Deeptools_bamCoverage/02_bin_100/' + sample + '.downstream_25K.txt'
    downstream_25K_Zscore= out_dir + '/' + sample + '/06_Deeptools_bamCoverage/02_bin_100/' + sample + '.downstream_25K_Zscore.txt'
    
    bin_100_log = out_dir + '/' + sample + '/06_Deeptools_bamCoverage/02_bin_100/' + sample + '.bin_100.log'
    
    bedGraph_bigwig_bin_100_str = '#/bin/sh\necho "Start deeptools bamCoverage of ' + sample + ' at bin size 100 bp !" && \n' \
    + bamCoverage + ' -p 24 -b ' + merge_bam + ' --outFileFormat bedgraph -o ' + deeptools_bin_100_bedGraph + ' --normalizeUsing RPKM --binSize 100 --ignoreForNormalization chrM chrY 2>> ' + bin_100_log + ' &&\n' \
    + bamCoverage + ' -p 24 -b ' + merge_bam + ' -o ' + deeptools_bin_100_BigWig + ' --normalizeUsing RPKM --binSize 100 --ignoreForNormalization chrM chrY 2> ' + bin_100_log + ' && \n' \
    + bedtools + '/bedtools intersect -a /work2/liugh/liuzunpeng/04_Softwares/bedtools2/genomes/human.hg19.genome.w100.bed -b '+deeptools_bin_100_bedGraph+' -wao |awk \'OFS="\\t" {print $1,$2,$3,$10}\'|sort -k1,1 -k2,2n >'+deeptools_bin_100_bedGraph_chr + ' && \n' \
    + "/work2/liugh/liuzunpeng/04_Softwares/python/local/python2/bin/python2 /work2/liugh/liuzunpeng/04_Softwares/source/Zscore.py -b "+deeptools_bin_100_bedGraph_chr +' -o '+ deeptools_bin_100_bedGraph_chr_zcore  + ' && \n' \
    + bedGraphToBigWig + ' ' + deeptools_bin_100_bedGraph_chr_zcore + ' /work2/liugh/liuzunpeng/04_Softwares/bedtools2/genomes/human.hg19.genome '+  deeptools_bin_100_BigWig_Zscore + ' && \n' \
    + bedtools + '/bedtools intersect -a /work2/liugh/liuzunpeng/03_database/00_genome/01_hg19/03_Ensembl/01_fa_gtf_bed_chr/00_bed/gtf_2_bed.20190908/Homo_sapiens.GRCh37.87.final_promoter_5000.bed -b '+ deeptools_bin_100_bedGraph_chr +' -wao |awk \'OFS="\\t" {print $1,$2,$3,$4,$5,$6,$10}\'| '+bedtools + '/bedtools groupby -g 1,2,3,4,5,6  -c 7 -o mean |awk \'OFS="\\t" {print $1,$2,$3,$4,$7,$6}\' >'+promoter_5K + ' && \n' \
    + bedtools + '/bedtools intersect -a /work2/liugh/liuzunpeng/03_database/00_genome/01_hg19/03_Ensembl/01_fa_gtf_bed_chr/00_bed/gtf_2_bed.20190908/Homo_sapiens.GRCh37.87.final_promoter_3000.bed -b '+ deeptools_bin_100_bedGraph_chr +' -wao |awk \'OFS="\\t" {print $1,$2,$3,$4,$5,$6,$10}\'| '+bedtools + '/bedtools groupby -g 1,2,3,4,5,6  -c 7 -o mean |awk \'OFS="\\t" {print $1,$2,$3,$4,$7,$6}\' >'+promoter_3K + ' && \n' \
    + bedtools + '/bedtools intersect -a /work2/liugh/liuzunpeng/03_database/00_genome/01_hg19/03_Ensembl/01_fa_gtf_bed_chr/00_bed/gtf_2_bed.20190908/Homo_sapiens.GRCh37.87.final_promoter_5000.bed -b '+ deeptools_bin_100_bedGraph_chr_zcore +' -wao |awk \'OFS="\\t" {print $1,$2,$3,$4,$5,$6,$10}\'| '+bedtools + '/bedtools groupby -g 1,2,3,4,5,6  -c 7 -o mean |awk \'OFS="\\t" {print $1,$2,$3,$4,$7,$6}\' >'+promoter_5K_Zscore + ' && \n' \
    + bedtools + '/bedtools intersect -a /work2/liugh/liuzunpeng/03_database/00_genome/01_hg19/03_Ensembl/01_fa_gtf_bed_chr/00_bed/gtf_2_bed.20190908/Homo_sapiens.GRCh37.87.final_promoter_3000.bed -b '+ deeptools_bin_100_bedGraph_chr_zcore +' -wao |awk \'OFS="\\t" {print $1,$2,$3,$4,$5,$6,$10}\'| '+bedtools + '/bedtools groupby -g 1,2,3,4,5,6  -c 7 -o mean |awk \'OFS="\\t" {print $1,$2,$3,$4,$7,$6}\' >'+promoter_3K_Zscore + ' && \n' \
    + bedtools + '/bedtools intersect -a /work2/liugh/liuzunpeng/04_Softwares/bedtools2/genomes/hg19.25k.upstream_2M.final.bed -b '+ deeptools_bin_100_bedGraph_chr +' -wao |awk \'OFS="\\t" {print $1,$2,$3,$4,$5,$6,$10}\'| '+bedtools + '/bedtools groupby -g 1,2,3,4,5,6  -c 7 -o mean |awk \'OFS="\\t" {print $1,$2,$3,$4,$7,$6}\' >'+upstream_25K + ' && \n' \
    + bedtools + '/bedtools intersect -a /work2/liugh/liuzunpeng/04_Softwares/bedtools2/genomes/hg19.25k.downstream_2M.final.bed -b '+ deeptools_bin_100_bedGraph_chr +' -wao |awk \'OFS="\\t" {print $1,$2,$3,$4,$5,$6,$10}\'| '+bedtools + '/bedtools groupby -g 1,2,3,4,5,6  -c 7 -o mean |awk \'OFS="\\t" {print $1,$2,$3,$4,$7,$6}\' >'+downstream_25K + ' && \n' \
    + bedtools + '/bedtools intersect -a /work2/liugh/liuzunpeng/04_Softwares/bedtools2/genomes/hg19.25k.upstream_2M.final.bed -b '+ deeptools_bin_100_bedGraph_chr_zcore +' -wao |awk \'OFS="\\t" {print $1,$2,$3,$4,$5,$6,$10}\'| '+bedtools + '/bedtools groupby -g 1,2,3,4,5,6  -c 7 -o mean |awk \'OFS="\\t" {print $1,$2,$3,$4,$7,$6}\' >'+upstream_25K_Zscore + ' && \n' \
    + bedtools + '/bedtools intersect -a /work2/liugh/liuzunpeng/04_Softwares/bedtools2/genomes/hg19.25k.downstream_2M.final.bed -b '+ deeptools_bin_100_bedGraph_chr_zcore +' -wao |awk \'OFS="\\t" {print $1,$2,$3,$4,$5,$6,$10}\'| '+bedtools + '/bedtools groupby -g 1,2,3,4,5,6  -c 7 -o mean |awk \'OFS="\\t" {print $1,$2,$3,$4,$7,$6}\' >'+downstream_25K_Zscore + ' && \n' \
    + 'echo "Finish deeptools bamCoverage of ' + sample + ' at bin size 100 bp!" '
    
    bedGraph_bigwig_bin_100_shell_file = open(qsu_dir + '/06_Deeptools_bamCoverage/02_bin_100/bedGraph_bigwig_bin_100.' + sample + '.sh', 'w')
    bedGraph_bigwig_bin_100_shell_file.write(bedGraph_bigwig_bin_100_str)
    bedGraph_bigwig_bin_100_shell_file.close()
    os.system('chmod 755 ' + qsu_dir + '/06_Deeptools_bamCoverage/02_bin_100/bedGraph_bigwig_bin_100.' + sample + '.sh')
    
    
    ''' 7.3  deeptools_bin_2000  '''
    ensure_dir(out_dir + '/' + sample + '/06_Deeptools_bamCoverage/03_bin_2000')
    ensure_dir(qsu_dir + '/06_Deeptools_bamCoverage/03_bin_2000')
    
    deeptools_bin_2000_bedGraph = out_dir + '/' + sample + '/06_Deeptools_bamCoverage/03_bin_2000/' + sample + '.no_chrY_unique.sort.rmdup.bin_2000.bedGraph'
    deeptools_bin_2000_BigWig = out_dir + '/' + sample + '/06_Deeptools_bamCoverage/03_bin_2000/' + sample + '.no_chrY_unique.sort.rmdup.bin_2000.bw'
    
    bin_2000_log = out_dir + '/' + sample + '/06_Deeptools_bamCoverage/03_bin_2000/' + sample + '.bin_2000.log'
    
    bedGraph_bigwig_bin_2000_str = '#/bin/sh\n\echo "Start deeptools bamCoverage of ' + sample + ' at bin size 2000 bp !" && \n' \
    + bamCoverage + ' -p 24  -b ' + merge_bam + ' --outFileFormat bedgraph -o ' + deeptools_bin_2000_bedGraph + ' --normalizeUsing RPKM --binSize 2000 --ignoreForNormalization chrX chrM chrY 2>> ' + bin_2000_log + ' &&\n' \
    + 'echo "Finish deeptools bamCoverage of ' + sample + ' at bin size 2000 bp!" '
    
    bedGraph_bigwig_bin_2000_shell_file = open(qsu_dir + '/06_Deeptools_bamCoverage/03_bin_2000/bedGraph_bigwig_bin_2000.' + sample + '.sh', 'w')
    bedGraph_bigwig_bin_2000_shell_file.write(bedGraph_bigwig_bin_2000_str)
    bedGraph_bigwig_bin_2000_shell_file.close()
    os.system('chmod 755 ' + qsu_dir + '/06_Deeptools_bamCoverage/03_bin_2000/bedGraph_bigwig_bin_2000.' + sample + '.sh')
    
    ''' 7.4  deeptools_bin_100000  '''
    ensure_dir(out_dir + '/' + sample + '/06_Deeptools_bamCoverage/04_bin_100000')
    ensure_dir(qsu_dir + '/06_Deeptools_bamCoverage/04_bin_100000')
    
    deeptools_bin_100000_bedGraph = out_dir + '/' + sample + '/06_Deeptools_bamCoverage/04_bin_100000/' + sample + '.no_chrY_unique.sort.rmdup.bin_100000.bedGraph'
    deeptools_bin_100000_BigWig = out_dir + '/' + sample + '/06_Deeptools_bamCoverage/04_bin_100000/' + sample + '.no_chrY_unique.sort.rmdup.bin_100000.bw'
    
    bin_100000_log = out_dir + '/' + sample + '/06_Deeptools_bamCoverage/04_bin_100000/' + sample + '.bin_100000.log'
    
    bedGraph_bigwig_bin_100000_str = '#/bin/sh\n\echo "Start deeptools bamCoverage of ' + sample + ' at bin size 100000 bp !" && \n' \
    + bamCoverage + ' -p 24  -b ' + merge_bam + ' --outFileFormat bedgraph -o ' + deeptools_bin_100000_bedGraph + ' --normalizeUsing RPKM --binSize 100000 --ignoreForNormalization chrX chrM chrY 2>> ' + bin_100000_log + ' &&\n' \
    + 'echo "Finish deeptools bamCoverage of ' + sample + ' at bin size 100000 bp!" '
    
    bedGraph_bigwig_bin_100000_shell_file = open(qsu_dir + '/06_Deeptools_bamCoverage/04_bin_100000/bedGraph_bigwig_bin_100000.' + sample + '.sh', 'w')
    bedGraph_bigwig_bin_100000_shell_file.write(bedGraph_bigwig_bin_100000_str)
    bedGraph_bigwig_bin_100000_shell_file.close()
    os.system('chmod 755 ' + qsu_dir + '/06_Deeptools_bamCoverage/04_bin_100000/bedGraph_bigwig_bin_100000.' + sample + '.sh')


###################### normlized to Dam ################

''' 20. bamcompare  '''
for sample in samples.keys():
    ensure_dir(out_dir + '/' + sample + '/07_bamcompare_norm_to_Dam')
    ensure_dir(qsu_dir + '/07_bamcompare_norm_to_Dam')
    
    merge_bam = out_dir + '/' + sample + '/01_merge_bam/' + sample + '.merge.bam'
    control_n="Dam"+re.findall(r"(-.*)",sample)[0]
    control_bam=out_dir + '/' + control_n + '/01_merge_bam/' + control_n + '.merge.bam'

    ''' 20.1  deeptools_bin_10  '''
    ensure_dir(out_dir + '/' + sample + '/07_bamcompare_norm_to_Dam/01_bin_10')
    ensure_dir(qsu_dir + '/07_bamcompare_norm_to_Dam/01_bin_10')
    
    deeptools_bin_10_BigWig = out_dir + '/' + sample + '/07_bamcompare_norm_to_Dam/01_bin_10/' + sample + '.bin_10.norm_Dam.bw'
    
    bin_10_log = out_dir + '/' + sample + '/07_bamcompare_norm_to_Dam/01_bin_10/' + sample + '.bin_10.log'
    
    bamcompare_bin_10_str = '#!/bin/sh\necho "Start deeptools bamCoverage of ' + sample + ' at bin size 10 bp !" && \n' \
    + bamCompare + ' -p 24  --bamfile1 ' + merge_bam   + ' --bamfile2 ' + control_bam + ' --outFileFormat bigwig --outFileName ' + deeptools_bin_10_BigWig + ' --normalizeUsing RPKM --ratio log2  --binSize 10 --ignoreForNormalization chrX chrM chrY 2> ' + bin_10_log + ' &&\n' \
    + 'echo "Finish deeptools bamCoverage of ' + sample + ' at bin size 10 bp!" '
    
    bamcompare_bin_10_shell_file = open(qsu_dir + '/07_bamcompare_norm_to_Dam/01_bin_10/bigwig_bin_10.' + sample + '.sh', 'w')
    bamcompare_bin_10_shell_file.write(bamcompare_bin_10_str)
    bamcompare_bin_10_shell_file.close()
    os.system('chmod 755 ' + qsu_dir + '/07_bamcompare_norm_to_Dam/01_bin_10/bigwig_bin_10.' + sample + '.sh')

    
    
    ''' 20.2  deeptools_bin_100  '''
    ensure_dir(out_dir + '/' + sample + '/07_bamcompare_norm_to_Dam/02_bin_100')
    ensure_dir(qsu_dir + '/07_bamcompare_norm_to_Dam/02_bin_100')
    
    deeptools_bin_100_BigWig = out_dir + '/' + sample + '/07_bamcompare_norm_to_Dam/02_bin_100/' + sample + '.bin_100.norm_Dam.bw'
    deeptools_bin_100_bedgraph = out_dir + '/' + sample + '/07_bamcompare_norm_to_Dam/02_bin_100/' + sample + '.bin_100.norm_Dam.bedgraph'

    bin_100_log = out_dir + '/' + sample + '/07_bamcompare_norm_to_Dam/02_bin_100/' + sample + '.bin_100.log'
    
    bamcompare_bin_100_str = '#!/bin/sh\necho "Start deeptools bamCoverage of ' + sample + ' at bin size 100 bp !" && \n' \
    + bamCompare + ' -p 24  --bamfile1 ' + merge_bam   + ' --bamfile2 ' + control_bam + ' --outFileFormat bigwig --outFileName ' + deeptools_bin_100_BigWig + ' --normalizeUsing RPKM --ratio log2 --binSize 100 --ignoreForNormalization chrM chrY 2> ' + bin_100_log + ' &&\n' \
    + bamCompare + ' -p 24  --bamfile1 ' + merge_bam   + ' --bamfile2 ' + control_bam + ' --outFileFormat bedgraph --outFileName ' + deeptools_bin_100_bedgraph + ' --normalizeUsing RPKM --ratio log2 --binSize 100 --ignoreForNormalization chrM chrY 2> ' + bin_100_log + ' &&\n' \
    + 'echo "Finish deeptools bamCoverage of ' + sample + ' at bin size 100 bp!" '
    
    bamcompare_bin_100_shell_file = open(qsu_dir + '/07_bamcompare_norm_to_Dam/02_bin_100/bigwig_bin_100.' + sample + '.sh', 'w')
    bamcompare_bin_100_shell_file.write(bamcompare_bin_100_str)
    bamcompare_bin_100_shell_file.close()
    os.system('chmod 755 ' + qsu_dir + '/07_bamcompare_norm_to_Dam/02_bin_100/bigwig_bin_100.' + sample + '.sh')

    ''' 20.3  deeptools_bin_2000  '''
    ensure_dir(out_dir + '/' + sample + '/07_bamcompare_norm_to_Dam/03_bin_2000')
    ensure_dir(qsu_dir + '/07_bamcompare_norm_to_Dam/03_bin_2000')
    
    deeptools_bin_2000_bedGraph = out_dir + '/' + sample + '/07_bamcompare_norm_to_Dam/03_bin_2000/' + sample + '.bin_500.norm_Dam.bw'
    
    bin_2000_log = out_dir + '/' + sample + '/07_bamcompare_norm_to_Dam/03_bin_2000/' + sample + '.bin_10.log'
    
    bamcompare_bin_2000_str = '#!/bin/sh\necho "Start deeptools bamCoverage of ' + sample + ' at bin size 2000 bp !" && \n' \
    + bamCompare + ' -p 24  --bamfile1 ' + merge_bam   + ' --bamfile2 ' + control_bam + ' --scaleFactorsMethod None  --outFileFormat bigwig  --outFileName ' + deeptools_bin_2000_bedGraph + ' --normalizeUsing RPKM --operation log2 --binSize 500 --ignoreForNormalization chrM chrY 2> ' + bin_2000_log + ' &&\n' \
    + 'echo "Finish deeptools bamCoverage of ' + sample + ' at bin size 2000 bp!" '
    
    bamcompare_bin_2000_shell_file = open(qsu_dir + '/07_bamcompare_norm_to_Dam/03_bin_2000/bedGraph_bin_2000.' + sample + '.sh', 'w')
    bamcompare_bin_2000_shell_file.write(bamcompare_bin_2000_str)
    bamcompare_bin_2000_shell_file.close()
    os.system('chmod 755 ' + qsu_dir + '/07_bamcompare_norm_to_Dam/03_bin_2000/bedGraph_bin_2000.' + sample + '.sh')
    
    ''' 20.4_bin_100000  '''
    ensure_dir(out_dir + '/' + sample + '/07_bamcompare_norm_to_Dam/04_bin_100000')
    ensure_dir(qsu_dir + '/07_bamcompare_norm_to_Dam/04_bin_100000')
    
    deeptools_bin_100000_bedGraph = out_dir + '/' + sample + '/07_bamcompare_norm_to_Dam/04_bin_100000/' + sample + '.bin_100000.norm_Dam.bed'
    
    bin_100000_log = out_dir + '/' + sample + '/07_bamcompare_norm_to_Dam/04_bin_100000/' + sample + '.bin_100000.log'
    
    bamcompare_bin_100000_str = '#!/bin/sh\necho "Start deeptools bamCoverage of ' + sample + ' at bin size 100000 bp !" && \n' \
    + bamCompare + ' -p 24  --bamfile1 ' + merge_bam   + ' --bamfile2 ' + control_bam + ' --outFileFormat bedgraph --outFileName ' + deeptools_bin_100000_bedGraph + ' --normalizeUsing RPKM --ratio log2 --binSize 100000 --ignoreForNormalization chrX chrM chrY 2> ' + bin_100000_log + ' &&\n' \
    + 'echo "Finish deeptools bamCoverage of ' + sample + ' at bin size 100000 bp!" '
    
    bamcompare_bin_100000_shell_file = open(qsu_dir + '/07_bamcompare_norm_to_Dam/04_bin_100000/bedGraph_bin_100000.' + sample + '.sh', 'w')
    bamcompare_bin_100000_shell_file.write(bamcompare_bin_100000_str)
    bamcompare_bin_100000_shell_file.close()
    os.system('chmod 755 ' + qsu_dir + '/07_bamcompare_norm_to_Dam/04_bin_100000/bedGraph_bin_100000.' + sample + '.sh')



