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
    
    python3 ChIP_seq.part2.py \
    -i /work1/liugh/liuzunpeng/05_Results/05_Lamin_KI_H3K9me3 \
    -o /work1/liugh/liuzunpeng/05_Results/05_Lamin_KI_H3K9me3_part2 \
    -q /work1/liugh/liuzunpeng/06_Qsub/05_Lamin_KI_H3K9me3_part2 \
    -c /work2/liugh/liuzunpeng/01_Script/01_pipeline/02_ChIP-seq/00_ChIP_seq_pipeline/config.txt \
    -g H3,H9,M3,M9,R3,R9,W3,W9 \
    -p 05_Lamin_KI_H3K9me3 \
    -b /work2/liugh/liuzunpeng/01_Script/01_pipeline/02_ChIP-seq/00_ChIP_seq_pipeline/bed.txt
    
    
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
    + samtools + ' merge -@ 20 ' + merge_bam_bf_select + ' ' + ' '.join(bams) + ' 2>' + merge_bam_log + '  && \n' \
    + bam_selection + ' -i ' +merge_bam_bf_select + ' -p ' + sample + ' -o ' + out_dir + '/' + sample + '/01_merge_bam  -n 210000000 2>> ' + merge_bam_log + ' && \n' \
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

''' 3. Calling-peak-by-MACS2 '''
for sample in samples.keys():
    merge_bam = out_dir + '/' + sample + '/01_merge_bam/' + sample + '.merge.bam'
    macs_control_bam='/work2/liugh/liuzunpeng/05_Results/06_Rloop_cell_type_input/WT-MSC/05_unique_mapped_reads/WT-MSC_no_chrY.sort.rmdup.bam'
    
    ensure_dir(out_dir+'/'+sample+'/03_MACS2_calling_peak')
    ensure_dir(qsu_dir+'/03_MACS2_calling_peak')
    
    peaks_xls=out_dir+'/'+sample+'/03_MACS2_calling_peak/'+sample+'_peaks.xls'
    peak_bed=out_dir+'/'+sample+'/03_MACS2_calling_peak/'+sample+'.macs2.peaks.bed'
    peak_enrich5_bed= out_dir+'/'+sample+'/03_MACS2_calling_peak/'+sample+'.macs2.peaks_enrich5.bed'
    
    macs2_enrich_txt= out_dir + '/' + sample + '/03_MACS2_calling_peak/' + sample + '.macs2.peaks.enrich.element.txt'
    macs2_enrich_p1 = out_dir + '/' + sample + '/03_MACS2_calling_peak/' + sample + '.macs2.peaks.enrich.element_p1.txt'
    macs2_enrich_p2 = out_dir + '/' + sample + '/03_MACS2_calling_peak/' + sample + '.macs2.peaks.enrich.element_p2.txt'
    
    MACS2_log=out_dir+'/'+sample+'/03_MACS2_calling_peak/'+sample+'.MACS2.log'
    peak_stat=out_dir+'/'+sample+'/03_MACS2_calling_peak/'+sample+'.MACS2_peak_stat.txt'
    peak_enrich5_stat=out_dir+'/'+sample+'/03_MACS2_calling_peak/'+sample+'.MACS2_peak_enrich5_stat.txt'
    
    MACS2_str='#/bin/sh\necho "Start calling peak of '+sample+' by MACS2!" && \n'\
    + macs2 + ' callpeak -t ' + merge_bam + ' -c ' + macs_control_bam + ' -f BAM -g hs -B --cutoff-analysis -q 0.01 -n ' + sample + ' --outdir ' + out_dir + '/' + sample + '/03_MACS2_calling_peak 2> ' + MACS2_log + ' && \n' \
    +'echo "Finish calling peak of '+sample+' by MACS2!" && \n'\
    + macs2_peak_reshape + ' -i ' + peaks_xls + ' -o ' + out_dir + '/' + sample + '/03_MACS2_calling_peak/ -p ' + sample + ' -s + ' + ' && \n' \
    + 'head -n 12 ' + macs2_enrich_txt + ' |grep -E "3UTR|TTS|Exon|Intron|Intergenic|Promoter|5UTR|Ann" |sed \'1c Annotation\\tPeak_num\\tTotal_size\\tLog2_Enrichment\' >' +macs2_enrich_p1 + ' && \n' \
    + 'tail -n 34 ' + macs2_enrich_txt + ' |grep -E "SINE|LINE|LTR|Satellite|rRNA|RC|Low_complexity|Simple_repeat|Other|DNA|RNA|snoRNA|AnnncRNA|srpRNA|tRNA|snRNA|scRNA|CpG-Island|Unknown" |grep -v ? |sed \'1c Annotation\\tPeak_num\\tTotal_size\\tLog2_Enrichment\' >' +macs2_enrich_p2 + ' && \n' \
    +'echo "Calling peak and evaluate peak state of '+sample+' has been done!" '
    
    MACS2_shell_file=open(qsu_dir+'/03_MACS2_calling_peak/MACS2.'+sample+'.sh','w')
    MACS2_shell_file.write(MACS2_str)
    MACS2_shell_file.close()
    os.system('chmod 755 ' + qsu_dir+'/03_MACS2_calling_peak/MACS2.'+sample+'.sh')



''' 4. Calling-peak-by-SICER '''
for sample in samples.keys():
    merge_bam = out_dir + '/' + sample + '/01_merge_bam/' + sample + '.merge.bam'
    ensure_dir(out_dir + '/' + sample + '/04_SICER_calling_peak')
    ensure_dir(qsu_dir + '/04_SICER_calling_peak')
    
    control_bam='/work2/liugh/liuzunpeng/05_Results/06_Rloop_cell_type_input/WT-MSC/05_unique_mapped_reads/WT-MSC_no_chrY.sort.rmdup.bam'
    SICER_log = out_dir + '/' + sample + '/04_SICER_calling_peak/' + sample + '.SICER.log'
    SICER_bed = out_dir + '/' + sample + '/04_SICER_calling_peak/' + sample + '.SICER.bed'
    SICER_FDR_001_bed = out_dir + '/' + sample + '/04_SICER_calling_peak/' + sample + '.SICER_FDR_001.bed'
    SICER_FDR_001_reshape=out_dir + '/' + sample + '/04_SICER_calling_peak/' + sample + '.SICER.peaks.bed'
    SICER_enrich_txt= out_dir + '/' + sample + '/04_SICER_calling_peak/' + sample + '.SICER.peaks.enrich.element.txt'
    SICER_enrich_p1 = out_dir + '/' + sample + '/04_SICER_calling_peak/' + sample + '.SICER.peaks.enrich.element_p1.txt'
    SICER_enrich_p2 = out_dir + '/' + sample + '/04_SICER_calling_peak/' + sample + '.SICER.peaks.enrich.element_p2.txt'
    SICER_peaks_d12000 = out_dir + '/' + sample + '/04_SICER_calling_peak/' + sample + '.SICER.peaks.d12000.bed'
    SICER_peaks_d12000_cov = out_dir + '/' + sample + '/04_SICER_calling_peak/' + sample + '.SICER.peaks.d12000.cov.txt'
    
    peak_stat = out_dir + '/' + sample + '/04_SICER_calling_peak/' + sample + '.SICER_peak_stat.txt'
    
    SICER_str = '#/bin/sh\necho "Start calling peak of ' + sample + ' by SICER!" && \n' \
    + SICER + '  -t ' + merge_bam + ' -c ' + control_bam + ' -w 500 -g 5  > ' + SICER_bed + ' 2> ' + SICER_log + ' && \n' \
    + 'awk \'$8 < 0.01 \' ' + SICER_bed + ' > ' + SICER_FDR_001_bed + ' && \n' \
    + 'echo "Finish calling peak of ' + sample + ' by SICER!" && \n' \
    + SICER_peak_reshape + ' -i ' + SICER_FDR_001_bed + ' -o ' + out_dir + '/' + sample + '/04_SICER_calling_peak/ -p ' + sample + ' -s + ' + ' && \n' \
    + 'head -n 12 ' + SICER_enrich_txt + ' |grep -E "3UTR|TTS|Exon|Intron|Intergenic|Promoter|5UTR|Ann" |sed \'1c Annotation\\tPeak_num\\tTotal_size\\tLog2_Enrichment\' >' +SICER_enrich_p1 + ' && \n' \
    + 'tail -n 34 ' + SICER_enrich_txt + ' |grep -E "SINE|LINE|LTR|Satellite|rRNA|RC|Low_complexity|Simple_repeat|Other|DNA|RNA|snoRNA|AnnncRNA|srpRNA|tRNA|snRNA|scRNA|CpG-Island|Unknown" |grep -v ? |sed \'1c Annotation\\tPeak_num\\tTotal_size\\tLog2_Enrichment\' >' +SICER_enrich_p2 + ' && \n' \
    + bedtools + '/bedtools merge -d 12000 -i ' + SICER_FDR_001_reshape + ' > ' +SICER_peaks_d12000 + ' && \n' \
    + bedtools + '/bedtools multicov -bams ' + merge_bam +' -bed ' + SICER_peaks_d12000 + ' > ' + SICER_peaks_d12000_cov + ' && \n' \
    + 'echo "Calling peak and evaluate peak state of ' + sample + ' has been done!" '
    
    SICER_shell_file = open(qsu_dir + '/04_SICER_calling_peak/SICER.' + sample + '.sh', 'w')
    SICER_shell_file.write(SICER_str)
    SICER_shell_file.close()
    os.system('chmod 755 ' + qsu_dir + '/04_SICER_calling_peak/SICER.' + sample + '.sh')



''' 5.Tag_dir  '''
for sample in samples.keys():
    ensure_dir(out_dir + '/' + sample + '/05_Tag_dir')
    ensure_dir(qsu_dir + '/05_Tag_dir')
    merge_bam = out_dir + '/' + sample + '/01_merge_bam/' + sample + '.merge.bam'
    out_prefix=out_dir + '/' + sample + '/05_Tag_dir/'+sample
    
    all_tag_log   = out_dir + '/' + sample + '/05_Tag_dir/' + sample + '.all_tag_log.txt'
    
    Tag_dir_str='#!/bin/sh\necho "makeTagDirectory' + sample + '!" && \n' \
+ makeTagDirectory + ' ' + out_dir + '/' + sample + '/05_Tag_dir -format sam ' + merge_bam + ' 2> ' + all_tag_log +' && \n' \
    + 'echo "makeTagDirectory of ' + sample + ' ALL done!" '
    
    Tag_dir_shell_file = open(qsu_dir + '/05_Tag_dir/Tag.' + sample + '.sh', 'w')
    Tag_dir_shell_file.write(Tag_dir_str)
    Tag_dir_shell_file.close()
    os.system('chmod 755 ' + qsu_dir + '/05_Tag_dir/Tag.' + sample + '.sh')

''' 6.Calling-peak-by-homer '''
for sample in samples.keys():
    ensure_dir(out_dir + '/' + sample + '/06_homer_peak')
    ensure_dir(qsu_dir + '/06_homer_peak')
    control_tag='/work2/liugh/liuzunpeng/05_Results/06_Rloop_cell_type_input/WT-MSC/12_Tag_dir'
    
    homer_bed = out_dir + '/' + sample + '/06_homer_peak/' + sample + '.homer.txt'
    homer_enrich_txt= out_dir + '/' + sample + '/06_homer_peak/' + sample + '.homer.peaks.anno.enrich.txt'
    homer_enrich_p1 = out_dir + '/' + sample + '/06_homer_peak/' + sample + '.homer.peaks.anno.enrich_p1.txt'
    homer_enrich_p2 = out_dir + '/' + sample + '/06_homer_peak/' + sample + '.homer.peaks.anno.enrich_p2.txt'
    
    homer_peak_log   = out_dir + '/' + sample + '/05_Tag_dir/' + sample + '.all_tag_log.txt'
    
    homer_str='#!/bin/sh\necho "homer call peaks of ' + sample + '!" && \n'\
+ findPeaks + ' ' + out_dir + '/' + sample + '/05_Tag_dir -i ' + control_tag + ' -style histone -size 1000 -minDist 1000 > ' + homer_bed +' 2> ' + homer_peak_log +' && \n' \
    + 'echo "Finish calling peak of ' + sample + ' by Homer!" && \n' \
    + homer_peak_reshape_and_anno + ' -i ' + homer_bed + ' -o ' + out_dir + '/' + sample + '/06_homer_peak/ -p ' + sample + ' -s + ' + ' && \n' \
    + 'head -n 12 ' + homer_enrich_txt + ' |grep -E "3UTR|TTS|Exon|Intron|Intergenic|Promoter|5UTR|Ann" |sed \'1c Annotation\\tPeak_num\\tTotal_size\\tLog2_Enrichment\' >' +homer_enrich_p1 + ' && \n' \
    + 'tail -n 34 ' + homer_enrich_txt + ' |grep -E "SINE|LINE|LTR|Satellite|rRNA|RC|Low_complexity|Simple_repeat|Other|DNA|RNA|snoRNA|AnnncRNA|srpRNA|tRNA|snRNA|scRNA|CpG-Island|Unknown" |grep -v ? |sed \'1c Annotation\\tPeak_num\\tTotal_size\\tLog2_Enrichment\' >' +homer_enrich_p2 + ' && \n' \
    + 'echo "Calling peak and evaluate peak state of ' + sample + ' has been done!" '
    
    homer_shell_file = open(qsu_dir + '/06_homer_peak/homer_' + sample + '.sh', 'w')
    homer_shell_file.write(homer_str)
    homer_shell_file.close()
    os.system('chmod 755 ' + qsu_dir + '/06_homer_peak/homer_' + sample + '.sh')


''' 07.  Deeptools bamCoverage and plotCoverage  '''
for sample in samples.keys():
    merge_bam = out_dir + '/' + sample + '/01_merge_bam/' + sample + '.merge.bam'
    
    ensure_dir(out_dir + '/' + sample + '/07_Deeptools_bamCoverage')
    ensure_dir(qsu_dir + '/07_Deeptools_bamCoverage')
    
    ''' 7.1  deeptools_bin_10  '''
    ensure_dir(out_dir + '/' + sample + '/07_Deeptools_bamCoverage/01_bin_10')
    ensure_dir(qsu_dir + '/07_Deeptools_bamCoverage/01_bin_10')
    
    deeptools_bin_10_bedGraph = out_dir + '/' + sample + '/07_Deeptools_bamCoverage/01_bin_10/' + sample + '.no_chrY_unique.sort.rmdup.bin_10.bedGraph'
    deeptools_bin_10_BigWig = out_dir + '/' + sample + '/07_Deeptools_bamCoverage/01_bin_10/' + sample + '.no_chrY_unique.sort.rmdup.bin_10.bw'
    
    bin_10_log = out_dir + '/' + sample + '/07_Deeptools_bamCoverage/01_bin_10/' + sample + '.bin_10.log'
    
    bedGraph_bigwig_bin_10_str = '#/bin/sh\necho "Start deeptools bamCoverage of ' + sample + ' at bin size 10 bp !" && \n' \
    + bamCoverage + ' -p 20  -b ' + merge_bam + ' -o ' + deeptools_bin_10_BigWig + ' --normalizeUsingRPKM --binSize 10 --ignoreForNormalization chrX chrM chrY 2> ' + bin_10_log + ' &&\n' \
    + 'echo "Finish deeptools bamCoverage of ' + sample + ' at bin size 10 bp!" '
    
    bedGraph_bigwig_bin_10_shell_file = open(qsu_dir + '/07_Deeptools_bamCoverage/01_bin_10/bedGraph_bigwig_bin_10.' + sample + '.sh', 'w')
    
    bedGraph_bigwig_bin_10_shell_file.write(bedGraph_bigwig_bin_10_str)
    bedGraph_bigwig_bin_10_shell_file.close()
    os.system('chmod 755 ' + qsu_dir + '/07_Deeptools_bamCoverage/01_bin_10/bedGraph_bigwig_bin_10.' + sample + '.sh')
    
    ''' 7.2  deeptools_bin_100 '''
    ensure_dir(out_dir + '/' + sample + '/07_Deeptools_bamCoverage/02_bin_100')
    ensure_dir(qsu_dir + '/07_Deeptools_bamCoverage/02_bin_100')
    
    deeptools_bin_100_bedGraph = out_dir + '/' + sample + '/07_Deeptools_bamCoverage/02_bin_100/' + sample + '.no_chrY_unique.sort.rmdup.bin_100.bedGraph'
    deeptools_bin_100_bedGraph_chr = out_dir + '/' + sample + '/07_Deeptools_bamCoverage/02_bin_100/' + sample + '.no_chrY_unique.sort.rmdup.bin_100_chr.bedGraph'
    deeptools_bin_100_bedGraph_chr_zcore = out_dir + '/' + sample + '/07_Deeptools_bamCoverage/02_bin_100/' + sample + '.no_chrY_unique.sort.rmdup.bin_100_chr_Zscore.bedGraph'
    
    deeptools_bin_100_BigWig = out_dir + '/' + sample + '/07_Deeptools_bamCoverage/02_bin_100/' + sample + '.no_chrY_unique.sort.rmdup.bin_100.bw'
    deeptools_bin_100_BigWig_Zscore = out_dir + '/' + sample + '/07_Deeptools_bamCoverage/02_bin_100/' + sample + '.no_chrY_unique.sort.rmdup.bin_100_Zscore.bw'
    
    promoter_5K = out_dir + '/' + sample + '/07_Deeptools_bamCoverage/02_bin_100/' + sample + '.promoter_5K.txt'
    promoter_3K = out_dir + '/' + sample + '/07_Deeptools_bamCoverage/02_bin_100/' + sample + '.promoter_3K.txt'
    promoter_5K_Zscore = out_dir + '/' + sample + '/07_Deeptools_bamCoverage/02_bin_100/' + sample + '.promoter_5K_Zscore.txt'
    promoter_3K_Zscore = out_dir + '/' + sample + '/07_Deeptools_bamCoverage/02_bin_100/' + sample + '.promoter_3K_Zscore.txt'
    
    upstream_25K= out_dir + '/' + sample + '/07_Deeptools_bamCoverage/02_bin_100/' + sample + '.upstream_25K.txt'
    upstream_25K_Zscore= out_dir + '/' + sample + '/07_Deeptools_bamCoverage/02_bin_100/' + sample + '.upstream_25K_Zscore.txt'
    downstream_25K= out_dir + '/' + sample + '/07_Deeptools_bamCoverage/02_bin_100/' + sample + '.downstream_25K.txt'
    downstream_25K_Zscore= out_dir + '/' + sample + '/07_Deeptools_bamCoverage/02_bin_100/' + sample + '.downstream_25K_Zscore.txt'

    bin_100_log = out_dir + '/' + sample + '/07_Deeptools_bamCoverage/02_bin_100/' + sample + '.bin_100.log'
    
    bedGraph_bigwig_bin_100_str = '#/bin/sh\necho "Start deeptools bamCoverage of ' + sample + ' at bin size 100 bp !" && \n' \
    + bedtools + '/bedtools intersect -a /work2/liugh/liuzunpeng/04_Softwares/bedtools2/genomes/hg19.25k.downstream_2M.final.bed -b '+ deeptools_bin_100_bedGraph_chr +' -wao |awk \'OFS="\\t" {print $1,$2,$3,$4,$5,$6,$10}\'| '+bedtools + '/bedtools groupby -g 1,2,3,4,5,6  -c 7 -o mean |awk \'OFS="\\t" {print $1,$2,$3,$4,$7,$6}\' >'+downstream_25K + ' & \n' \
    + bedtools + '/bedtools intersect -a /work2/liugh/liuzunpeng/04_Softwares/bedtools2/genomes/hg19.25k.downstream_2M.final.bed -b '+ deeptools_bin_100_bedGraph_chr_zcore +' -wao |awk \'OFS="\\t" {print $1,$2,$3,$4,$5,$6,$10}\'| '+bedtools + '/bedtools groupby -g 1,2,3,4,5,6  -c 7 -o mean |awk \'OFS="\\t" {print $1,$2,$3,$4,$7,$6}\' >'+downstream_25K_Zscore + '  \n' \
    + 'echo "Finish deeptools bamCoverage of ' + sample + ' at bin size 100 bp!" '
    
    bedGraph_bigwig_bin_100_shell_file = open(qsu_dir + '/07_Deeptools_bamCoverage/02_bin_100/bedGraph_bigwig_bin_100.' + sample + '.sh', 'w')
    bedGraph_bigwig_bin_100_shell_file.write(bedGraph_bigwig_bin_100_str)
    bedGraph_bigwig_bin_100_shell_file.close()
    os.system('chmod 755 ' + qsu_dir + '/07_Deeptools_bamCoverage/02_bin_100/bedGraph_bigwig_bin_100.' + sample + '.sh')
    
    
