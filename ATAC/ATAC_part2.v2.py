#!/usr/bin/env python3
# coding: utf-8

from optparse import OptionParser
import os
import sys
import re

'''
    Title:        ATAC-seq pipeline
    Description:  This pipeline is used for ChIP-seq analysis and downstream analysis.
    
    Author:       Zunpeng Liu
    Email:        zunpengliuAT163.com
    Date:         1/6/2018
    Version:      1.1
    Python:       based on python 3.6
    Citation:     please
    '''

MY_USAGE = '''
    
    python3 ATAC_part2.v2.py \
    -i /work1/liugh/liuzunpeng/05_Results/11_Lamin_KI_ATAC \
    -o /work1/liugh/liuzunpeng/05_Results/11_Lamin_KI_ATAC_part2 \
    -q /work1/liugh/liuzunpeng/06_Qsub/11_Lamin_KI_ATAC_part2 \
    -c /work2/liugh/liuzunpeng/01_Script/01_pipeline/02_ChIP-seq/00_ChIP_seq_pipeline/config.txt \
    -g H3,H9,M3,M9,R3,R9,W3,W9 \
    -p Lamin_KI_ATAC \
    -b /work2/liugh/liuzunpeng/01_Script/01_pipeline/02_ChIP-seq/00_ChIP_seq_pipeline/bed.txt
    
    
    Author:       Zunpeng Liu
    Email:        zunpengliu@163.com
    Date:         1/6/2018
    
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
in_dir = str(options.in_dir)
group  = str(options.group)
out_dir = options.out_dir
qsu_dir = str(options.qsub_dir)
log_dir = '/home/liuzunpeng/../..' + qsu_dir

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

all_bam_lst=[]

for sample in samples.keys():
    ensure_dir(out_dir + '/' + sample + '/01_merge_bam')
    merge_bam_log = out_dir + '/' + sample + '/01_merge_bam/' + sample + '.bmerge.bam.log.txt'
    merge_bam = out_dir + '/' + sample + '/01_merge_bam/' + sample + '.merge.bam'
    bams=[in_dir + '/' + i + '/05_unique_mapped_reads/' + i + '_no_chrM_chrY.sort.rmdup.bam' for i in samples[sample]]
    merge_bam_bf_select = out_dir + '/' + sample + '/01_merge_bam/' + sample + '.merge.bf_select.bam'
    
    for i in samples[sample]:
        all_bam_lst.append(in_dir + '/' + i +'/05_unique_mapped_reads/' + i + '_no_chrM_chrY.sort.rmdup.bam')
    print(bams)

    merge_bam_str = '#!/bin/sh\necho "Start merging bam of ' + sample + '!" && \n' \
    + samtools + ' merge -@ 20 ' + merge_bam_bf_select + ' ' + ' '.join(bams) + ' 2>' + merge_bam_log + '  && \n' \
    + bam_selection + ' -i ' +merge_bam_bf_select + ' -p ' + sample + ' -o ' + out_dir + '/' + sample + '/01_merge_bam  -n 24600000 && \n' \
    + 'echo "Merge bam ' + sample + ' done! \"'
    
    merge_bam_shell_file = open(qsu_dir + '/01_merge_bam/merge_bam.' + sample + '.sh', 'w')
    merge_bam_shell_file.write(merge_bam_str)
    merge_bam_shell_file.close()
    os.system('chmod 755 ' + qsu_dir + '/01_merge_bam/merge_bam.' + sample + '.sh')


''' 7. Bam to bed and remove blacklisted-intervals '''
for sample in samples.keys():
    merge_bam = out_dir + '/' + sample + '/01_merge_bam/' + sample + '.merge.bam'
    
    ensure_dir(out_dir+'/'+sample+'/07_Bam2bed_remove_black_list')
    ensure_dir(qsu_dir+'/07_Bam2bed_remove_black_list')

    srt_rmdup_bed=out_dir+'/'+sample+'/07_Bam2bed_remove_black_list/'+sample+'_merge.bed'
    intersect_BL_bed=out_dir+'/'+sample+'/07_Bam2bed_remove_black_list/'+sample+'_merge.intersect_BL.bed'
    remove_BL_bed=out_dir+'/'+sample+'/07_Bam2bed_remove_black_list/'+sample+'_merge.remove_BL.bed'
    bam2bed_log=out_dir+'/'+sample+'/07_Bam2bed_remove_black_list/'+sample+'.merge.bam2bed.log.txt'

    Bam2bed_remove_black_list_str='#!/bin/sh\necho "Start converting the bam file of '+sample+' to bed file!" && \n'\
    +bamToBed + ' -i ' + merge_bam + ' > ' + srt_rmdup_bed +' && \n'\
    +'echo "Start remove black list of '+sample+'! " && \n'\
    +bedtools+'/bedtools intersect -u -a ' + srt_rmdup_bed + ' -b ' + black_list + ' > ' + intersect_BL_bed +' && \n'\
    +bedtools+'/bedtools intersect -v -a ' + srt_rmdup_bed + ' -b ' + black_list + ' > ' + remove_BL_bed +' && \n'\
    +'echo "Convert the bam file of '+sample+' to bed file and remove black list has been done!" '

    Bam2bed_shell_file=open(qsu_dir+'/07_Bam2bed_remove_black_list/Bam2bed_remove_black_list.'+sample+'.sh','w')
    Bam2bed_shell_file.write(Bam2bed_remove_black_list_str)
    Bam2bed_shell_file.close()
    os.system('chmod 755 ' + qsu_dir+'/07_Bam2bed_remove_black_list/Bam2bed_remove_black_list.'+sample+'.sh')

''' 8. Bedtools_genomeCoverageBed bedGraph, bedGraphToBigWig bigwig '''

for sample in samples.keys():
    merge_bam = out_dir + '/' + sample + '/01_merge_bam/' + sample + '.merge.bam'
    ensure_dir(out_dir+'/'+sample+'/08_Bedtools_bedGraph_bigwig')
    ensure_dir(qsu_dir+'/08_Bedtools_bedGraph_bigwig')

    Bedtools_bedGraph_file=out_dir+'/'+sample+'/08_Bedtools_bedGraph_bigwig/'+sample+'.merge.bedGraph'
    bedGraphToBigWig_file =out_dir+'/'+sample+'/08_Bedtools_bedGraph_bigwig/'+sample+'merge.bw'

    bedGraph_bigwig_str='#!/bin/sh\necho "Start converting the bam file of '+sample+' to bedGraph file by Bedtools_genomeCoverageBed !" && \n'\
    +genomeCoverageBed + ' -bg -split -ibam ' + merge_bam + ' -g ' + chromosomesize+' > '+Bedtools_bedGraph_file+' && \n'\
    +'echo "bedGraphToBigWig convert bedGraph of ' + sample+'!" && \n'\
    +bedGraphToBigWig + ' '+Bedtools_bedGraph_file + ' ' + chromosomesize + ' ' + bedGraphToBigWig_file + ' && \n'\
    +'echo "bedGraphToBigWig make bedGraph of '+sample+' has been done!" '

    bedGraph_bigwig_shell_file=open(qsu_dir+'/08_Bedtools_bedGraph_bigwig/bedGraph_bigwig.'+sample+'.sh','w')
    bedGraph_bigwig_shell_file.write(bedGraph_bigwig_str)
    bedGraph_bigwig_shell_file.close()
    os.system('chmod 755 ' + qsu_dir+'/08_Bedtools_bedGraph_bigwig/bedGraph_bigwig.'+sample+'.sh')

''' 9. Calling-peak-by-MACS2 '''
for sample in samples.keys():
    remove_BL_bed=out_dir+'/'+sample+'/07_Bam2bed_remove_black_list/'+sample+'_merge.remove_BL.bed'
    ensure_dir(out_dir+'/'+sample+'/09_MACS2_calling_peak')
    ensure_dir(qsu_dir+'/09_MACS2_calling_peak')

    peaks_xls=out_dir+'/'+sample+'/09_MACS2_calling_peak/'+sample+'_peaks.xls'
    peak_bed=out_dir+'/'+sample+'/09_MACS2_calling_peak/'+sample+'.macs2.peaks.bed'
    peak_enrich5_bed= out_dir+'/'+sample+'/09_MACS2_calling_peak/'+sample+'.macs2.peaks_enrich5.bed'
    
    macs2_enrich_txt= out_dir + '/' + sample + '/09_MACS2_calling_peak/' + sample + '.macs2.peaks.enrich.element.txt'
    macs2_enrich_p1 = out_dir + '/' + sample + '/09_MACS2_calling_peak/' + sample + '.macs2.peaks.enrich.element_p1.txt'
    macs2_enrich_p2 = out_dir + '/' + sample + '/09_MACS2_calling_peak/' + sample + '.macs2.peaks.enrich.element_p2.txt'
    
    MACS2_log=out_dir+'/'+sample+'/09_MACS2_calling_peak/'+sample+'.MACS2.log'
    peak_stat=out_dir+'/'+sample+'/09_MACS2_calling_peak/'+sample+'.MACS2_peak_stat.txt'
    peak_enrich5_stat=out_dir+'/'+sample+'/09_MACS2_calling_peak/'+sample+'.MACS2_peak_enrich5_stat.txt'

    MACS2_str='#!/bin/sh\necho "Start calling peak of '+sample+' by MACS2!" && \n'\
    +macs2+' callpeak -t ' +remove_BL_bed +' -f BED -g hs -B  --shift 100 --extsize 200 -q 0.01 -n '+ sample + ' --outdir '+out_dir+'/'+sample+'/09_MACS2_calling_peak 2> '+MACS2_log+' && \n'\
    +'echo "Finish calling peak of '+sample+' by MACS2!" && \n'\
    + macs2_peak_reshape + ' -i ' + peaks_xls + ' -o ' + out_dir + '/' + sample + '/09_MACS2_calling_peak/ -p ' + sample + ' -s + ' + ' && \n' \
    + 'head -n 12 ' + macs2_enrich_txt + ' |grep -E "3UTR|TTS|Exon|Intron|Intergenic|Promoter|5UTR|Ann" |sed \'1c Annotation\\tPeak_num\\tTotal_size\\tLog2_Enrichment\' >' +macs2_enrich_p1 + ' && \n' \
    + 'tail -n 34 ' + macs2_enrich_txt + ' |grep -E "SINE|LINE|LTR|Satellite|rRNA|RC|Low_complexity|Simple_repeat|Other|DNA|RNA|snoRNA|AnnncRNA|srpRNA|tRNA|snRNA|scRNA|CpG-Island|Unknown" |grep -v ? |sed \'1c Annotation\\tPeak_num\\tTotal_size\\tLog2_Enrichment\' >' +macs2_enrich_p2 + ' && \n' \
    +'echo "Calling peak and evaluate peak state of '+sample+' has been done!" '

    MACS2_shell_file=open(qsu_dir+'/09_MACS2_calling_peak/MACS2.'+sample+'.sh','w')
    MACS2_shell_file.write(MACS2_str)
    MACS2_shell_file.close()
    os.system('chmod 755 ' + qsu_dir+'/09_MACS2_calling_peak/MACS2.'+sample+'.sh')


''' 10.  Deeptools bamCoverage and plotCoverage  '''
for sample in samples.keys():
    merge_bam = out_dir + '/' + sample + '/01_merge_bam/' + sample + '.merge.bam'
    ensure_dir(out_dir+'/'+sample+'/10_Deeptools_bamCoverage')
    ensure_dir(qsu_dir+'/10_Deeptools_bamCoverage')

    ''' 10.1  deeptools_bin_10_extend_250  '''
    ensure_dir(out_dir+'/'+sample+'/10_Deeptools_bamCoverage/01_bin_10_ext250')
    ensure_dir(qsu_dir+'/10_Deeptools_bamCoverage/01_bin_10_ext250')

    deeptools_bin_10_extend_250_bedGraph=out_dir+'/'+sample+'/10_Deeptools_bamCoverage/01_bin_10_ext250/'+sample+'.bin_10_extend_250.bedGraph'
    deeptools_bin_10_extend_250_BigWig=out_dir+'/'+sample+'/10_Deeptools_bamCoverage/01_bin_10_ext250/'+sample+'.bin_10_extend_250.bw'

    bin_10_extend_250_log=out_dir+'/'+sample+'/10_Deeptools_bamCoverage/01_bin_10_ext250/'+sample+'.bin_10.extend_250.log'

    bedGraph_bigwig_bin_10_extend_250_str='#!/bin/sh\necho "Start deeptools bamCoverage of '+sample+' at bin size 10 bp !" && \n'\
    +bamCoverage+' -b '+ merge_bam +' -o '+deeptools_bin_10_extend_250_BigWig+' --normalizeUsingRPKM --binSize 10 --extendReads 250 --ignoreForNormalization chrX chrM chrY 2> '+bin_10_extend_250_log+' && \n'\
    +bamCoverage+' -b '+ merge_bam +' --outFileFormat bedgraph -o '+deeptools_bin_10_extend_250_bedGraph+' --normalizeUsingRPKM --binSize 10 --extendReads 250 --ignoreForNormalization chrX chrM chrY 2>> '+bin_10_extend_250_log+' && \n'\
    +'echo "Finish deeptools bamCoverage of '+sample+' at bin size 10 bp!" '

    bedGraph_bigwig_bin_10_extend_250_shell_file=open(qsu_dir+'/10_Deeptools_bamCoverage/01_bin_10_ext250/bedGraph_bigwig_bin_10_ext250.'+sample+'.sh','w')
    bedGraph_bigwig_bin_10_extend_250_shell_file.write(bedGraph_bigwig_bin_10_extend_250_str)
    bedGraph_bigwig_bin_10_extend_250_shell_file.close()
    os.system('chmod 755 ' + qsu_dir+'/10_Deeptools_bamCoverage/01_bin_10_ext250/bedGraph_bigwig_bin_10_ext250.'+sample+'.sh')

    ''' 7.2  deeptools_bin_100 '''
    ensure_dir(out_dir + '/' + sample + '/10_Deeptools_bamCoverage/02_bin_100')
    ensure_dir(qsu_dir + '/10_Deeptools_bamCoverage/02_bin_100')
    
    deeptools_bin_100_bedGraph = out_dir + '/' + sample + '/10_Deeptools_bamCoverage/02_bin_100/' + sample + '.no_chrY_unique.sort.rmdup.bin_100.bedGraph'
    deeptools_bin_100_bedGraph_chr = out_dir + '/' + sample + '/10_Deeptools_bamCoverage/02_bin_100/' + sample + '.no_chrY_unique.sort.rmdup.bin_100_chr.bedGraph'
    deeptools_bin_100_bedGraph_chr_zcore = out_dir + '/' + sample + '/10_Deeptools_bamCoverage/02_bin_100/' + sample + '.no_chrY_unique.sort.rmdup.bin_100_chr_Zscore.bedGraph'
    
    deeptools_bin_100_BigWig = out_dir + '/' + sample + '/10_Deeptools_bamCoverage/02_bin_100/' + sample + '.no_chrY_unique.sort.rmdup.bin_100.bw'
    deeptools_bin_100_BigWig_Zscore = out_dir + '/' + sample + '/10_Deeptools_bamCoverage/02_bin_100/' + sample + '.no_chrY_unique.sort.rmdup.bin_100_Zscore.bw'
    
    promoter_5K = out_dir + '/' + sample + '/10_Deeptools_bamCoverage/02_bin_100/' + sample + '.promoter_5K.txt'
    promoter_3K = out_dir + '/' + sample + '/10_Deeptools_bamCoverage/02_bin_100/' + sample + '.promoter_3K.txt'
    promoter_5K_Zscore = out_dir + '/' + sample + '/10_Deeptools_bamCoverage/02_bin_100/' + sample + '.promoter_5K_Zscore.txt'
    promoter_3K_Zscore = out_dir + '/' + sample + '/10_Deeptools_bamCoverage/02_bin_100/' + sample + '.promoter_3K_Zscore.txt'
    
    bin_100_log = out_dir + '/' + sample + '/10_Deeptools_bamCoverage/02_bin_100/' + sample + '.bin_100.log'
    
    bedGraph_bigwig_bin_100_str = '#/bin/sh\necho "Start deeptools bamCoverage of ' + sample + ' at bin size 100 bp !" && \n' \
    + bamCoverage + ' -p 24 -b ' + merge_bam + ' --outFileFormat bedgraph -o ' + deeptools_bin_100_bedGraph + ' --normalizeUsingRPKM --binSize 100 --extendReads 250 --ignoreForNormalization chrM chrY 2>> ' + bin_100_log + ' &&\n' \
    + bamCoverage + ' -p 24 -b ' + merge_bam + ' -o ' + deeptools_bin_100_BigWig + ' --normalizeUsingRPKM --binSize 100 --extendReads 250 --ignoreForNormalization chrM chrY 2> ' + bin_100_log + ' && \n' \
    + bedtools + '/bedtools intersect -a /work2/liugh/liuzunpeng/04_Softwares/bedtools2/genomes/human.hg19.genome.w100.bed -b '+deeptools_bin_100_bedGraph+' -wao |awk \'OFS="\\t" {print $1,$2,$3,$10}\'|sort -k1,1 -k2,2n >'+deeptools_bin_100_bedGraph_chr + ' && \n' \
    + "/work2/liugh/liuzunpeng/04_Softwares/python/local/python2/bin/python2 /work2/liugh/liuzunpeng/04_Softwares/source/Zscore.py -b "+deeptools_bin_100_bedGraph_chr +' -o '+ deeptools_bin_100_bedGraph_chr_zcore  + ' && \n' \
    + bedGraphToBigWig + ' ' + deeptools_bin_100_bedGraph_chr_zcore + ' /work2/liugh/liuzunpeng/04_Softwares/bedtools2/genomes/human.hg19.genome '+  deeptools_bin_100_BigWig_Zscore + ' && \n' \
    + bedtools + '/bedtools intersect -a /work2/liugh/liuzunpeng/03_database/00_genome/01_hg19/03_Ensembl/01_fa_gtf_bed_chr/00_bed/gtf_2_bed.20190908/Homo_sapiens.GRCh37.87.final_promoter_5000.bed -b '+ deeptools_bin_100_bedGraph_chr +' -wao |awk \'OFS="\\t" {print $1,$2,$3,$4,$5,$6,$10}\'| '+bedtools + '/bedtools groupby -g 1,2,3,4,5,6  -c 7 -o mean |awk \'OFS="\\t" {print $1,$2,$3,$4,$7,$6}\' >'+promoter_5K + ' && \n' \
    + bedtools + '/bedtools intersect -a /work2/liugh/liuzunpeng/03_database/00_genome/01_hg19/03_Ensembl/01_fa_gtf_bed_chr/00_bed/gtf_2_bed.20190908/Homo_sapiens.GRCh37.87.final_promoter_3000.bed -b '+ deeptools_bin_100_bedGraph_chr +' -wao |awk \'OFS="\\t" {print $1,$2,$3,$4,$5,$6,$10}\'| '+bedtools + '/bedtools groupby -g 1,2,3,4,5,6  -c 7 -o mean |awk \'OFS="\\t" {print $1,$2,$3,$4,$7,$6}\' >'+promoter_3K + ' && \n' \
    + bedtools + '/bedtools intersect -a /work2/liugh/liuzunpeng/03_database/00_genome/01_hg19/03_Ensembl/01_fa_gtf_bed_chr/00_bed/gtf_2_bed.20190908/Homo_sapiens.GRCh37.87.final_promoter_5000.bed -b '+ deeptools_bin_100_bedGraph_chr_zcore +' -wao |awk \'OFS="\\t" {print $1,$2,$3,$4,$5,$6,$10}\'| '+bedtools + '/bedtools groupby -g 1,2,3,4,5,6  -c 7 -o mean |awk \'OFS="\\t" {print $1,$2,$3,$4,$7,$6}\' >'+promoter_5K_Zscore + ' && \n' \
    + bedtools + '/bedtools intersect -a /work2/liugh/liuzunpeng/03_database/00_genome/01_hg19/03_Ensembl/01_fa_gtf_bed_chr/00_bed/gtf_2_bed.20190908/Homo_sapiens.GRCh37.87.final_promoter_3000.bed -b '+ deeptools_bin_100_bedGraph_chr_zcore +' -wao |awk \'OFS="\\t" {print $1,$2,$3,$4,$5,$6,$10}\'| '+bedtools + '/bedtools groupby -g 1,2,3,4,5,6  -c 7 -o mean |awk \'OFS="\\t" {print $1,$2,$3,$4,$7,$6}\' >'+promoter_3K_Zscore + ' && \n' \
    + 'echo "Finish deeptools bamCoverage of ' + sample + ' at bin size 100 bp!" '
    
    bedGraph_bigwig_bin_100_shell_file = open(qsu_dir + '/10_Deeptools_bamCoverage/02_bin_100/bedGraph_bigwig_bin_100.' + sample + '.sh', 'w')
    bedGraph_bigwig_bin_100_shell_file.write(bedGraph_bigwig_bin_100_str)
    bedGraph_bigwig_bin_100_shell_file.close()
    os.system('chmod 755 ' + qsu_dir + '/10_Deeptools_bamCoverage/02_bin_100/bedGraph_bigwig_bin_100.' + sample + '.sh')
    


    ''' 10.3  deeptools_bin2000  '''
    ensure_dir(out_dir+'/'+sample+'/10_Deeptools_bamCoverage/03_bin_2000')
    ensure_dir(qsu_dir+'/10_Deeptools_bamCoverage/03_bin_2000')

    deeptools_bin_2000_bedGraph=out_dir+'/'+sample+'/10_Deeptools_bamCoverage/03_bin_2000/'+sample+'.bin_2000.bedGraph'
    deeptools_bin_2000_BigWig=out_dir+'/'+sample+'/10_Deeptools_bamCoverage/03_bin_2000/'+sample+'.bin_2000.bw'

    bin_2000_log=out_dir+'/'+sample+'/10_Deeptools_bamCoverage/03_bin_2000/'+sample+'.bin_2000.log'

    bedGraph_bigwig_bin_2000_str='#!/bin/sh\necho "Start deeptools bamCoverage of '+sample+' at bin size 2000 bp !" && \n'\
    +bamCoverage+' -b '+ merge_bam + ' -o '+deeptools_bin_2000_BigWig+' --normalizeUsingRPKM --binSize 2000 --ignoreForNormalization chrX chrM chrY 2> '+bin_2000_log+' && \n'\
    +bamCoverage+' -b '+ merge_bam + ' --outFileFormat bedgraph -o '+deeptools_bin_2000_bedGraph+' --normalizeUsingRPKM --binSize 2000 --ignoreForNormalization chrX chrM chrY 2> '+bin_2000_log+' && \n'\
    +'echo "Finish deeptools bamCoverage of '+sample+' at bin size 2000 bp!" '

    bedGraph_bigwig_bin_2000_shell_file=open(qsu_dir+'/10_Deeptools_bamCoverage/03_bin_2000/bedGraph_bigwig_bin_2000.'+sample+'.sh','w')
    bedGraph_bigwig_bin_2000_shell_file.write(bedGraph_bigwig_bin_2000_str)
    bedGraph_bigwig_bin_2000_shell_file.close()
    os.system('chmod 755 ' + qsu_dir+'/10_Deeptools_bamCoverage/03_bin_2000/bedGraph_bigwig_bin_2000.'+sample+'.sh')
