#!/usr/bin/env python3
# coding: utf-8

from optparse import OptionParser
import os
import sys
import re

'''
Title:        DamID.seq
Description:  This pipeline is used for ChIP-seq analysis and downstream analysis.
    
Author:       Zunpeng Liu
Email:        zunpengliuAT163.com
Date:         22/3/2019
Version:      2.0
Python:       based on python 3.6
Citation:
'''


MY_USAGE='''
Author:       Zunpeng Liu
Email:        zunpengliu@163.com
Date:         1/6/2018
    
Blessings:    If you have any question or query of this pipeline, please feel free to contact with me.
'''

parser=OptionParser(MY_USAGE,version='Version 1.1')
parser.add_option('-i','--in_dir',   dest='in_dir',   type='string', help='You should sepecify the input file path, wich stored all raw data.')
parser.add_option('-p','--prefix',   dest='prefix',   type='string', help='You should specifty the prefix of samples, which will be used to label all processed results. For example, if the prefix was defined as "M91_HHVK7CCXY_L1", the prefix of file names will be precessed as "M91".')
parser.add_option('-o','--out_dir',  dest='out_dir',  type='string', help='You should specifty the output directory, which will store all processed results.')
parser.add_option('-q','--qsub_dir', dest='qsub_dir', type='string', help='You should specifty the qub directory, which will store all shell script derived from this pipeline.')
parser.add_option('-c','--configure_file',dest='configure_file',type='string', default='/pfs1/liuguanghui/liuzunpeng/01_Script/03_ATAC-seq/07_ATAC_seq_pipeline/config.txt',help='The configure lib file should include the absolute path of softwares used in this pipeline')
(options,args)=parser.parse_args()


''' Predefine variables '''
in_dir  = str(options.in_dir)
prefix  = str(options.prefix)
sample=re.findall(r'([a-zA-Z0-9\-]+)?_.*',prefix)[0]
out_dir = str(options.out_dir)
qsu_dir = str(options.qsub_dir)
log_dir = '/home/liuzunpeng/../..'+qsu_dir

''' ensure all directory has been created, if not exist,then creat the new ones. '''
def ensure_dir(directory):
    #directory = os.path.dirname(f)
    if not os.path.exists(directory):
        os.makedirs(directory)

''' make new directions '''
ensure_dir(out_dir)
ensure_dir(qsu_dir)

''' Src '''

''' Configure all softwares and library '''

sft={}
config= options.configure_file
with open(config,'r') as config:
    for line in config:
        #        if len(line)!=0 and not line.startswith('#'):
        if len(re.findall('=',line))==1:
            sft[re.findall(r'(\w+?)=.*',line)[0]]=re.findall(r'\w+=(.*)',line)[0]
globals().update(sft)

''' Raw data files '''
fq1=in_dir+'/'+prefix+'_1.fq.gz'
fq2=in_dir+'/'+prefix+'_2.fq.gz'


''' 1. FastQC rawdata '''
ensure_dir(str(out_dir+'/'+sample+'/01_fastQC_raw'))
ensure_dir(str(qsu_dir+'/01_fastQC_raw'))

fqv1=out_dir+'/'+sample+'/01_fastQC_raw/'+sample+'.R1.basic.xls'
fqv2=out_dir+'/'+sample+'/01_fastQC_raw/'+sample+'.R2.basic.xls'

FastQC_str='\
#/bin/sh\n\
\
echo "Start fastQC and fqvalue raw data of '+sample+'" && \n'\
+fastqc+' --outdir '+out_dir+'/'+sample+'/01_fastQC_raw -f fastq '+fq1+' '+fq2+' && \n'\
+fqvalue+' -q 33 '+fq1+' > '+fqv1+' && \n'\
+fqvalue+' -q 33 '+fq2+' > '+fqv2+' && \n'\
+'echo "FastQC and fqvalue raw data of '+sample+' has been done! " '

FastQC_shell_file=open(qsu_dir+'/01_fastQC_raw/fastQC_raw.'+sample+'.sh','w')
FastQC_shell_file.write(FastQC_str)
FastQC_shell_file.close()
os.system('chmod 755 ' + qsu_dir+'/01_fastQC_raw/fastQC_raw.'+sample+'.sh')


''' 2. TrimGalore clean '''
ensure_dir(out_dir+'/'+sample+'/02_TrimGalore')
ensure_dir(qsu_dir+'/02_TrimGalore')

clean_fq1=out_dir+'/'+sample+'/02_TrimGalore/'+prefix+'_1_val_1.fq.gz'
clean_fq2=out_dir+'/'+sample+'/02_TrimGalore/'+prefix+'_2_val_2.fq.gz'

trim_galore_log=out_dir+'/'+sample+'/02_TrimGalore/'+sample+'.trim_galore.log'

trim_galore_str='\
#/bin/sh\n\
\
echo "Start trim raw data of '+sample+'" && \n'\
+trim_galore+' --fastqc --path_to_cutadapt '+cutadapt+' --stringency 3 --paired --output_dir '+out_dir+'/'+sample+'/02_TrimGalore '\
+fq1+' '+fq2+' 2>'+trim_galore_log+' && \n\
echo "Clean and FastQC of '+sample+' done! All done! \"'


trim_galore_shell_file=open(qsu_dir+'/02_TrimGalore/trim_galore.'+sample+'.sh','w')
trim_galore_shell_file.write(trim_galore_str)
trim_galore_shell_file.close()
os.system('chmod 755 ' + qsu_dir+'/02_TrimGalore/trim_galore.'+sample+'.sh')


''' 3. FastQC clean '''
ensure_dir(out_dir+'/'+sample+'/03_fastQCclean')
ensure_dir(qsu_dir+'/03_fastQCclean')

fqvC1=out_dir+'/'+sample+'/03_fastQCclean/'+sample+'.R1.clean.basic.xls'
fqvC2=out_dir+'/'+sample+'/03_fastQCclean/'+sample+'.R2.clean.basic.xls'

Basicresult=out_dir+'/'+sample+'/03_fastQCclean/'+sample+'.BasicInfo.xls'

FastQC_clean_str='\
#/bin/sh\n\
\
echo "Start FastQC Clean data of '+sample+'" && \n'\
+fqvalue+' -q 33 '+clean_fq1+' > '+fqvC1+' && \n'\
+fqvalue+' -q 33 '+clean_fq2+' > '+fqvC2+' && \n\
echo "Fqvalue raw data of '+sample+' has been done! \"'

FastQC_clean_shell_file=open(qsu_dir+'/03_fastQCclean/fastQCclean.'+sample+'.sh','w')
FastQC_clean_shell_file.write(FastQC_clean_str)
FastQC_clean_shell_file.close()
os.system('chmod 755 ' + qsu_dir+'/03_fastQCclean/fastQCclean.'+sample+'.sh')


''' 4. Alignment to reference by bowtie2 '''
ensure_dir(out_dir+'/'+sample+'/04_bowtie2')
ensure_dir(qsu_dir+'/04_bowtie2')

bowtie2_Align_log=out_dir+'/'+sample+'/04_bowtie2/'+sample+'.bowtie2_Align_log.txt'

fragment_size    =out_dir+'/'+sample+'/04_bowtie2/'+sample+'.fragment_size.txt'
fragment_size_tsv=out_dir+'/'+sample+'/04_bowtie2/'+sample+'.fragment_size.tsv'

sam=out_dir+'/'+sample+'/04_bowtie2/'+sample+'.sam'

bowtie2_str='\
#/bin/sh\n\
\
echo "Start the alignment of '+sample+'to reference by bowtie2!" && \n'\
+bowtie2+' -p 24 -x '+bowtie2_index+' --no-mixed --no-discordant -t -1 '+clean_fq1+' -2 '+clean_fq2+' -S ' +sam+ ' 2>'+bowtie2_Align_log+' && \n\
echo "Alignment to reference of '+sample+' by bowtie2 has been done! \"'

bowtie2_shell_file=open(qsu_dir+'/04_bowtie2/bowtie2.'+sample+'.sh','w')
bowtie2_shell_file.write(bowtie2_str)
bowtie2_shell_file.close()
os.system('chmod 755 ' + qsu_dir+'/04_bowtie2/bowtie2.'+sample+'.sh')


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

unique_reads_str='\
#/bin/sh\n\
\
echo "Start transfering the sam file of '+sample+' to bam file and remove chrY DNA!" && \n'\
+'awk \'$3!="chrY" \' '+sam+' | '+samtools+' view -@ 24 -S -b -q 10 > '+bam+' && \n'\
+'echo "Start filtering unique mapped reads of '+sample+'!" && \n'\
+'awk \'$3!="chrY" \' '+sam+' | grep -v "XS:i:" | '+samtools+' view -@ 24 -S -b -q 10 > '+unique_bam+' && \n'\
+'echo "Start sorting no_chrY_unique_rmdup.bam of '+sample+' on position by samtools!" && \n'\
+samtools+' sort -@ 24 -l 9 '+bam+' -o '+srt_bam+' && \n'\
+'echo "Start removing duplicates of '+sample+' by picard MarkDuplicates tools!" && \n'\
+java+' -jar '+picard+' REMOVE_DUPLICATES=true INPUT='+srt_bam+' METRICS_FILE='+METRICS_FILE+' OUTPUT='+srt_rmdup_bam+' && \n'\
+samtools+' index -@ 24 '+srt_rmdup_bam+' && \n'\
+samtools+' flagstat -@ 24  '+srt_rmdup_bam+' > '+flagstat_txt+' && \n'\
+'echo "Samtools and Picard of '+sample+' done!" '

unique_reads_shell_file=open(qsu_dir+'/05_unique_mapped_reads/unique_mapped.'+sample+'.sh','w')
unique_reads_shell_file.write(unique_reads_str)
unique_reads_shell_file.close()
os.system('chmod 755 ' + qsu_dir+'/05_unique_mapped_reads/unique_mapped.'+sample+'.sh')



''' 6. Read length distribution and plot density '''
ensure_dir(out_dir+'/'+sample+'/06_fragment_distribution')
ensure_dir(qsu_dir+'/06_fragment_distribution')

fragment_size=out_dir+'/'+sample+'/06_fragment_distribution/'+sample+'.fragment.size.txt'
fragment_size_tsv=out_dir+'/'+sample+'/06_fragment_distribution/'+sample+'.fragment.size.tsv'

fragment_str='\
#/bin/sh\n\
\
echo "Start fragment_distribution analysis of '+sample+' !" && \n'\
+'awk \'$3!="chrY" && $9!="0" {print $9} \' '+out_dir+'/'+sample+'/04_bowtie2/'+sample+'.sam |grep -v "-" > '+fragment_size+' && \n\
sed \'1i fragment_size\' '+fragment_size+' > '+fragment_size_tsv+' && \n'\
+Rscript+' '+Fragment_analysis+' '+fragment_size_tsv+' '+out_dir+'/'+sample+'/06_fragment_distribution/'+sample+' && \n'\
+'echo "Analysis read length distribution of '+sample+' done!" '

fragment_shell_file=open(qsu_dir+'/06_fragment_distribution/fragment_analysis.'+sample+'.sh','w')
fragment_shell_file.write(fragment_str)
fragment_shell_file.close()
os.system('chmod 755 ' + qsu_dir+'/06_fragment_distribution/fragment_analysis.'+sample+'.sh')




''' 7. Bam to bed and remove blacklisted-intervals '''
ensure_dir(out_dir+'/'+sample+'/07_Bam2bed_remove_black_list')
ensure_dir(qsu_dir+'/07_Bam2bed_remove_black_list')

srt_rmdup_bed=out_dir+'/'+sample+'/07_Bam2bed_remove_black_list/'+sample+'_no_chrM_chrY_unique.sort.rmdup.bed'
interdect_BL_bed=out_dir+'/'+sample+'/07_Bam2bed_remove_black_list/'+sample+'_no_chrM_chrY_unique.sort.rmdup.interdect_BL.bed'
remove_BL_bed=out_dir+'/'+sample+'/07_Bam2bed_remove_black_list/'+sample+'_no_chrM_chrY_unique.sort.rmdup.remove_BL.bed'
bam2bed_log=out_dir+'/'+sample+'/07_Bam2bed_remove_black_list/'+sample+'.no_chrM_chrY_unique.sort.rmdup.bam2bed.log.txt'

Bam2bed_remove_black_list_str='\
#/bin/sh\n\
echo "Start converting the bam file of '+sample+' to bed file!" && \n'\
+bamToBed + ' -i ' + srt_rmdup_bam + ' > ' + srt_rmdup_bed +' && \n'\
+'echo "Start remove black list of '+sample+'! " && \n'\
+bedtools+'/bedtools intersect -u -a ' + srt_rmdup_bed + ' -b ' + black_list + ' > ' + interdect_BL_bed +' && \n'\
+bedtools+'/bedtools intersect -v -a ' + srt_rmdup_bed + ' -b ' + black_list + ' > ' + remove_BL_bed +' && \n'\
+'echo "Convert the bam file of '+sample+' to bed file and remove black list has been done!" '

Bam2bed_shell_file=open(qsu_dir+'/07_Bam2bed_remove_black_list/Bam2bed_remove_black_list.'+sample+'.sh','w')
Bam2bed_shell_file.write(Bam2bed_remove_black_list_str)
Bam2bed_shell_file.close()
os.system('chmod 755 ' + qsu_dir+'/07_Bam2bed_remove_black_list/Bam2bed_remove_black_list.'+sample+'.sh')



''' 7. Bedtools_genomeCoverageBed bedGraph, bedGraphToBigWig bigwig '''
ensure_dir(out_dir+'/'+sample+'/07_Bedtools_bedGraph_bigwig')
ensure_dir(qsu_dir+'/07_Bedtools_bedGraph_bigwig')

Bedtools_bedGraph_file=out_dir+'/'+sample+'/07_Bedtools_bedGraph_bigwig/'+sample+'.no_chrY.sort.rmdup.bedGraph'
bedGraphToBigWig_file =out_dir+'/'+sample+'/07_Bedtools_bedGraph_bigwig/'+sample+'.no_chrY.sort.rmdup.bw'

bedGraph_bigwig_str='\
#/bin/sh\n\
\
echo "Start converting the bam file of '+sample+' to bedGraph file by Bedtools_genomeCoverageBed !" && \n'\
+genomeCoverageBed+' -bg -split -ibam '+srt_rmdup_bam+' -g '+chromosomesize+' > '+Bedtools_bedGraph_file+' && \n'\
+'echo "Make bedGraph of '+sample+'has been done!" && \n'\
+'echo "bedGraphToBigWig convert bedGraph of '+sample+'!" && \n'\
+bedGraphToBigWig+' '+Bedtools_bedGraph_file+' '+ chromosomesize + ' '+bedGraphToBigWig_file+' && \n'\
+'echo "bedGraphToBigWig make bedGraph of '+sample+' has been done!" '

bedGraph_bigwig_shell_file=open(qsu_dir+'/07_Bedtools_bedGraph_bigwig/bedGraph_bigwig.'+sample+'.sh','w')
bedGraph_bigwig_shell_file.write(bedGraph_bigwig_str)
bedGraph_bigwig_shell_file.close()
os.system('chmod 755 ' + qsu_dir+'/06_fragment_distribution/fragment_analysis.'+sample+'.sh')


######

''' 9.  Deeptools bamCoverage and plotCoverage  '''
ensure_dir(out_dir+'/'+sample+'/09_Deeptools_bamCoverage')
ensure_dir(qsu_dir+'/09_Deeptools_bamCoverage')

''' 9.1  deeptools_bin_10  '''
ensure_dir(out_dir+'/'+sample+'/09_Deeptools_bamCoverage/01_bin_10')
ensure_dir(qsu_dir+'/09_Deeptools_bamCoverage/01_bin_10')

deeptools_bin_10_bedGraph=out_dir+'/'+sample+'/09_Deeptools_bamCoverage/01_bin_10/'+sample+'.no_chrY_unique.sort.rmdup.bin_10.bedGraph'
deeptools_bin_10_BigWig=out_dir+'/'+sample+'/09_Deeptools_bamCoverage/01_bin_10/'+sample+'.no_chrY_unique.sort.rmdup.bin_10.bw'

bin_10_log=out_dir+'/'+sample+'/09_Deeptools_bamCoverage/01_bin_10/'+sample+'.bin_10.log'

bedGraph_bigwig_bin_10_str='\
#/bin/sh\n\
\
echo "Start deeptools bamCoverage of '+sample+' at bin size 10 bp !" && \n'\
+bamCoverage+' -p 24  -b '+srt_rmdup_bam+' -o '+deeptools_bin_10_BigWig+' --normalizeUsingRPKM --binSize 10 --ignoreForNormalization chrX chrM chrY 2> '+bin_10_log+' &&\n'\
+bamCoverage+' -p 24  -b '+srt_rmdup_bam+' --outFileFormat bedgraph -o '+deeptools_bin_10_bedGraph+' --normalizeUsingRPKM --binSize 10 --ignoreForNormalization chrX chrM chrY 2>> '+bin_10_log+' &&\n'\
+'echo "Finish deeptools bamCoverage of '+sample+' at bin size 10 bp!" '

bedGraph_bigwig_bin_10_shell_file=open(qsu_dir+'/09_Deeptools_bamCoverage/01_bin_10/bedGraph_bigwig_bin_10.'+sample+'.sh','w')
bedGraph_bigwig_bin_10_shell_file.write(bedGraph_bigwig_bin_10_str)
bedGraph_bigwig_bin_10_shell_file.close()
os.system('chmod 755 ' + qsu_dir+'/09_Deeptools_bamCoverage/01_bin_10/bedGraph_bigwig_bin_10.'+sample+'.sh')



''' 9.2  deeptools_bin_100 '''
ensure_dir(out_dir+'/'+sample+'/09_Deeptools_bamCoverage/02_bin_100')
ensure_dir(qsu_dir+'/09_Deeptools_bamCoverage/02_bin_100')

deeptools_bin_100_bedGraph=out_dir+'/'+sample+'/09_Deeptools_bamCoverage/02_bin_100/'+sample+'.no_chrY_unique.sort.rmdup.bin_100.bedGraph'
deeptools_bin_100_BigWig=out_dir+'/'+sample+'/09_Deeptools_bamCoverage/02_bin_100/'+sample+'.no_chrY_unique.sort.rmdup.bin_100.bw'

bin_100_log=out_dir+'/'+sample+'/09_Deeptools_bamCoverage/02_bin_100/'+sample+'.bin_100.log'

bedGraph_bigwig_bin_100_str='\
#/bin/sh\n\
\
echo "Start deeptools bamCoverage of '+sample+' at bin size 100 bp !" && \n'\
+bamCoverage+' -p 24 -b '+srt_rmdup_bam+' -o '+deeptools_bin_100_BigWig+' --normalizeUsingRPKM --binSize 100 --ignoreForNormalization chrX chrM chrY 2> '+bin_100_log+' && \n'\
+bamCoverage+' -p 24 -b '+srt_rmdup_bam+' --outFileFormat bedgraph -o '+deeptools_bin_100_bedGraph+' --normalizeUsingRPKM --binSize 10 --ignoreForNormalization chrX chrM chrY 2>> '+bin_100_log+' &&\n'\
+'echo "Finish deeptools bamCoverage of '+sample+' at bin size 100 bp!" '

bedGraph_bigwig_bin_100_shell_file=open(qsu_dir+'/09_Deeptools_bamCoverage/02_bin_100/bedGraph_bigwig_bin_100.'+sample+'.sh','w')
bedGraph_bigwig_bin_100_shell_file.write(bedGraph_bigwig_bin_100_str)
bedGraph_bigwig_bin_100_shell_file.close()
os.system('chmod 755 ' + qsu_dir+'/09_Deeptools_bamCoverage/02_bin_100/bedGraph_bigwig_bin_100.'+sample+'.sh')


''' 9.3  deeptools_bin_2000  '''
ensure_dir(out_dir+'/'+sample+'/09_Deeptools_bamCoverage/03_bin_2000')
ensure_dir(qsu_dir+'/09_Deeptools_bamCoverage/03_bin_2000')

deeptools_bin_2000_bedGraph=out_dir+'/'+sample+'/09_Deeptools_bamCoverage/03_bin_2000/'+sample+'.no_chrY_unique.sort.rmdup.bin_2000.bedGraph'
deeptools_bin_2000_BigWig=out_dir+'/'+sample+'/09_Deeptools_bamCoverage/03_bin_2000/'+sample+'.no_chrY_unique.sort.rmdup.bin_2000.bw'

bin_2000_log=out_dir+'/'+sample+'/09_Deeptools_bamCoverage/03_bin_2000/'+sample+'.bin_2000.log'

bedGraph_bigwig_bin_2000_str='\
#/bin/sh\n\
\
echo "Start deeptools bamCoverage of '+sample+' at bin size 2000 bp !" && \n'\
+bamCoverage+' -p 24  -b '+srt_rmdup_bam+' -o '+deeptools_bin_2000_BigWig+' --normalizeUsingRPKM --binSize 2000 --ignoreForNormalization chrX chrM chrY 2> '+bin_2000_log+' &&\n'\
+bamCoverage+' -p 24  -b '+srt_rmdup_bam+' --outFileFormat bedgraph -o '+deeptools_bin_2000_bedGraph+' --normalizeUsingRPKM --binSize 2000 --ignoreForNormalization chrX chrM chrY 2>> '+bin_2000_log+' &&\n'\
+'echo "Finish deeptools bamCoverage of '+sample+' at bin size 2000 bp!" '

bedGraph_bigwig_bin_2000_shell_file=open(qsu_dir+'/09_Deeptools_bamCoverage/03_bin_2000/bedGraph_bigwig_bin_2000.'+sample+'.sh','w')
bedGraph_bigwig_bin_2000_shell_file.write(bedGraph_bigwig_bin_2000_str)
bedGraph_bigwig_bin_2000_shell_file.close()
os.system('chmod 755 ' + qsu_dir+'/09_Deeptools_bamCoverage/03_bin_2000/bedGraph_bigwig_bin_2000.'+sample+'.sh')


''' 10. computeMatrix plotHeatmap plotProfile binSize 100 '''

ensure_dir(out_dir+'/'+sample+'/10_computeMatrix_plotHeatmap')
ensure_dir(out_dir+'/'+sample+'/10_computeMatrix_plotHeatmap/01_TSS')
ensure_dir(out_dir+'/'+sample+'/10_computeMatrix_plotHeatmap/02_Genebody')
ensure_dir(qsu_dir+'/10_computeMatrix_plotHeatmap')


TSS_bin_100_matrix=out_dir+'/'+sample+'/10_computeMatrix_plotHeatmap/01_TSS/'+sample+'.TSS.RPKM.matrix.mat.gz'
TSS_FileSortedRegions=out_dir+'/'+sample+'/10_computeMatrix_plotHeatmap/01_TSS/'+sample+'.TSS.RPKM.Sorted.Regions.bed'

Genebody_bin_100_matrix=out_dir+'/'+sample+'/10_computeMatrix_plotHeatmap/02_Genebody/'+sample+'.Genebody.RPKM.matrix.mat.gz'
Genebody_FileSortedRegions=out_dir+'/'+sample+'/10_computeMatrix_plotHeatmap/02_Genebody/'+sample+'.Genebody.RPKM.Sorted.Regions.bed'

computeMatrix_log=out_dir+'/'+sample+'/10_computeMatrix_plotHeatmap/'+sample+'.computeMatrix.log'

computeMatrix_bin_100_str='\
#/bin/sh\n\
\
echo "Start deeptools computeMatrix of '+sample+' at bin size 100 bp !" && \n'\
+computeMatrix+' reference-point --referencePoint TSS -b 3000 -a 3000 -R '+RefSeq_bed+' -S '+deeptools_bin_100_BigWig+' --skipZeros -o '+TSS_bin_100_matrix+' --outFileSortedRegions '+TSS_FileSortedRegions+' 2> '+computeMatrix_log+' && \n'\
+computeMatrix+' scale-regions --beforeRegionStartLength 21000 --regionBodyLength 35000 --afterRegionStartLength 21000 -R '+RefSeq_bed+' -S '+deeptools_bin_100_BigWig+' --skipZeros -o '+Genebody_bin_100_matrix+' --outFileSortedRegions '+Genebody_FileSortedRegions+' 2>> '+computeMatrix_log+' && \n'\
+'echo "Finish computeMatrix of '+sample+' at bin size 100 bp!"  && \n'\
+plotHeatmap+' -m '+TSS_bin_100_matrix+' --colorMap RdPu -out '+out_dir+'/'+sample+'/10_computeMatrix_plotHeatmap/01_TSS/'+sample+'.TSS.heatmap.pdf --heatmapHeight 26 --heatmapWidth 4 --plotTitle '+sample+' && \n'\
+plotProfile+' -m '+TSS_bin_100_matrix+' --colors blue -out '+out_dir+'/'+sample+'/10_computeMatrix_plotHeatmap/01_TSS/'+sample+'.TSS.profile.pdf --plotTitle '+sample+' && \n'\
+plotHeatmap+' -m '+Genebody_bin_100_matrix+' --colorMap RdPu -out '+out_dir+'/'+sample+'/10_computeMatrix_plotHeatmap/02_Genebody/'+sample+'.TSS.heatmap.pdf --heatmapHeight 26 --heatmapWidth 4 --plotTitle '+sample+' && \n'\
+plotProfile+' -m '+Genebody_bin_100_matrix+' --colors blue -out '+out_dir+'/'+sample+'/10_computeMatrix_plotHeatmap/02_Genebody/'+sample+'.TSS.profile.pdf --plotTitle '+sample+' && \n'\
+'echo "Finish deeptools plotHeatmap and plotProfile TSS of '+sample+' at bin size 100 bp!" '


computeMatrix_plotHeatmap_shell_file=open(qsu_dir+'/10_computeMatrix_plotHeatmap/computeMatrix.'+sample+'.sh','w')
computeMatrix_plotHeatmap_shell_file.write(computeMatrix_bin_100_str)
computeMatrix_plotHeatmap_shell_file.close()
os.system('chmod 755 ' + qsu_dir+'/10_computeMatrix_plotHeatmap/computeMatrix.'+sample+'.sh')


''' 11_damidseq_pipeline '''

ensure_dir(out_dir+'/'+sample+'/11_damidseq_pipeline')
ensure_dir(qsu_dir+'/11_damidseq_pipeline')

control_name='Dam-'+re.findall(r'-(.*)',sample)[0]
control_bam=out_dir+'/'+control_name+'/05_unique_mapped_reads/'+control_name+'_no_chrY.sort.rmdup.bam'
damidseq_pipeline='/work2/liugh/liuzunpeng/01_Script/01_pipeline/07_DamID_seq/DamID_pipeline-master/damidseq_pipeline'

damidseq_pipeline_str='\
#!/bin/sh\n\
echo "Start damidseq_pipeline of '+sample+'!" && \n'\
+damidseq_pipeline+' --load_defaults=hg19 --threads=24 --bamfiles ' +srt_rmdup_bam+ '--dam='+control_bam+' \n'\
+'echo "Finish damidseq_pipeline of '+sample+' !" '


damidseq_pipeline_shell_file=open(qsu_dir+'/11_damidseq_pipeline/damidseq.'+sample+'.sh','w')
damidseq_pipeline_shell_file.write(damidseq_pipeline_str)
damidseq_pipeline_shell_file.close()
os.system('chmod 755 ' + qsu_dir+'/11_damidseq_pipeline/damidseq.'+sample+'.sh')

