# !/usr/bin/env python3
# coding: utf-8

from __future__ import division

__author__ = "Zunpeng Liu"
__copyright__ = "Copyright 2018, IBP, CAS."
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Zunpeng Liu"
__email__ = "zunpengliuAT163.com"



from optparse import OptionParser
import os,sys,pandas


def prepare_optparser():
    usage = """usage: python3  -b Rloop_fr.peak.bed -B  Rloop_rev.peak.bed   --o ./ -r hg19.fa 
        """

    description = "This script can be used to get sequence corresponding to the bed region."

    optparser = OptionParser(
            version="%s v1.0 20180710" % (sys.argv[0]),
            description=description,
            usage=usage,
            add_help_option=False
    )

    optparser.add_option('-b', '--bed1', dest='bed1', default=None,help="\nGiven the bed file path.")
    optparser.add_option('-B', '--bed2', dest='bed2', default=None,help="\nGiven the bed file path.")
    optparser.add_option('-o', '--out', dest='out', default=None,help="\nGiven the out directory of out fa path.")
    optparser.add_option('-r', '--ref', dest='ref', default='/work2/liugh/liuzunpeng/03_database/00_genome/01_hg19/03_Ensembl/01_fa_gtf_bed_chr/00_bed/Homo_sapiens.GRCh37.87_Rloop_all_genes',help="\nGiven the reference fa.")
    optparser.add_option('-p', '--prefix', dest='prefix',default='',help="\nGiven the reference fa.")
    optparser.add_option("-h", '--help', action="help",help="\nShow this help message and exit.")
    return optparser


def main(opt):
    """
    Calculate ATGC content and GC skew of input bedgraph
    """
    forward_bed=opt.bed1
    reverse_bed=opt.bed2
    ref_bed    =opt.ref
    out_path   =opt.out
    prefix     =opt.prefix

    all_ref   = ref_bed + '.bed'
    forwd_ref = ref_bed + '_fw.bed'
    rev_ref   = ref_bed + '_rev.bed'


    #fw_fw_sense1         = out_path + '/' + prefix + '.fw_fw_sense1.txt'
    #fw_rev_antisense1    = out_path + '/' + prefix + '.fw_rev_antisense1.txt'
    #rev_rev_sense2       = out_path + '/' + prefix + '.rev_rev_sense2.txt'
    #rev_fw_antisense2    = out_path + '/' + prefix + '.rev_fw_antisense2.txt'
    sense                 = out_path + '/' + prefix + '.sense.txt'
    antisense             = out_path + '/' + prefix + '.antisense.txt'
    #sense_antisense = out_path + '/' + prefix + '.sense_antisense.overlap_txt'

    ''' Get sense_Rloop & antisense_Rloop '''
    cmd1='bedtools intersect -a ' + forward_bed + ' -b ' + forwd_ref + ' -wo  > ' + sense + ' && ' + \
    'bedtools intersect -a ' + forward_bed + ' -b ' + rev_ref + ' -wo  > ' + antisense + ' && ' + \
    'bedtools intersect -a ' + reverse_bed + ' -b ' + rev_ref + ' -wo  >> ' + sense + ' && ' + \
    'bedtools intersect -a ' + reverse_bed + ' -b ' + forwd_ref + ' -wo  >> ' + antisense
    os.system(cmd1)

    sense_unique = out_path + '/' + prefix + '.sense.unique.txt'
    antisense_unique = out_path + '/' + prefix + '.antisense_unique.txt'

    cmd2='bedtools intersect -a ' + forward_bed + ' -b ' + forwd_ref + ' -u  > ' + sense_unique + ' && ' + \
    'bedtools intersect -a ' + forward_bed + ' -b ' + rev_ref + ' -u  > ' + antisense_unique + ' && ' + \
    'bedtools intersect -a ' + reverse_bed + ' -b ' + rev_ref + ' -u  >> ' + sense_unique + ' && ' + \
    'bedtools intersect -a ' + reverse_bed + ' -b ' + forwd_ref + ' -u  >> ' + antisense_unique
    os.system(cmd2)

    # ''' Get sense and antisense overlaped R-loop peak '''
    # cmd2='bedtools intersect -a ' + sense +  ' -b ' +  antisense + ' > ' + sense_antisense
    # os.system(cmd2)

    ''' Genes without R - loops '''
    Gene_with_sense     = out_path + '/' + prefix + '.genelist_with_sense_peak.txt'
    Gene_with_antisense = out_path + '/' + prefix + '.genelist_with_antisense_peak.txt'
    Gene_with_sense_anti= out_path + '/' + prefix + '.genelist_with_sense_and_anti_peak.txt'
    Gene_without_Rloop  = out_path + '/' + prefix + '.genelist_without_sense_and_anti_peak.txt'

    cmd3 = 'bedtools intersect -a ' + forwd_ref + ' -b ' + forward_bed + ' -u  > ' + Gene_with_sense + ' && ' + \
           'bedtools intersect -a ' + rev_ref + ' -b ' + reverse_bed + ' -u  >> ' + Gene_with_sense + ' && ' + \
           'bedtools intersect -a ' + forwd_ref + ' -b ' + reverse_bed + ' -u  > ' + Gene_with_antisense + ' && ' + \
           'bedtools intersect -a ' + rev_ref + ' -b ' + forward_bed + ' -u  >> ' + Gene_with_antisense + ' && ' + \
           'bedtools intersect -a ' + Gene_with_sense + ' -b ' + Gene_with_antisense + ' -u  >' + Gene_with_sense_anti + ' && ' + \
           'bedtools intersect -v -a ' + all_ref + ' -b ' + Gene_with_sense + ' |bedtools intersect -v  -b ' + Gene_with_antisense  +' > ' + Gene_without_Rloop


    os.system(cmd3)




if __name__ == '__main__':
    prepare_optparser()
    (options, args) = prepare_optparser().parse_args()
    main(options)
    sys.exit()
