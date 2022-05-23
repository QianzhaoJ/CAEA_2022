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
import sys,pysam
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


def all_peak_to_fa(A,o,r,p):
    fa = pysam.FastaFile(r)
    all_out = open(o + '/' + p + '.all.peak.fa','w')
    with open(A,'r') as bed_f1:
         for lines in bed_f1:
            lines=lines.strip()
            if lines.startswith('chr') and ("start" not in lines):
                line=lines.split('\t')
                if float(line[7])>=5:
                    seq=fa.fetch(str(line[0]), int(line[1])-1, int(line[2]))
                    all_out.write('>'+'_'.join(line)+'\n')
                    all_out.write(seq+'\n')
    all_out.close()

def fw_peak_to_fa(B,o,r,p):
    fa = pysam.FastaFile(r)
    fw_out = open(o + '/' + p + '.forward.peak.fa', 'w')
    with open(B,'r') as bed_f1:
         for lines in bed_f1:
            lines=lines.strip()
            if lines.startswith('chr') and ("start" not in lines):
                line=lines.split('\t')
                if float(line[7])>=5:
                    seq=fa.fetch(str(line[0]), int(line[1])-1, int(line[2]))
                    fw_out.write('>'+'_'.join(line)+'\n')
                    fw_out.write(seq+'\n')
    fw_out.close()

def rev_region_to_fa(C , o, r, p):
    fa = pysam.FastaFile(r)
    rev_out = open(o + '/' + p + '.reverse.peak.fa', 'w')
    with open(C,'r') as bed_f2:
         for lines in bed_f2:
            lines=lines.strip()
            if lines.startswith('chr') and ("start" not in lines):
                line=lines.split('\t')
                if float(line[7])>=5 :
                    seq=fa.fetch(str(line[0]), int(line[1])-1, int(line[2]))
                    my_seq = Seq(seq, IUPAC.unambiguous_dna)
                    my_seq2= str(my_seq.reverse_complement())
                    rev_out.write('>'+'_'.join(line)+'\n')
                    rev_out.write(my_seq2+'\n')


    rev_out.close()


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

    optparser.add_option('-A', '--bed1', dest='bed1', default=None,help="\nGiven the bed file path.")
    optparser.add_option('-B', '--bed2', dest='bed2', default=None,help="\nGiven the bed file path.")
    optparser.add_option('-C', '--bed3', dest='bed3', default=None, help="\nGiven the bed file path.")
    optparser.add_option('-o', '--out', dest='out', default=None,help="\nGiven the out directory of out fa path.")
    optparser.add_option('-r', '--ref', dest='ref', default='/work2/liugh/liuzunpeng/03_database/00_genome/01_hg19/03_Ensembl/01_fa_gtf_bed_chr/hg19.fa',help="\nGiven the reference fa.")
    optparser.add_option('-p', '--prefix', dest='prefix',default='',help="\nGiven the reference fa.")
    optparser.add_option("-h", '--help', action="help",help="\nShow this help message and exit.")
    return optparser


def main(opt):
    """
    Calculate ATGC content and GC skew of input bedgraph
    """
    All_bed    = opt.bed1
    forward_bed=opt.bed2
    reverse_bed= opt.bed1
    fasta      =opt.ref
    out_path   =opt.out
    prefix     =opt.prefix

    all_peak_to_fa(All_bed, out_path, fasta, prefix)
    fw_peak_to_fa(forward_bed, out_path, fasta, prefix)
    rev_region_to_fa(reverse_bed, out_path, fasta, prefix)


if __name__ == '__main__':
    prepare_optparser()
    (options, args) = prepare_optparser().parse_args()

    main(options)
    sys.exit()
