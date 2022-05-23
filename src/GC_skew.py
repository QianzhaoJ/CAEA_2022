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
import pandas as pd
import sys,pysam
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


def region_to_fa(b,B,o,r,p):
    fa = pysam.FastaFile(r)
    out2={}
    fr_out = open(o+'/'+p+'.forward.fa','w')
    rev_out= open(o+'/'+p+'.reverse.fa','w')
    with open(b,'r') as bed_f1:
         for lines in bed_f1:
            lines=lines.strip()
            if lines.startswith('chr') and ("start" not in lines):
                line=lines.split('\t')
                if float(line[7])>=5 and (int(line[4]) - 501) >= 0:
                    seq=fa.fetch(str(line[0]), int(line[4])-501, int(line[4])+500)
                    fr_out.write('>'+'_'.join(line)+'\n')
                    fr_out.write(seq+'\n')
                    if 0 in out2.keys():
                        for i in range(0, len(seq)):
                            out2[i].append(seq[i])
                    else:
                        out2 = {i: [seq[i]] for i in range(0, len(seq))}

    with open(B,'r') as bed_f2:
         for lines in bed_f2:
            lines=lines.strip()
            if lines.startswith('chr') and ("start" not in lines):
                line=lines.split('\t')
                if float(line[7])>=5 and int(line[4])-501 >=0 :
                    seq=fa.fetch(str(line[0]), int(line[4])-501, int(line[4])+500)
                    my_seq = Seq(seq, IUPAC.unambiguous_dna)
                    my_seq2= str(my_seq.reverse_complement())
                    rev_out.write('>'+'_'.join(line)+'\n')
                    rev_out.write(my_seq2+'\n')
                    if 0 in out2.keys():
                        for i in range(0, len(my_seq2)):
                            out2[i].append(my_seq2[i])
                    else:
                        out2 = {i: [my_seq2[i]] for i in range(0, len(my_seq2))}

    fr_out.close()
    rev_out.close()
    return out2



def GC_content(out2,out_path,prefix):
    """
    Return GC,GC skew, of input sequence
    """
    dic={}
    dic.fromkeys(out2.keys())
    for k in out2.keys():
        l = len(out2[k])
        A = round(sum(out2[k].count(x) for x in ['A', 'a'])/l, 4)
        G = round(sum(out2[k].count(x) for x in ['G', 'g'])/l, 4)
        C = round(sum(out2[k].count(x) for x in ['C', 'c'])/l, 4)
        T = round(sum(out2[k].count(x) for x in ['T', 't'])/l, 4)
        GC =round(sum(out2[k].count(x) for x in ['G', 'g','C','c','S','s'])/l, 4)
        AT =round(sum(out2[k].count(x) for x in ['A', 'a','T','t','W','w'])/l, 4)
        AT_skew = round((A - T) / (A + T), 4)
        GC_skew = round((G - C) / (G + C), 4)

        dic[int(k)-500]= {'A':A,'G' :G,'C' :C,'T' :T,'GC': GC,'AT': AT,'GC_skew': GC_skew,'AT_skew':AT_skew}

    pd.DataFrame(dic).T.to_csv(out_path+'/'+prefix+'.AGCT_information.csv', sep='\t', index=True,
                                       index_label="Site",columns=['A', 'G', 'C', 'T', 'GC','AT','GC_skew','AT_skew'])



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
    optparser.add_option('-r', '--ref', dest='ref', default='/work2/liugh/liuzunpeng/03_database/00_genome/01_hg19/03_Ensembl/01_fa_gtf_bed_chr/hg19.fa',help="\nGiven the reference fa.")
    optparser.add_option('-p', '--prefix', dest='prefix',default='',help="\nGiven the reference fa.")
    optparser.add_option("-h", '--help', action="help",help="\nShow this help message and exit.")
    return optparser


def main(opt):
    """
    Calculate ATGC content and GC skew of input bedgraph
    """
    forward_bed=opt.bed1
    reverse_bed=opt.bed2
    fasta      =opt.ref
    out_path   =opt.out
    prefix     =opt.prefix

    res1=region_to_fa(forward_bed,reverse_bed,out_path,fasta,prefix)

    GC_content(res1, out_path, prefix)


if __name__ == '__main__':
    prepare_optparser()
    (options, args) = prepare_optparser().parse_args()

    main(options)
    sys.exit()
