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

    optparser.add_option('-i', '--input', dest='input', default=None,help="\nGiven the bed file path.")
    optparser.add_option('-b', '--bed',    dest='bed',    default='/work2/liugh/liuzunpeng/03_database/00_genome/01_hg19/03_Ensembl/01_fa_gtf_bed_chr/00_bed/Homo_sapiens.GRCh37.87_Rloop_all_genes.bed',help="\nGiven the bed file path.")
    optparser.add_option('-o', '--out', dest='out',
                         default='',
                         help="\nGiven the bed file path.")
    optparser.add_option('-p', '--prefix', dest='prefix', default=' ',help="\nGiven the out directory of out fa path.")
    optparser.add_option('-t', '--thread', dest='thread', default=' ',help="\nGiven the thread of multiBamSummary.")
    optparser.add_option("-h", '--help', action="help",help="\nShow this help message and exit.")
    return optparser


def main(opt):
    input    =opt.input
    print(input)
    bed      =opt.bed
    out_path =opt.out
    prefix   =opt.prefix
    thread   =opt.thread

    sense    =pandas.DataFrame()
    antisense=pandas.DataFrame()

    forwd_bed=pandas.DataFrame()
    rev_bed  =pandas.DataFrame()

    out_file_path=out_path + '/' + prefix

    ref_df=pandas.DataFrame(pandas.read_csv(bed,header=0,names=["chr","start","end","gene","score","strand"],sep="\t"))

    for strand in ["fw","rev"]:
        bamlist=[]
        bamlist.append(input + '/' + prefix + '.merge.select.sort.' + strand + '.bam')
        print(bamlist)
        count_xls   = out_path + '/' + prefix + '.merge.select.sort.' + strand + '.counts.xls'
        counts_deal_xls = out_path + '/' + prefix + '.merge.select.sort.' + strand + '.counts_deal.xls'
        addgenename_counts_xls = out_path + '/' + prefix + '.merge.select.sort.' + strand + '.addgenename_counts.xls'
        npz         = out_path + '/' + prefix + '.merge.select.sort.' + strand + '.npz'

        cmd="multiBamSummary BED-file --BED %s -p %s --bamfiles %s --label %s --outRawCounts %s -o %s"%(bed,thread," ".join(bamlist),prefix,count_xls,npz)
        os.system(cmd)
        print(cmd)
        os.system("sed \"1s/[#|']//g\" %s >%s"%(count_xls,counts_deal_xls))

        df = pandas.read_csv(counts_deal_xls, header=0, names=["chr", "start", "end", prefix], sep="\t")

        df3=pandas.merge(df,ref_df,on=["chr","start","end"],how="outer").fillna(0)

        df3.to_csv(addgenename_counts_xls, sep="\t",index=False)
        print("err")

        if strand=="fw":
            print("fw")
            sense=sense.append(df3[df3["strand"]=="+"].ix[:,["gene",prefix]],ignore_index=True)
            antisense=antisense.append(df3[df3["strand"]=="-"].ix[:,["gene",prefix]],ignore_index=True)
        elif strand=="rev":
            sense=sense.append(df3[df3["strand"]=="-"].ix[:,["gene",prefix]],ignore_index=True)
            antisense=antisense.append(df3[df3["strand"]=="+"].ix[:,["gene",prefix]],ignore_index=True)

    sense.to_csv("%s_sense_counts.xls"%out_file_path,sep="\t",index=False)
    antisense.to_csv("%s_antisense_counts.xls"%out_file_path,sep="\t",index=False)



if __name__ == '__main__':

    (options, args) = prepare_optparser().parse_args()
    main(options)

