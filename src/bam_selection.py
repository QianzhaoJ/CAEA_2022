#!/usr/bin/env python3
#-*-coding:utf-8-*-

from __future__ import division

__author__ = "Zunpeng Liu"
__copyright__ = "Copyright 2019, IBP, CAS."
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Zunpeng Liu"
__email__ = "zunpengliuAT163.com"



from optparse import OptionParser
import os,sys,re
import random

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
    optparser.add_option('-o', '--out', dest='out',default='',help="\nGiven the bed file path.")
    optparser.add_option('-p', '--prefix', dest='prefix', default=' ',help="\nGiven the out directory of out fa path.")
    optparser.add_option('-s', '--samtools', dest='samtools', default='/work2/liugh/liuzunpeng/04_Softwares/bin/samtools',help="\nGiven the out directory of out fa path.")
    optparser.add_option('-n', '--num', dest='num', default=' ',
                         help="\nGiven the number of reads for selection.")
    optparser.add_option("-h", '--help', action="help",help="\nShow this help message and exit.")
    return optparser


def main(opt):

    inbam     = opt.input
    prefix    = opt.prefix
    out       = opt.out
    samtools  = opt.samtools
    num       = opt.num

    sam                   = out+'/'+prefix+'.tmp.sam'
    merge_sam_select_p1   = out+'/'+prefix+'.merge.select.p1.sam'
    merge_sam_select_p2   = out+'/'+prefix+'.merge.select.p2.sam'
    merge_sam_select_p2_o = open(out+'/'+prefix+'.merge.select.p2.sam','w')
    merge_sam_select  = out+'/'+prefix+'.merge.select.sam'
    merge_bam_select  = out+'/'+prefix+'.merge.select.bam'
    bam_selected_sort = out+'/'+prefix + '.merge.bam'
    flagstat_txt      = out + '/' + prefix + '.merge.flagstat.txt'

    os.system(samtools + ' view -@ 24 -H ' + inbam + ' > ' + merge_sam_select_p1)

    os.system(samtools +' view -@ 24 ' + inbam + ' > ' + sam )

    print(' bamtosam done! ')

    reads={}

    with open(sam,'r') as s:
        for lines in s:
            lines = lines.strip()
            line=lines.split('\t')
            if line[0] not in reads.keys():
                reads[line[0]] = [lines]
            else:
                reads[line[0]].append(lines)

    random.seed(1)

    ID = random.sample(reads.keys(),int(int(num)/2))

    for k in ID:
        merge_sam_select_p2_o.write('\n'.join(reads[k])+'\n')

    merge_sam_select_p2_o.close()
    print(' selection done! ')


    os.system('cat ' + merge_sam_select_p1 + ' ' + merge_sam_select_p2 + ' > ' + merge_sam_select)
    os.system(samtools +' view -@ 24 -bS -1 ' + merge_sam_select + ' > ' + merge_bam_select)
    os.system(samtools +' sort -@ 24 -l 9 ' + merge_bam_select + ' -o ' + bam_selected_sort)
    print(' samtobam done! ')


    os.system(samtools + ' index -@ 24 ' + bam_selected_sort)
    os.system(samtools + ' flagstat -@ 24  ' + bam_selected_sort + ' > ' + flagstat_txt)


    #os.system('rm -rf '+ sam)
    #os.system('rm -rf '+ merge_sam_select)


if __name__ == '__main__':

    (options, args) = prepare_optparser().parse_args()
    main(options)


