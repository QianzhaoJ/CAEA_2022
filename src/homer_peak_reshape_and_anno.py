#!/work2/liugh/liuzunpeng/04_Softwares/bin/python3
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
    optparser.add_option('-s', '--strand', dest='strand', default='+',help="\nGiven the thread of multiBamSummary.")
    optparser.add_option('-a', '--anno', dest='anno', default='/work2/liugh/liuzunpeng/04_Softwares/homer/bin/annotatePeaks.pl', help="\nGiven the path of annotatePeaks.pl")
    optparser.add_option('-t', '--state', dest='state',
                         default='/work2/liugh/liuzunpeng/04_Softwares/ATAC_seq_src/macs_stat.pl',
                         help="\nGiven the path of annotatePeaks.pl'.")
    optparser.add_option("-h", '--help', action="help",help="\nShow this help message and exit.")
    return optparser


def main(opt):
    inpeak    = opt.input
    out_path  = opt.out
    prefix    = opt.prefix
    out_peak    = out_path + '/' + prefix + '.homer.peaks.bed'
    out_peak_f  = open(out_peak,'w')
    out_peak5   = out_path + '/' + prefix + '.homer.peaks_enrich5.bed'
    out_peak5_f  = open(out_peak5,'w')
    #out_peak10  = out_path + '/' + prefix + '.homer.peaks_enrich10.bed'
    #out_peak10_f  = open(out_peak10,'w')
    #out_peak30  = out_path + '/' + prefix + '.homer.peaks_enrich30.bed'
    #out_peak30_f  = open(out_peak30,'w')
    annotatePeaks = opt.anno
    state = opt.state

    state_txt   = out_path + '/' + prefix + '.homer.peaks.state.txt'
    state_txt5  = out_path + '/' + prefix + '.homer.peaks.state_enrich5.txt'
    #state_txt10 = out_path + '/' + prefix + '.homer.peaks.state_enrich10.txt'
    #state_txt30 = out_path + '/' + prefix + '.homer.peaks.state_enrich30.txt'

    anno_txt    = out_path + '/' + prefix + '.homer.peaks.anno.txt'
    anno_enrich_txt  = out_path + '/' + prefix + '.homer.peaks.anno.enrich.txt'

    strand=opt.strand

    chrs=['chr'+str(x) for x in range(1,23)]+['chrX']

    with open(inpeak,'r') as inf:
        for lines in inf:
            if lines.startswith("chr"):
                lines=lines.strip()
                line=lines.split('\t')
                if line[1] in chrs:
                    #out_peak_f.write('chr\tstart\tend\tpeak_name\tenrichment\tstrand\n')
                    #out_peak5_f.write('chr\tstart\tend\tpeak_name\tenrichment\tstrand\n')
                    #out_peak10_f.write('chr\tstart\tend\tpeak_name\tenrichment\tstrand\n')
                    #out_peak30_f.write('chr\tstart\tend\tpeak_name\tenrichment\tstrand\n')
                    s=line[1]+'\t'+line[2]+'\t'+line[3]+'\t'+line[0]+'\t'+line[10]+'\t'+strand+'\n'
                    out_peak_f.write(s)
                    if float(line[10]) >= 5:
                        out_peak5_f.write(s)
                    #if float(line[10]) >= 10:
                    #    out_peak10_f.write(s)
                    #if float(line[10]) >= 30:
                    #    out_peak30_f.write(s)

    out_peak_f.close()
    out_peak5_f.close()
    #out_peak10_f.close()
    #out_peak30_f.close()

    os.system('perl '+ state + ' ' + out_peak   + ' > ' + state_txt)
    os.system('perl '+ state + ' ' + out_peak5  + ' > ' + state_txt5)
    #os.system('perl '+ state + ' ' + out_peak10 + ' > ' + state_txt10)
    #os.system('perl '+ state + ' ' + out_peak30 + ' > ' + state_txt30)
    os.system(annotatePeaks + ' ' + out_peak + ' hg19 -annStats ' + anno_enrich_txt  + ' > ' + anno_txt)


if __name__ == '__main__':

    (options, args) = prepare_optparser().parse_args()
    main(options)
