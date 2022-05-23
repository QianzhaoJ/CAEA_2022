#!/usr/bin/env python3
# coding: utf-8

import sys
import re
import os

in_bam=sys.argv[1]
sample=sys.argv[2]

fw_bam1=re.findall(r'(.*)/.*',in_bam)[0]+'/'+sample+'.fw1.bam'
fw_bam2=re.findall(r'(.*)/.*',in_bam)[0]+'/'+sample+'.fw2.bam'
fw_bam=re.findall(r'(.*)/.*',in_bam)[0]+'/'+sample+'.fw.bam'

rev_bam1=re.findall(r'(.*)/.*',in_bam)[0]+'/'+sample+'.rev1.bam'
rev_bam2=re.findall(r'(.*)/.*',in_bam)[0]+'/'+sample+'.rev2.bam'
rev_bam=re.findall(r'(.*)/.*',in_bam)[0]+'/'+sample+'.rev.bam'

samtools="/work2/liugh/liuzunpeng/04_Softwares/samtools/samtools-1.6/samtools"

# include reads that are 2nd in a pair (128);
# exclude reads that are mapped to the reverse strand (16)
str1=samtools+' view -@ 24 -b -f 128 -F 16 '+in_bam+' > '+fw_bam1

# exclude reads that are mapped to the reverse strand (16) and
# first in a pair (64): 64 + 16 = 80
str2=samtools+' view -@ 24 -b -f 80 '+in_bam+' > '+fw_bam2

# combine the temporary files
str3=samtools+' merge -@ 24 -f '+fw_bam+' '+fw_bam1+' '+fw_bam2

os.system(str1)
os.system(str2)
os.system(str3)

# include reads that map to the reverse strand (128)
# and are second in a pair (16): 128 + 16 = 144
str4=samtools+' view -@ 24 -b -f 144 '+in_bam+' > '+rev_bam1

# include reads that are first in a pair (64), but
# exclude those ones that map to the reverse strand (16)
str5=samtools+' view -@ 24 -b -f 64 -F 16 '+in_bam+' > '+rev_bam2

# merge the temporary files
str6=samtools+' merge -@ 24 -f '+rev_bam+' '+rev_bam1+' '+rev_bam2

os.system(str4)
os.system(str5)
os.system(str6)

str7='rm '+fw_bam1+' && rm '+fw_bam2+' && rm '+rev_bam1+' && rm '+rev_bam2

os.system(str7)




