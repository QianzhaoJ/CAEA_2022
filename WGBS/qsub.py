#!/usr/bin/python
#coding:utf-8

import os
import re
import sys

input=sys.argv[1]
lsd = os.listdir(input)
lsd.sort()
lst=[]

py='/work2/liugh/liuzunpeng/01_Script/01_pipeline/09_WGBS/02_progerin/WGBS.v3.py'
ii=input
out=re.findall(r'(.*/)',py)[0]+'step1.qsub.sh'
oo='/work1/liugh/liuzunpeng/05_Results/33_WGBS_progerin'
cc='/work2/liugh/liuzunpeng/01_Script/01_pipeline/09_WGBS/02_progerin/config.txt'
qq='/work1/liugh/liuzunpeng/06_Qsub/33_WGBS_progerin'

for i in lsd:
    if i.endswith('fq.gz'):
        pre=re.findall(r'(.*)_[12].clean.fq.gz',i)
        lst.append(pre[0])
uni_lst=list(set(lst))

f1=open(out,'w')
for x in uni_lst:
    tempstr='python3 '+py+' -i '+input+' -o '+oo+' -q '+qq+' -c '+cc+' -p '+x+'\n'
    f1.write(tempstr)

f1.close()


