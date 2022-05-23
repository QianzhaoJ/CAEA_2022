#!/usr/bin/python
#coding:utf-8

import os
import re
import sys

input=sys.argv[1]
lsd = os.listdir(input)
lsd.sort()
lst=[]

py='/work2/liugh/liuzunpeng/01_Script/01_pipeline/02_ChIP-seq/00_ChIP_seq_pipeline/30_Lamin_KI_DamID/DamID.seq.py'
ii=input
out=re.findall(r'(.*/)',py)[0]+'step1.qsub.sh'
oo='/work2/liugh/liuzunpeng/05_Results/59_Lamin_KI_DamID'
cc='/work2/liugh/liuzunpeng/01_Script/01_pipeline/02_ChIP-seq/00_ChIP_seq_pipeline/config.txt'
qq='/work2/liugh/liuzunpeng/06_Qsub/59_Lamin_KI_DamID'

for i in lsd:
    pre=re.findall(r'(.*)_1.fq.gz',i)
    if len(pre)==1:
        lst.append(pre[0])

f1=open(out,'w')
for x in lst:
    tempstr='python3 '+py+' -i '+input+' -o '+oo+' -q '+qq+' -c '+cc+' -p '+x+'\n'
    f1.write(tempstr)

f1.close()

