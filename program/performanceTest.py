#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 09:53:21 2019

@author: weizhong
"""

import numpy as np
from util import PfamClan,pfamcmp
"""
calculate pfam distance
"""
clanDict = {}
with open('Pfam-A.clans.tsv','r',encoding='utf-8') as fo:
    for line in fo:
        line = line.replace("\n","")
        s = line.split("\t")
        pc = PfamClan(s[1],s[2],s[3],s[4])
        clanDict[s[0]] = pc
        
import scipy.io as sci
pfmat = sci.loadmat('7329pfams_np.mat')
pfs=pfmat['pfs']
pfam_dist = np.zeros((7329,7329))

for i in range(7329):
    for j in range(7329):
        pfam_dist[i,j] = pfamcmp(pfs[0,i][0].split(','),pfs[0,j][0].split(','),clanDict)  

"""
calculate grey incidence degree
"""
from greymodel import GreyIncidenceDegree
gia_dist = np.zeros((7329,7329))
with open('trainPseAAC_7329.npy','rb') as f:
    data = np.load(f)
for i in range(7329):
    gia_dist[i,:] = GreyIncidenceDegree(data[i,:],data)
        
