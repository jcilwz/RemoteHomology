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
with open('trainPseAAC_7329.npy','rb') as f:
    data = np.load(f)
"""
from greymodel import GreyIncidenceDegree
gia_dist = np.zeros((7329,7329))

for j in range(20,30,2):
    data[:,j] = data[:,j] * 50
    data[:,j+1] = data[:,j+1] * 0.1   
for i in range(7329):
    gia_dist[i,:] = GreyIncidenceDegree(data[i,:],data)
np.savez('7329_pfam_weightedgrePseAAC_dist',pfam=pfam_dist, gia=gia_dist)

# jack knife test   
w = [0.5, 0.5]
count = 0
familymat = sci.loadmat('7329sequencesFamily.mat')
fms = familymat['familyId']       
for i in range(7329):
    d = pfam_dist[i] * w[0] + gia_dist[i] * w[1]
    indx_arr = np.argsort(d)
    
    if fms[0,indx_arr[-2]][0] == fms[0,i][0]:
        count += 1
print("predicted result: in {0} samples corrected predict {1}: {2:.2%}".format(7329,count,count/7329))    
"""
eucl_dist = np.zeros((7329,7329))
for i in range(7329):
    print("\r Calculating Eucl dist {.2f}%".format(i/7329 *100), end="")
    for j in range(7329):
        eucl_dist[i,j] = np.linalg.norm(data[i] - data[j])
np.savez('7329_weightedgrePseAAC_eucldist',pfam=pfam_dist, eucl=eucl_dist)

# jack knife test
w = [0.5, 0.5]
count = 0
familymat = sci.loadmat('7329sequencesFamily.mat')
fms = familymat['familyId']       
for i in range(7329):
    print("\r Testing {.2f}%".format(i/7329 *100), end="")
    d = (1-pfam_dist[i]) * w[0] + eucl_dist[i] * w[1]
    indx_arr = np.argsort(d)
    
    if fms[0,indx_arr[1]][0] == fms[0,i][0]:
        count += 1
print("predicted result: in {0} samples corrected predict {1}: {2:.2%}".format(7329,count,count/7329))    