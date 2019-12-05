#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 09:53:21 2019

@author: weizhong
"""

import numpy as np
from util import PfamClan,pfamcmp
import scipy.io as sci
from SeqFormulate import greyPseAAC
#calculate pfam caln
clanDict = {}
with open('Pfam-A.clans.tsv','r',encoding='utf-8') as fo:
    for line in fo:
        line = line.replace("\n","")
        s = line.split("\t")
        pc = PfamClan(s[1],s[2],s[3],s[4])
        clanDict[s[0]] = pc
        

pfm_dist= sci.loadmat('pfam_dist.mat')

#calculate distance
from Bio import SeqIO
test_records = list(SeqIO.parse('scope_independent.fa','fasta'))
num_test = len(test_records)
codeTypes = ['MolecularWeight','Hydrophobicity','PK1','PK2','PI']
test_pseAAC = []
for seq_record in test_records:
    test_pseAAC.append(greyPseAAC(str(seq_record.seq), codeTypes))
test_pseAAC = np.array(test_pseAAC)

train_pseAAC = np.load('trainPseAAC.npy')
num_train = 145187
eucl_dist = np.zeros((num_test, num_train))

print("calculating euclidean distance between test samples and train samples......")
for i  in range(num_test):
    print("\r{}/{}==>{:.2%}".format(i, num_test, i/num_test), end="")
    for j in range(num_train):
        eucl_dist[i,j] = np.linalg.norm(test_pseAAC[i] - train_pseAAC[j])
        
with open('scope_test_train_euclidean_distance.npy','wb') as f:
    np.save(f, eucl_dist)
    
# jack knife test   
w = [0.5, 0.5]
count = 0
familymat = sci.loadmat('train_family.mat')

print("calculating prediction accuracy......")  
for i in range(num_test):
    d = (1-pfm_dist[i]) * w[0] + eucl_dist[i] * w[1]
    indx_arr = np.argsort(d)
    seq = test_records[i]
    desc = seq.descriptioin
    faml = (desc.split(" "))[1]

    if familymat[str(indx_arr[1])][0] == faml:
        count += 1
    print("\r{}/{}==>{:.2%}".format(i, num_test, i/num_test), end="")
print("predicted result: in {0} samples corrected predict {1}: {2:.2%}".format(num_test,count,count/num_test))    


