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
from greymodel import GreyIncidenceDegree
# load trains family
mat_family = sci.loadmat(r'E:\360yunpan\MyProgram\RemoteHomo\3187sequencesFamily.mat')
familyId = mat_family['familyId']

# load trains pfam distance 
mat_dist= sci.loadmat(r'E:\360yunpan\MyProgram\RemoteHomo\3187FunctionDomainMatch_0.2.mat')
pfam_dist = mat_dist['v']

#calculate euclidean distance
from Bio import SeqIO
train_records = list(SeqIO.parse('3187seqs.fasta','fasta'))
num_train = len(train_records)
codeTypes = ['MolecularWeight','Hydrophobicity','PK1','PK2','PI']
train_pseAAC = []
weight = np.ones((30,))
weight[20:30:2] = 40
weight[21:30:2] = 0.1
for seq_record in train_records:
    v = greyPseAAC(str(seq_record.seq), codeTypes) * weight
    train_pseAAC.append(v)
train_pseAAC = np.array(train_pseAAC)

#eucl_dist = np.zeros((num_train, num_train))
gid_dist = np.zeros((num_train, num_train))

print("calculating grey incidence degree distance between train samples......")
for i  in range(num_train):
    print("\r{}/{}==>{:.2%}".format(i, num_train, i/num_train), end="")
    gid_dist[i] = GreyIncidenceDegree(train_pseAAC[i], train_pseAAC)
    #for j in range(num_train):
        #eucl_dist[i,j] = np.linalg.norm(train_pseAAC[i] * weight- train_pseAAC[j] * weight)
        
#with open('3187_euclidean_weight_distance.npy','wb') as f:
#   np.save(f, eucl_dist)
    
# jack knife test   
w = [1, 10]
count = 0
print("calculating prediction accuracy......")  
for i in range(num_train):
    seq = train_records[i]
    desc = seq.description
    head = desc.split(" ")
    sid = head[0]
    faml = head[1]
    #d = (1-pfam_dist[i]) * w[0] + eucl_dist[i] * w[1]
    d = pfam_dist[i] * w[0] + gid_dist[i] * w[1]
    indx_arr = np.argsort(d)   

    if familyId[0,indx_arr[-2]] == faml:
        count += 1
    print("\r{}/{}==>{:.2%}".format(i, num_train, i/num_train), end="")
print("predicted result: in {0} samples corrected predict {1}: {2:.2%}".format(num_train,count,count/num_train))    


