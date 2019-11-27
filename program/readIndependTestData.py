# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 10:27:19 2019

@author: Administrator
"""
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
scopeSeqDict={}
for record in SeqIO.parse(r'E:\Repoes\jcilwz\RemoteHomology\program\independ_test.fa','fasta'):
    scopeSeqDict[record.id]=str(record.seq).upper()

file = r'E:\Repoes\jcilwz\RemoteHomology\program\independ_test.txt'
protFamilyDict = {}
protSeqDict = {}
flag = False

with open(file,'r') as fo:
    for line in fo.readlines():
        if len(line) == 1:
            continue
        if line.startswith('Family'):
            k = line.index('(')
            j = line.index(')')
            family = line[k+1:j]
            continue
        elif line.startswith("Positive"):
            flag = True
            continue
        elif line.startswith("Negative"):
            flag = False
            continue
        if flag:
            prots = line.split()
            for p in prots:
                protSeqDict[p] = scopeSeqDict[p]
                protFamilyDict[p]=family

prots = []
for p in protFamilyDict.keys():
    record = SeqRecord(Seq(protSeqDict[p]), id=p, description=protFamilyDict[p])
    prots.append(record)
SeqIO.write(prots,'scope_independent.fa','fasta')      
'''
trainProtCount={}
trainFile = r'E:\Repoes\jcilwz\RemoteHomology\program\SCOP167_pos_train.fa'
from Bio import SeqIO
for record in SeqIO.parse(trainFile,'fasta'):
    trainProtCount[record.id] = trainProtCount.get(record.id,0)+1
items = list(trainProtCount.items())
items.sort(key=lambda x:x[1],reverse=True)    
'''