# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 09:59:18 2019

@author: Administrator
"""
"""
 读取全部的Pfam的家族信息
 形如：PF00001	CL0192	GPCR_A	7tm_1	7 transmembrane receptor (rhodopsin family)
 构建一个字典，key是Pfam的编号，后面clan的ID，name，abb和描述作为一个结构体PfamClan
"""
"""
计算测试集每个样本与训练集样本的pfam集合距离
2019-11-14
"""
"""
from util import PfamClan
import scipy.io as sio
clanDict = {}

with open('Pfam-A.clans.tsv','r',encoding='utf-8') as fo:
    for line in fo:
        line = line.replace("\n","")
        s = line.split("\t")
        pc = PfamClan(s[1],s[2],s[3],s[4])
        clanDict[s[0]] = pc


from util import pfamcmp
testdata_pfams=sio.loadmat('independent_pfams.mat')
traindata_pfams=sio.loadmat('train_pfams.mat')
pfam_dist={}
testkeys = list(testdata_pfams.keys())
trainkeys = list(traindata_pfams.keys())
totalTest = len(testkeys[3:])
totalTrain = len(trainkeys[3:])
count = 1
for k1 in testkeys[3:]:
    print("\n        {}/{}==>{:.2f}%".format(count,totalTest,count/totalTest * 100),end="")
    dist=[]
    pf1=testdata_pfams[k1]
    
    n = 1
    for k2 in trainkeys[3:]:
        print("\r{:.2f}%".format(n/totalTrain * 100),end="")
        pf2=traindata_pfams[k2]
        d=pfamcmp(pf1,pf2,clanDict)
        dist.append(d)
        n = n + 1
    pfam_dist[k1] = dist
    count = count + 1
"""

"""
基于PSSM矩阵计算测试集每个样本与训练集的灰色关联度
"""  
"""     
from util import readPSSMFromFile
from greymodel import greyPsePSSM
import os
import numpy as np
import scipy.io as sio
dirname = r'E:\Repoes\jcilwz\RemoteHomology\program\scope_independent_PSSM'
files = os.listdir(dirname)
psePssmAAC = {}
for file in files:
    pid = file[:-4]
    path = os.path.join(dirname, file)
    pssm = readPSSMFromFile(path)
    psePssmAAC[file[:-4]] = greyPsePSSM(np.array(pssm,dtype=np.float), model=2)
sio.savemat('independent_psePssmAAC.mat',psePssmAAC)
"""

"""
读入训练集家族类别
"""
"""
from Bio import SeqIO
import os
import scipy.io as sio 
files = os.listdir('SCOP167-superfamily')
posFiles = []
for file in files:
    if file.startswith('pos-train'):
        posFiles.append(file)
train_records = list(SeqIO.parse('SCOP167_pos_train.fa','fasta'))
#family = {}
family = sio.loadmat('train_family_15000.mat')
k = 0
totalTrain = 145187
for record in train_records:
    fam = []
    k = k + 1
    if k <= 15000:
        continue
    for file in posFiles:
        name = file[10:-6]
        path = os.path.join('SCOP167-superfamily',file)
        for r in SeqIO.parse(path,'fasta'):
            if str(record.seq) == str(r.seq):
                fam.append(name)
    family[str(k)] = fam
    if k%5000 == 0:
        sio.savemat("".join(['train_family_', str(k), '.mat']),family)
    print("\n        {}/{}==>{:.2f}%".format(k, totalTrain, k/totalTrain * 100),end="")        
sio.savemat('train_family.mat',family) 
"""
 
"""
计算序列SCOP167 Dataset greyPseAAC
"""
"""
from SeqFormulate import greyPseAAC
from Bio import SeqIO
import numpy as np
# training data

trainFile = 'SCOP167_pos_train.fa'
codeTypes = ['MolecularWeight','Hydrophobicity','PK1','PK2','PI']
PseAAC = []
k = 1
totalTrain = 145187
for seq_record in SeqIO.parse(trainFile,'fasta'):
     PseAAC.append(greyPseAAC(str(seq_record.seq),codeTypes)) 
     print("\r{}/{}==>{:.2f}%".format(k, totalTrain, k/totalTrain * 100),end="")
     k += 1
train_pseAAC = np.array(PseAAC)
with open('trainPseAAC.npy','wb') as f:
    np.save(f, train_pseAAC)  

# testing data
testFile = 'scope_independent.fa'
#codeTypes = ['MolecularWeight','Hydrophobicity','PK1','PK2','PI']
PseAAC = []
k = 1
for seq_record in SeqIO.parse(testFile,'fasta'):
     PseAAC.append(greyPseAAC(str(seq_record.seq),codeTypes))
test_pseAAC = np.array(PseAAC)
with open('testPseAAC.npy','wb') as f:
    np.save(f, test_pseAAC)  
"""

"""
calculate 7329 dataset greyPseAAC GM(1,1)
"""
from SeqFormulate import greyPseAAC
from Bio import SeqIO
import numpy as np
trainFile = '7329seqs.fasta'
codeTypes = ['MolecularWeight','Hydrophobicity','PK1','PK2','PI']
PseAAC = []
k = 1
totalTrain = 7329
for seq_record in SeqIO.parse(trainFile,'fasta'):
     PseAAC.append(greyPseAAC(str(seq_record.seq),codeTypes)) 
     print("\r{}/{}==>{:.2f}%".format(k, totalTrain, k/totalTrain * 100),end="")
     k += 1
train_pseAAC = np.array(PseAAC)
with open('trainPseAAC_7329.npy','wb') as f:
    np.save(f, train_pseAAC)  