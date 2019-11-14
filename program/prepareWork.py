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
        
