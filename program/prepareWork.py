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
clanDict = {}

with open(r'e:\360Downloads\Pfam-A.clans.tsv','r',encoding='utf-8') as fo:
    for line in fo:
        line = line.replace("\n","")
        s = line.split("\t")
        pc = PfamClan(s[1],s[2],s[3],s[4])
        clanDict[s[0]] = pc

import scipy.io as sio
sio.savemat('pfamClan.mat', clanDict)

