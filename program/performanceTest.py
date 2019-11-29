#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 09:53:21 2019

@author: weizhong
"""

import numpy as np
with open('trainPseAAC_7329.npy','rb') as f:
    data = np.load(f)
import scipy.io as sci
pfmat = sci.loadmat('7329pfams_np.mat')
pfs=pfmat['pfs']

pfam_dist = np.zeros((7329,7329))
for i in range(7329):
    for j in range(7329):
        
