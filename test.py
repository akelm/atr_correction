#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 19 09:34:22 2017

@author: ania
"""

import time
import numpy as np
from atr import kramers as kramers_pyt
import kramers
import cProfile
#import sys

#sys.path.insert(1,'/home/ania/Desktop/atr_correction')
#print(sys.path)

size=8000
wavenumbers=np.arange(size,dtype=np.float64)+1
k=np.random.rand(size)

#t0=time.time()
#for l in range(5):
cProfile.run('kramers_pyt(k,wavenumbers)', sort='time')
#t1=time.time()
#print('python ',t1-t0)

#t2=time.time()
#for l in range(5):
cProfile.run('kramers.kramers(k,wavenumbers)', sort='time')
#t3=time.time()
#print('cython ',(t3-t2)/(t1-t0)*100)