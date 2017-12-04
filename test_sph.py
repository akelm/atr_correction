#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 15:11:21 2017

@author: ania
"""

import numpy as np
#cimport numpy as np
from scipy import special
#cimport scipy.special.cython_special as special
#cimport cython


nmax=5
z=np.expand_dims(1j*np.arange(0.5,1.2,0.1)+np.arange(0.3,1,0.1),axis=1)

n=np.expand_dims(np.arange(1,nmax+1),axis=1).swapaxes(0,1)

yn=special.spherical_yn(n, z, 0)
#print(yn.shape)
#print(yn)
yv=np.sqrt(np.pi/2/z)*special.yv(n+0.5,z)
#print(yv.shape)
print(yv-yn)