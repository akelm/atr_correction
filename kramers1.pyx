#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 14:18:09 2017

@author: ania
"""

import numpy as np
cimport numpy as np
cimport cython

#cdef int i=1
#print("running kramers1")
@cython.boundscheck(False)
cdef np.ndarray[np.float64_t] kramers(np.ndarray[np.float64_t] k_in, np.ndarray[np.float64_t] wavenumbers_in):
    cdef np.ndarray[np.float64_t] wavenumbers, delta_waven, n, licznik
    cdef np.ndarray[np.int64_t] sort_ind
    cdef np.ndarray[np.float64_t,ndim=1] mianownik, calka
    cdef int l, wavenumbers_l
#    k - imaginary part of refractive index
#    wavenumbers - in cm-1
#
#    print("running kramers1")
    wavenumbers_l=wavenumbers_in.size
#   sorting data
    sort_ind=np.argsort(wavenumbers_in)
    wavenumbers=np.zeros(wavenumbers_l)
    k=np.zeros(wavenumbers_l)
    for l in range(wavenumbers_l):
        wavenumbers[l]=wavenumbers_in[sort_ind[l]]
        k[l]=k_in[sort_ind[l]]      

#   vector with dataspacing
    delta_waven=np.zeros(wavenumbers_l)
    for l in range(1,wavenumbers_l-1):
        delta_waven[l]=0.5*(wavenumbers[l+1]-wavenumbers[l-1])
    delta_waven[0]=wavenumbers[1]-wavenumbers[0]
    delta_waven[wavenumbers_l-1]=wavenumbers[wavenumbers_l-1]-wavenumbers[wavenumbers_l-2]
    licznik=wavenumbers * k * delta_waven
#    denominator
    n=np.zeros(wavenumbers_l)
    for l in range(wavenumbers_l):
        mianownik=wavenumbers**2 - wavenumbers[l]**2
        mianownik[l]=1
    #    matrix for integration along 0th dimension
        calka=licznik/mianownik
        calka[l]=0
#        calka[np.isinf(calka)]=0
        #    result - the refractive index
        #    value 1.485 is the shift for toluene taken from "Determination of infrared     optical constants for single-component hydrocarbon fuels"
#        calka[calka==np.NaN]=0
#        with rearangement
        n[sort_ind[l]]=2/np.pi*np.sum(calka) + 1.485    
#==============================================================================
# koniec zmian
#==============================================================================
#   rearranging back
#    n_out=np.zeros(wavenumbers_l)
#    for l in range(wavenumbers_l):
#        n_out[sort_ind[l]]=n[l]        
    return n