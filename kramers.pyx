#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 18:50:16 2017

@author: ania
"""

import numpy as np
cimport numpy as np

#cdef int i=1

def kramers(np.ndarray[np.float64_t] k, wavenumbers_in):
    cdef np.ndarray[np.float64_t] wavenumbers, d_waven, delta_waven, n
    cdef np.ndarray[np.int64_t] sort_ind, new_range, new_ind
    cdef np.ndarray[np.float64_t,ndim=1] mianownik, calka
    cdef int l
#    k - imaginary part of refractive index
#    wavenumbers - in cm-1
#
#   sorting data
    sort_ind=np.argsort(wavenumbers_in)
    wavenumbers=wavenumbers_in[sort_ind]
#    print(np.sum(np.abs(wavenumbers_in-wavenumbers[sort_ind])))
    k=k[sort_ind]
#   vector with dataspacing
    d_waven=wavenumbers[1:]-wavenumbers[0:-1]
#    vector for integration
    delta_waven=np.zeros(wavenumbers.size)
    delta_waven[0:-1]=delta_waven[0:-1]+d_waven*0.5
    delta_waven[1:]=delta_waven[1:]+d_waven*0.5
    delta_waven[[0,-1]]=2*delta_waven[[0,-1]]
#==============================================================================
#     stare, to sie zmieni
#==============================================================================
##    denominator
#    mianownik=wavenumbers[:,None]**2 - wavenumbers[None,:]**2
#    mianownik[np.diag_indices(wavenumbers.size)]=np.finfo(float).eps
##    matrix for integration along 0th dimension
#    calka=(wavenumbers[:,None] * k[:,None] * delta_waven[:,None] )/mianownik
#    calka[np.diag_indices(wavenumbers.size)]=0
#    calka[np.isinf(calka)]=0
##    result - the refractive index
##    value 1.485 is the shift for toluene taken from "Determination of infrared optical constants for single-component hydrocarbon fuels"
#    n=2/np.pi*np.nansum(calka,axis=0) + 1.485
#==============================================================================
#     koniec starego
#==============================================================================
#==============================================================================
# poczatek zmian
#==============================================================================
    n=np.zeros(wavenumbers.size)
#    denominator
    for l in range(wavenumbers.size):
        mianownik=wavenumbers**2 - wavenumbers[l]**2
        mianownik[l]=np.finfo(float).eps
    #    matrix for integration along 0th dimension
        calka=(wavenumbers * k * delta_waven )/mianownik
        calka[l]=0
        calka[np.isinf(calka)]=0
        #    result - the refractive index
        #    value 1.485 is the shift for toluene taken from "Determination of infrared     optical constants for single-component hydrocarbon fuels"
        calka[calka==np.NaN]=0
        n[l]=2/np.pi*np.sum(calka) + 1.485    
#==============================================================================
# koniec zmian
#==============================================================================
#   rearranging back
    new_range=np.arange(wavenumbers.size,dtype=np.int64)
    new_ind=new_range[sort_ind]
    return n[new_ind]