#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 10:59:40 2017

@author: ania
"""
import numpy as np
from scipy import constants
import matplotlib.pyplot as plt
import os
from sklearn.metrics import r2_score


def kramers(k,wavenumbers_in):
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
    delta_waven=np.zeros(wavenumbers.shape)
    delta_waven[0:-1]=delta_waven[0:-1]+d_waven*0.5
    delta_waven[1:]=delta_waven[1:]+d_waven*0.5
    delta_waven[[0,-1]]=2*delta_waven[[0,-1]]
#    denominator
    mianownik=wavenumbers[:,None]**2 - wavenumbers[None,:]**2
    mianownik[np.diag_indices(wavenumbers.size)]=np.finfo(float).eps
#    matrix for integration along 0th dimension
    calka=(wavenumbers[:,None] * k[:,None] * delta_waven[:,None] )/mianownik
    calka[np.diag_indices(wavenumbers.size)]=0
    calka[np.isinf(calka)]=0
#    result - the refractive index
#    value 1.485 is the shift for toluene taken from "Determination of infrared optical constants for single-component hydrocarbon fuels"
    n=2/np.pi*np.nansum(calka,axis=0) + 1.485
#   rearranging back
    new_range=np.arange(wavenumbers.size)
    new_ind=new_range[sort_ind]
    return n[new_ind]
    
    

def abs_exact(n,k,wavenumbers,ncryst,angle,nrefl):
    theta=angle*np.pi/180
    n2=n+k*1j
    rs12= (ncryst*np.cos(theta) - 1j*np.sqrt(ncryst**2 * np.sin(theta)**2 - n2**2) ) / \
        (ncryst*np.cos(theta) + 1j*np.sqrt(ncryst**2 * np.sin(theta)**2 - n2**2) )
    rp12= (ncryst**2 * np.cos(theta) - 1j*ncryst*np.sqrt(ncryst**2 * np.sin(theta)**2 - n2**2) ) / \
        (ncryst**2 * np.cos(theta) + 1j*ncryst*np.sqrt(ncryst**2 * np.sin(theta)**2 - n2**2) ) 
    Rs=(rs12 * rs12.conj()).real
    Rp=(rp12 * rp12.conj()).real
    R=0.5*(Rs**nrefl+Rp**nrefl)
    Aexact=-np.log10(R)
    return Aexact

def plotting(x,ys,filename):
    if type(ys) is not list:
        ys=[ys]
    plt.ioff()
    fig = plt.figure(figsize=(8,8*9/16))
    fig.add_subplot(111)
    colors=['k','b','r','g','m']
    colors=colors*(len(ys)//5) + colors[slice(len(ys)%5)]
    
    for y,color in zip(ys,colors):
       plt.plot(x,y/np.max(y),'-'+color) 
        
    plt.savefig(os.path.splitext(filename)[0]+'.svg')
    plt.cla()
    plt.close('all')
    

def atr_exact(filename,ncryst=2.4,angle=45,nrefl=1,r2=1):
    
    widmo=np.genfromtxt(filename,delimiter=' ')
    Aexp=widmo[:,1]
    Aexp[Aexp==0]=np.finfo(float).eps
    wavenumbers=widmo[:,0]
    theta=angle*np.pi/180
#==============================================================================
#     # initial guesses for n and k
#==============================================================================
    n=1.5*np.ones(Aexp.shape)
    # k from approx formula
    k=Aexp * (ncryst**2-n**2) * np.sqrt(ncryst**2 * np.sin(theta)**2 - n**2) / (0.434 * 3/2 * 4 *n*ncryst*np.cos(theta))/nrefl
    n=kramers(k,wavenumbers)
#==============================================================================
# #    atr signal from exact formulas
#==============================================================================
    Aexact=abs_exact(n,k,wavenumbers,ncryst,angle,nrefl)
    cond = np.abs(1-r2_score(Aexp,Aexact))
    it=0
    while cond>0.01:
        # from approx formula from:
        # Milosevic, Milan. Internal reflection and ATR spectroscopy. Vol. 262. John Wiley &amp; Sons, 2012.
        delta_A=Aexp-Aexact
        delta_k=delta_A*(ncryst**2-n**2) * np.sqrt(ncryst**2 * np.sin(theta)**2 - n**2) / (0.434 * 3/2 * 4 *n*ncryst*np.cos(theta)) / nrefl
        delta_k[np.isnan(delta_k)]=np.finfo(float).eps
        k=k+delta_k
        n=kramers(k,wavenumbers)
        Aexact=abs_exact(n,k,wavenumbers,ncryst,angle,nrefl)
        cond = np.abs(1-r2_score(Aexp,Aexact))
        it+=1
    print('=====exact ',it)

#==============================================================================
# #    outputs    
#==============================================================================
    #    this is not molar ext coefficient, but a value proportional to molar ext coeff
    molar_ext=k * wavenumbers
    plotting(wavenumbers,[Aexp,Aexact,molar_ext],filename)
    np.savetxt(os.path.splitext(filename)[0]+'_ref_ind.txt',np.concatenate( (wavenumbers[:,None],n[:,None]), axis=1) )
    np.savetxt(os.path.splitext(filename)[0]+'_molar_ext.txt',np.concatenate( (wavenumbers[:,None],molar_ext[:,None]), axis=1) )
    
    
