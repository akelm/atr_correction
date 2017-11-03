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
#    matrix for integration along 0th dimension
    calka=(wavenumbers[:,None] * k[:,None] * delta_waven[:,None] )/mianownik
    calka[np.isinf(calka)]=0
#    result - the refractive index
#    value 1.485 is the shift for toluene taken from "Determination of infrared optical constants for single-component hydrocarbon fuels"
    n=2/np.pi*np.nansum(calka,axis=0) + 1.485
#
    new_range=np.arange(wavenumbers.size)
    new_ind=new_range[sort_ind]
#    print(np.sum(np.abs(wavenumbers_in-wavenumbers[new_ind])))
    return n[new_ind]
    
    
#    function [lambda,wynik,k]=kramers(lambda,k)
#        if size(lambda,2)>size(lambda,1); lambda=transpose(lambda); end
#        if size(k,2)>size(k,1); k=transpose(k); end;
#        omega=1./lambda;
#        % ~~~~~~~~~~
#        %ta sekcja interpoluje dane. aczkolwiek chyba nie jest to niezbedne.
#        m=omega(1);
#        M=omega(end);
#        if m<M; s=5*10^(-6); else s=-10^(-6); end
#        k=interp1(omega,k,transpose(m:s:M),'spline');
#        omega=transpose(m:s:M);
#        % ~~~~~~~~~~
#        d1=abs(omega(2:end)-omega(1:(end-1)));
#        deltaomega=zeros(size(omega));
#        deltaomega(1:(end-1))=deltaomega(1:(end-1))+0.5*d1;
#        deltaomega(2:end)=deltaomega(2:end)+0.5*d1;
#        deltaomega([1,end])=2*deltaomega([1,end]);
#        mianownik=bsxfun(@minus,omega.^2,omega'.^2);
#        calka=bsxfun(@times,omega.*k.*deltaomega,1./mianownik);
#        calka(abs(calka)==Inf)=0;
#        wynik=2/pi*sum(calka,1);
#        lambda=1./omega;
#    end

def abs_exact(n,k,wavenumbers,ncryst=2.4,angle=45):
    theta=angle*np.pi/180
    n2=n+k*1j
    rs12= (ncryst*np.cos(theta) - 1j*np.sqrt(ncryst**2 * np.sin(theta)**2 - n2**2) ) / \
        (ncryst*np.cos(theta) + 1j*np.sqrt(ncryst**2 * np.sin(theta)**2 - n2**2) )
    rp12= (ncryst**2 * np.cos(theta) - 1j*ncryst*np.sqrt(ncryst**2 * np.sin(theta)**2 - n2**2) ) / \
        (ncryst**2 * np.cos(theta) + 1j*ncryst*np.sqrt(ncryst**2 * np.sin(theta)**2 - n2**2) ) 
    Rs=(rs12 * rs12.conj()).real
    Rp=(rp12 * rp12.conj()).real
    R=0.5*(Rs+Rp)
    Aexact=-np.log10(R)
    return Aexact

def atr_exact(n,k,Aexp,wavenumbers,ncryst=2.4,angle=45):
    theta=angle*np.pi/180
    Aexact=abs_exact(n,k,Aexp,wavenumbers)
    delta_A=Aexp-Aexact
    
    
    cc=np.nansum((delta_A**2) / np.abs(Aexp))
#    print(cc)
    it=0
    while cc>0.01:
        delta_k1=delta_A*(ncryst**2-n**2) * np.sqrt(ncryst**2 * np.sin(theta)**2 - n**2) / (0.434 * 3/2 * 4 *n*ncryst*np.cos(theta))
        delta_k=delta_k1
        k=k+delta_k
        n=kramers(k,wavenumbers)
        Aexact=abs_exact(n,k,wavenumbers)
        delta_A=Aexp-Aexact
        cc=np.nansum((delta_A**2) / np.abs(Aexp))
        it+=1
    print('=====exact ',it)
    
    return (n,k)    
    

def atr_lowabs(filename,ncryst=2.4,angle=45):
    

    
    widmo=np.genfromtxt(filename,delimiter=' ')
#    print(widmo)
    Aexp=widmo[:,1]
    wavenumbers=widmo[:,0]
    # ref index of diamond
    #ncryst=2.4
    theta=angle*np.pi/180
    # initial guesses for n and k
    n=1.5*np.ones(Aexp.shape)
#    k_old=0
    k=Aexp * (ncryst**2-n**2) * np.sqrt(ncryst**2 * np.sin(theta)**2 - n**2) / (0.434 * 3/2 * 4 *n*ncryst*np.cos(theta))
    n=kramers(k,wavenumbers)
    Acalc = (0.434 * 3/2 * 4 *n*ncryst*k*np.cos(theta)) / ( (ncryst**2-n**2) * np.sqrt(ncryst**2 * np.sin(theta)**2 - n**2) )
    delta_A=Aexp-Acalc
    cc=np.nansum((delta_A**2) / np.abs(Aexp))
#    print(cc)
#    print(n.shape)
#    print(wavenumbers.shape)
#    sort_ind=np.argsort(wavenumbers)
    
#    imgplot4 = plt.plot(wavenumbers,n,'-k')
#    imgplot5 = plt.plot(wavenumbers,k,'-b')
    it=0
#    plt.plot(wavenumbers,k * wavenumbers,'-k')
#    cc=np.sum(delta_A**2 / Aexp)
    while cc > 0.01:
        k=Aexp * (ncryst**2-n**2) * np.sqrt(np.abs(ncryst**2 * np.sin(theta)**2 - n**2)) / (0.434 * 3/2 * 4 *n*ncryst*np.cos(theta))
        n=kramers(k,wavenumbers)
#        k_old=k.copy()
        Acalc = (0.434 * 3/2 * 4 *n*ncryst*k*np.cos(theta)) / ( (ncryst**2-n**2) * np.sqrt(ncryst**2 * np.sin(theta)**2 - n**2) )
        delta_A=Aexp-Acalc
        cc=np.nansum((delta_A**2) / np.abs(Aexp))
        it+=1
#        plt.plot(wavenumbers,k * wavenumbers,'-r')
    print('=====lowabs ',it)
    
    
    n,k=atr_exact(n,k,Aexp,wavenumbers)
#    this is not molar ext coefficient, but a value proportional to molar ext coefficient
    molar_ext=k * wavenumbers
#    atr signal from exact formulas
    Aexact=abs_exact(n,k,wavenumbers)
    
    
    
#==============================================================================
#     plots
#==============================================================================
    plt.ioff()
    fig = plt.figure(figsize=(8,8*9/16))
    ax = fig.add_subplot(111)
    imgplot1 = plt.plot(wavenumbers,Aexp/np.max(Aexp),'-k')
    imgplot2 = plt.plot(wavenumbers,molar_ext/np.max(molar_ext),'-b')
#    imgplot3 = plt.plot(wavenumbers,Acalc,'-r')
#    imgplot3 = plt.plot(wavenumbers,molar_ext/np.max(molar_ext),'-r')
#    imgplot4 = plt.plot(wavenumbers,(n-1.485)/np.max(n-1.485),'-g')
#    ax.set_title(obj[1])
    ax.set_xlabel('wavenumber, 1/cm')
#    figname='atr'
    plt.savefig(os.path.splitext(filename)[0]+'.svg')
    plt.cla()
#                    print(plt.get_fignums())
    plt.close('all')
    
    np.savetxt(os.path.splitext(filename)[0]+'_ref_ind.txt',np.concatenate( (wavenumbers[:,None],n[:,None]), axis=1) )
    np.savetxt(os.path.splitext(filename)[0]+'_molar_ext.txt',np.concatenate( (wavenumbers[:,None],molar_ext[:,None]), axis=1) )
    
    
#filename='/home/ania/Pulpit/ohpc23/atr/Kuba Ostapko/ohpc23H_new.txt'
#atr_lowabs(filename)