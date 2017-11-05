#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 14:12:48 2017

@author: ania
"""
from atr import atr_lowabs
import numpy as np
import os
import glob
import matplotlib.pyplot as plt

def corr_spctr(filename):
#    rozszerza widmo do zakresu 1:4500 i odejmuje tlo
#    dodatkowo usuwa artefakty z 1700-2300 cm-1
    widmo=np.genfromtxt(filename,delimiter=' ')
#==============================================================================
#     plots
#==============================================================================
    plt.ioff()
    fig = plt.figure(figsize=(8,8*9/16))
    ax = fig.add_subplot(111)
    imgplot1 = plt.plot(widmo[:,0],widmo[:,1],'-k')
#==============================================================================
# ##    wyciecie tla
#==============================================================================
    # constant term
    slice1=np.array([*range(170),*range(1460,1529)])
    widmo[:,1]=widmo[:,1]-np.min(widmo[slice1,1])
    # square term
    slice1=np.array([*range(200),*range(1479,1529)])
    poly1=np.polyfit(slice1,widmo[slice1,1],2)
    imgplot4=plt.plot(widmo[:,0],np.polyval(poly1,np.arange(widmo[:,0].size)),'-b')
    if poly1[-3] >= 0:
        widmo[:,1]=widmo[:,1]-np.polyval(poly1,np.arange(widmo[:,1].size))
    else:
        slice1=np.arange(11)
        slice2=np.arange(160,171)
        slice3=np.arange(1460,1529)
        poly1=np.polyfit([5,165,1563],[np.mean(widmo[slice1,1]),np.mean(widmo[slice2,1]),np.mean(widmo[slice3,1])],2)
        imgplot4=plt.plot(widmo[:,0],np.polyval(poly1,np.arange(widmo[:,0].size)),'-b')
        if poly1[-3] >= 0:
            widmo[:,1]=widmo[:,1]-np.polyval(poly1,np.arange(widmo[:,1].size))
            
#==============================================================================
# #    usuniecie srodka z ujemnymi wartosciami
#==============================================================================

#    slice2=np.array([*range(1079,1124),*range(1489,1529)])
#    poly2=np.polyfit(slice2,widmo[slice2,1],2)
#    zakres=np.arange(1079,1529)
#    widmo[zakres,1]=np.polyval(poly2,zakres)
    
#    widmo(:,2)=widmo(:,2)-widmo(1,2);
#    b_vect=[[4500:-1:4000]',zeros(501,1)];
#    zakres1=[widmo(end,1),1];
#    poly1=polyfit(zakres1,[widmo(end,2),0],2);
#    e_vect=[[398:-1:1]', polyval(poly1,zakres1(1):-1:zakres1(end))'];

#==============================================================================
# #    dorobienie czesci pow 4000 cm-1
#==============================================================================
    b_vect=np.ones([501,2])*np.mean(widmo[0:11,1])
    b_vect[:,0]=np.arange(4500,3999,-1)
#==============================================================================
# #    dorobienie czesci od 1 do 398 cm-1
#==============================================================================
    poly1=np.polyfit(np.array([widmo[-1,0],1]) ,\
                              np.array([widmo[-1,1],0]) , 2)
    e_vect=np.zeros([398,2])
    e_vect[:,0]=np.arange(398,0,-1)
    e_vect[:,1]=np.polyval(poly1,np.arange(398,0,-1))
#==============================================================================
#   full vector
#==============================================================================
    new=np.concatenate( (b_vect,widmo,e_vect),axis=0)
#==============================================================================
#     saving files
#==============================================================================
    filename_new=os.path.splitext(filename)[0]+'_corr'+os.path.splitext(filename)[1]
    np.savetxt(filename_new,new)
#==============================================================================
#     plots
#==============================================================================
    imgplot2 = plt.plot(new[:,0],new[:,1],'-r')
    ax.set_xlabel('wavenumber, 1/cm')
    plt.savefig(os.path.splitext(filename)[0]+'.svg')
    plt.cla()
    plt.close('all')    


files=glob.glob('/home/ania/Pulpit/ohpc23/atr/Kuba Ostapko/'+'/*[0-2].txt')
for filename in files:
    print('==',filename)
    corr_spctr(filename)
    filename2=os.path.splitext(filename)[0]+'_corr'+os.path.splitext(filename)[1]
    atr(filename2)
    
    
    
    