# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 13:31:03 2019

@author: Matthew
"""
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib import patches
from scipy.optimize import curve_fit
import scipy.special as scp

import os
cdir=os.getcwd()
os.chdir('../..')
from plot_colors import *
os.chdir(cdir)

def blendClr(color1,color2,alphaVal):
    return np.array(color1)*alphaVal + np.array(color2)*(1-alphaVal)


def fit_func_m0(t,v,a):
    return a*(scp.jv(0,(t)*v))**2
def fit_func_m1(t,v,a):
    return a*(scp.jv(1,(t)*v))**2
    
#    return a*(1-b*np.exp(-((x-d)/np.sqrt(2*(c**2)))**2 ) - \
#              b*np.exp(-((x-e)/np.sqrt(2*(c**2)))**2 ))
    
def wtAvgParam(s,f):
    w=1/(s**2)
    return np.divide(np.sum(np.multiply(w,f),0),np.sum(w,0))

def lin_fit(x, m, b):
    return m*x+b

QWDat = np.genfromtxt('data/QWCal_Data.csv', delimiter = ',')
QWThry = np.genfromtxt('data/QWThry.csv', delimiter = ',')
QWThryNN = np.genfromtxt('data/QSThry_Jall.csv', delimiter = ',')
#QWThryNN = np.genfromtxt('data/QSThry_Jall_BlochOsc.csv', delimiter = ',')

Er=1.24

lw=2.5
mksz=6
   
FX=3
FY=FX/1.5
fig = plt.figure(figsize=(FX,FY))

ax1 = fig.add_subplot(111, aspect='auto')

# 0th order of kapitza dirac
ax1.errorbar(QWDat[::,0],QWDat[::,1],yerr=QWDat[:,2]*2, fmt='o',\
             markersize=mksz, linewidth = lw, zorder = 1, label = 'Data',\
             color=matica99B)

ax1.plot(QWDat[::,0],QWDat[::,1],'o', color='white', markersize=2, zorder = 5 )

ax1.plot(QWThryNN[::,0],QWThryNN[::,1]*1000, zorder = -2,ls='--', \
         color=blendClr(matica99B,clrBlack,0.8), linewidth=lw, label='Fit NN')
ax1.plot(QWThryNN[::,0],QWThryNN[::,5]*1000, zorder = -1,ls='-', \
         color=blendClr(matica99B,clrWhite,0.5), linewidth=lw, label='Fit All')
#ax1.plot(QWThryNN[::,0],QWThryNN[::,2]*1000,'-')
#ax1.plot(QWThryNN[::,0],QWThryNN[::,3]*1000,'-')
ax1.set_yscale('linear')
ax1.set_xlim([1.5,6.5])
ax1.set_ylim([50,200])
#ax1.set_xtick(np.linspace(0,0.022,0.004))
#ax1.set_xticklabels(np.linspace(0,0.022,0.004))
ax1.set_xlabel('Lattice Depth (E_r)')
ax1.set_ylabel('J (Hz)')
ax1.set_title('Tunneling Calibration ')
ax1.grid(True, which="both", ls="-")
#ax1.legend(frameon=False, loc=(1), ncol=2)


plt.savefig('QWCalv2.pdf')
