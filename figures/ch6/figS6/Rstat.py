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
import matplotlib as mpl

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




L5Rn = np.genfromtxt('RStat_L5.csv', delimiter = ',')
L6Rn = np.genfromtxt('RStat_L6.csv', delimiter = ',')

L7Rn = np.genfromtxt('RStat_L7.csv', delimiter = ',')
L8Rn = np.genfromtxt('RStat_L8.csv', delimiter = ',')


lw=2.5
lws=2
mksz=6
   
FX=6.5*2/3*3
FY=FX/4
fig = plt.figure(figsize=(FX,FY))
ax1 = fig.add_subplot(121, aspect='auto')

#ax1.errorbar(ls,np.array(NeelPkVal[::])*4,yerr=np.array(NeelPkPerr[::])*4, fmt='o',\
#             markersize=mksz, linewidth = lw, zorder = 80, label = 'Data', \
#             color = matica99H)
#ax1.plot(ls,np.array(NeelPkVal[::])*4,'o', color='white', markersize=2, zorder = 81 )
#ax1.plot(L5Rn[::,0],L5Rn[::,1], linewidth=lw, label = 'U=2.7J\n  L=8',\
#         linestyle = '-', color = blendClr(matica99D,clrWhite,0.6))
ax1.plot(L6Rn[::,0]/2,L6Rn[::,1], linewidth=lw, label = 'L=6',\
         linestyle = '-', color = blendClr(matica99A,clrWhite,0.75))
ax1.plot(L7Rn[::,0]/2,L7Rn[::,1], linewidth=lw, label = 'L=7',\
         linestyle = '-', color = blendClr(matica99I,clrWhite,0.75))
ax1.plot(L8Rn[::,0]/2,L8Rn[::,1], linewidth=lw, label = 'L=8',\
         linestyle = '-', color = blendClr(matica99B,clrWhite,0.75))
xi=np.linspace(0,12,5)
ax1.plot(xi,np.ones(xi.shape)*0.39,'k--',linewidth=1.5,zorder=-5)
ax1.plot(xi,np.ones(xi.shape)*0.53,'k--',linewidth=1.5,zorder=-5)
ax1.set_ylim([0.38,0.55])
ax1.set_yticks(np.linspace(0.4,0.55,4))
ax1.set_xlim([0.2,11])

ax1.legend(frameon=False, loc=(0.24,-0.25), ncol=3)

plt.savefig('RStat.pdf')
