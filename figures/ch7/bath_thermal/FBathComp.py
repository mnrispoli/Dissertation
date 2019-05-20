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

F8Thry = np.genfromtxt('W9.1ThryF8s.csv', delimiter = ',')
F8Data = np.genfromtxt('W9.1DataF8s.csv', delimiter = ',')
F8Merr = np.genfromtxt('W9.1DataF8sMerr.csv', delimiter = ',')

F12Thry = np.genfromtxt('W9.1ThryF12s.csv', delimiter = ',')
F12Data = np.genfromtxt('W9.1DataF12s.csv', delimiter = ',')
F12Merr = np.genfromtxt('W9.1DataF12sMerr.csv', delimiter = ',')

FDifDat1 = np.genfromtxt('W9.1DifDat1.csv', delimiter = ',')
FDifDat2 = np.genfromtxt('W9.1DifDat2.csv', delimiter = ',')

FDifThry1 = np.genfromtxt('W9.1DifThry1.csv', delimiter = ',')
FDifThry2 = np.genfromtxt('W9.1DifThry2.csv', delimiter = ',')

lw=2.5
mksz=5
   
FX=6.5*2
FY=FX/3/2
fig = plt.figure(figsize=(FX,FY))


ax2 = fig.add_subplot(161, aspect='auto')
ax2.errorbar(F8Data[::,0],F8Data[::,2],yerr=F8Merr[:,2], fmt='o',\
             markersize=mksz, linewidth = lw, zorder = 40, label = 'Data', color = matica99B)

ax2.plot(F8Data[::,0],F8Data[::,2],'o', color='white', markersize=2, zorder = 41 )

ax2.plot(F8Thry[::,0],F8Thry[::,-2], zorder = 4,ls='-', color = matica99B)

ax2.errorbar(F12Data[::,0],F12Data[::,2],yerr=F12Merr[:,2], fmt='o',\
             markersize=mksz, linewidth = lw, zorder = 40, label = 'Data', color = matica99A)
ax2.plot(F12Data[::,0],F12Data[::,2],'o', color='white', markersize=2, zorder = 41 )
ax2.plot(F12Thry[::,0],F12Thry[::,-2], zorder = 4,ls='-', color = matica99A)


ax2.set_yscale('linear')
ax2.set_xscale('log')
ax2.set_xlim([0.01,250])
ax2.set_ylim([-0.05,1.5])


#ax2.set_xlabel('g (E-U)/J')
#ax2.set_ylabel('Pop(state)')
#ax2.set_title('Spin States Populations')
#ax2.grid(True, which="both", ls="-")
#ax1.legend(frameon=False, loc=(1), ncol=2)

ax3 = fig.add_subplot(162, aspect='auto')
ax3.errorbar(F8Data[::,0],F8Data[::,4],yerr=F8Merr[:,4], fmt='o',\
             markersize=mksz, linewidth = lw, zorder = 40, label = 'Data', color = matica99B)

ax3.plot(F8Data[::,0],F8Data[::,4],'o', color='white', markersize=2, zorder = 41 )

ax3.plot(F8Thry[::,0],F8Thry[::,-4], zorder = 4,ls='-', color = matica99B)

ax3.errorbar(F12Data[::,0],F12Data[::,4],yerr=F12Merr[:,4], fmt='o',\
             markersize=mksz, linewidth = lw, zorder = 40, label = 'Data', color = matica99A)
ax3.plot(F12Data[::,0],F12Data[::,4],'o', color='white', markersize=2, zorder = 41 )
ax3.plot(F12Thry[::,0],F12Thry[::,-4], zorder = 4,ls='-', color = matica99A)


ax3.set_yscale('linear')
ax3.set_xscale('log')
ax3.set_xlim([0.01,250])
ax3.set_ylim([-0.05,1.5])

ax4 = fig.add_subplot(163, aspect='auto')
ax4.errorbar(F8Data[::,0],F8Data[::,6],yerr=F8Merr[:,6], fmt='o',\
             markersize=mksz, linewidth = lw, zorder = 40, label = 'Data', color = matica99B)

ax4.plot(F8Data[::,0],F8Data[::,6],'o', color='white', markersize=2, zorder = 41 )

ax4.plot(F8Thry[::,0],F8Thry[::,-6], zorder = 4,ls='-', color = matica99B)


ax4.errorbar(F12Data[::,0],F12Data[::,6],yerr=F12Merr[:,6], fmt='o',\
             markersize=mksz, linewidth = lw, zorder = 40, label = 'Data', color = matica99A)
ax4.plot(F12Data[::,0],F12Data[::,6],'o', color='white', markersize=2, zorder = 41 )
ax4.plot(F12Thry[::,0],F12Thry[::,-6], zorder = 4,ls='-', color = matica99A)


ax4.set_yscale('linear')
ax4.set_xscale('log')
ax4.set_xlim([0.01,250])
ax4.set_ylim([-0.05,1.5])


ax5 = fig.add_subplot(164, aspect='auto')
ax5.errorbar(FDifDat1[::,0],FDifDat1[::,1],yerr=FDifDat1[:,2], fmt='o',\
             markersize=mksz, linewidth = lw, zorder = 40, label = 'Data', color = matica99E)

ax5.plot(FDifDat1[::,0],FDifDat1[::,1],'o', color='white', markersize=2, zorder = 41 )

ax5.plot(FDifDat1[::,0],FDifThry1[::], zorder = 4,ls='-', color = matica99E)


ax5.errorbar(FDifDat2[::,0],FDifDat2[::,1],yerr=FDifDat2[:,2], fmt='o',\
             markersize=mksz, linewidth = lw, zorder = 40, label = 'Data', color = matica99D)

ax5.plot(FDifDat2[::,0],FDifDat2[::,1],'o', color='white', markersize=2, zorder = 41 )
ax5.plot(FDifDat2[::,0],FDifThry2[::], zorder = 4,ls='-', color = matica99D)


ax5.set_yscale('linear')
ax5.set_xscale('linear')
ax5.set_xlim([0.5,6.5])
ax5.set_ylim([-0.1,0.75])


#ax2.set_xlabel('g (E-U)/J')
#ax2.set_ylabel('Pop(state)')
#ax2.set_title('Spin States Populations')
#ax4.grid(True, which="both", ls="-")

plt.savefig('BathSSFig.pdf')
