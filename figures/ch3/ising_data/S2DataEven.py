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

S2N4Thry = np.genfromtxt('4SiteData/S2_Thry_N4.csv', delimiter = ',')
S2N4Data = np.genfromtxt('4SiteData/S2Data_N4.csv', delimiter = ',')

S2N5Thry = np.genfromtxt('5SiteData/S2_Thry_N5.csv', delimiter = ',')
S2N5Data = np.genfromtxt('5SiteData/S2Data_N5.csv', delimiter = ',')

S2N6Thry = np.genfromtxt('6SiteData/S2_Thry_N6.csv', delimiter = ',')
S2N6Data = np.genfromtxt('6SiteData/S2Data_N6.csv', delimiter = ',')

S2N7Thry = np.genfromtxt('7SiteData/S2_Thry_N7.csv', delimiter = ',')
S2N7Data = np.genfromtxt('7SiteData/S2Data_N7.csv', delimiter = ',')

S2N8Thry = np.genfromtxt('8SiteData/S2_Thry_N8.csv', delimiter = ',')
S2N8Data = np.genfromtxt('8SiteData/S2Data_N8.csv', delimiter = ',')
#QWThryNN = np.genfromtxt('data/QSThry_Jall_BlochOsc.csv', delimiter = ',')

s2ind=2
s2thind=270
S2PkVal=[S2N4Data[s2ind,1],S2N5Data[s2ind,1],S2N6Data[s2ind,1],S2N7Data[s2ind,1],S2N8Data[s2ind,1]]
S2PkThry=[S2N4Thry[s2thind,1],S2N5Thry[s2thind,1],S2N6Thry[s2thind,1],S2N7Thry[s2thind,1],S2N8Thry[s2thind,1]]
S2PkMerr=[S2N4Data[s2ind,2],S2N5Data[s2ind,2],S2N6Data[s2ind,2],S2N7Data[s2ind,2],S2N8Data[s2ind,2]]
S2PkPerr=[S2N4Data[s2ind,3],S2N5Data[s2ind,3],S2N6Data[s2ind,3],S2N7Data[s2ind,3],S2N8Data[s2ind,3]]

Er=1.24

lw=2.5
mksz=5
   
FX=6.5*2/3*3
FY=FX/2
fig = plt.figure(figsize=(FX,FY))

ax1 = fig.add_subplot(121, aspect='auto')

# 0th order of kapitza dirac
ax1.errorbar(S2N4Data[::,0],S2N4Data[::,1],yerr=S2N4Data[:,2], fmt='o',\
             markersize=mksz, linewidth = lw, zorder = 40, label = 'Data', color = matica97A)
ax1.plot(S2N4Data[::,0],S2N4Data[::,1],'o', color='white', markersize=2, zorder = 41 )

ax1.plot(S2N4Thry[::,0],S2N4Thry[::,1], zorder = 4,ls='-', color = matica97A)


#ax1.plot(S2N5Thry[::,0],S2N5Thry[::,1], zorder = 5,ls='-')

ax1.errorbar(S2N6Data[::,0],S2N6Data[::,1],yerr=S2N6Data[:,2], fmt='o',\
             markersize=mksz, linewidth = lw, zorder = 60, label = 'Data', color = matica97B)
ax1.plot(S2N6Data[::,0],S2N6Data[::,1],'o', color='white', markersize=2, zorder = 61 )

ax1.plot(S2N6Thry[::,0],S2N6Thry[::,1], zorder = 6, ls='-', color = matica97B)


#ax1.plot(S2N7Thry[::,0],S2N7Thry[::,1], zorder = 7,ls='-')

ax1.errorbar(S2N8Data[::,0],S2N8Data[::,1],yerr=S2N8Data[:,2], fmt='o',\
             markersize=mksz, linewidth = lw, zorder = 80, label = 'Data', color = matica97C)
ax1.plot(S2N8Data[::,0],S2N8Data[::,1],'o', color='white', markersize=2, zorder = 81 )
ax1.plot(S2N8Thry[::,0],S2N8Thry[::,1], zorder = 8, ls='-', color = matica97C)
#ax1.plot(QWThryNN[::,0],QWThryNN[::,2]*1000,'-')
#ax1.plot(QWThryNN[::,0],QWThryNN[::,3]*1000,'-')

ax1.set_yscale('linear')
ax1.set_xlim([-10,255])
ax1.set_ylim([-0.1,1.1])
#ax1.set_xtick(np.linspace(0,0.022,0.004))
#ax1.set_xticklabels(np.linspace(0,0.022,0.004))
ax1.set_xlabel('time (ms)')
ax1.set_ylabel('S_2(t)')
ax1.set_title('Single-Site Renyi Entropy: L_even')
ax1.grid(True, which="both", ls="-")


#
ax2 = fig.add_subplot(122, aspect='auto')
ls=[4,5,6,7,8]
ax2.errorbar(ls,S2PkVal[::],yerr=S2PkPerr[::], fmt='o',\
             markersize=mksz, linewidth = lw, zorder = 80, label = 'Data', color = matica97C)
ax2.plot(ls,S2PkVal[::],'o', color='white', markersize=2, zorder = 81 )
ax2.plot(ls,S2PkThry[::])
ax2.set_xlim([3,9])
ax2.set_ylim([0,1])
#
#ax2.errorbar(S2N5Data[::,0],S2N5Data[::,1],yerr=S2N5Data[:,2], fmt='o',\
#             markersize=mksz, linewidth = lw, zorder = 50, label = 'Data', color = matica97E)
#ax2.plot(S2N5Data[::,0],S2N5Data[::,1],'o', color='white', markersize=2, zorder = 51 )
#ax2.plot(S2N5Thry[::,0],S2N5Thry[::,1], zorder = 5, ls='-', color = matica97E)
#
#ax2.errorbar(S2N7Data[::,0],S2N7Data[::,1],yerr=S2N7Data[:,2], fmt='o',\
#             markersize=mksz, linewidth = lw, zorder = 70, label = 'Data', color = matica97F)
#ax2.plot(S2N7Data[::,0],S2N7Data[::,1],'o', color='white', markersize=2, zorder = 71 )
#ax2.plot(S2N7Thry[::,0],S2N7Thry[::,1], zorder = 7, ls='-', color = matica97F)
#
#ax1.set_yscale('linear')
#ax1.set_xlim([-10,255])
#ax1.set_ylim([-0.1,1.1])
##ax1.set_xtick(np.linspace(0,0.022,0.004))
##ax1.set_xticklabels(np.linspace(0,0.022,0.004))
#ax1.set_xlabel('time (ms)')
#ax1.set_ylabel('S_2(t)')
#ax1.set_title('Single-Site Renyi Entropy: L_even')
#ax1.grid(True, which="both", ls="-")

#ax1.legend(frameon=False, loc=(1), ncol=2)


plt.savefig('S2_Even.pdf')
