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



CABErr = np.genfromtxt('CABEnd.csv', delimiter = ',')
Xdec = np.genfromtxt('allXDecohere.csv', delimiter = ',')
Ydec = np.genfromtxt('allYDecohere.csv', delimiter = ',')

lw=2.5
lws=2
mksz=5
   
FX=6.5*2/3*3
FY=FX/4
fig = plt.figure(figsize=(FX,FY))
ax1 = fig.add_subplot(121, aspect='auto')

#ax1.errorbar(ls,np.array(NeelPkVal[::])*4,yerr=np.array(NeelPkPerr[::])*4, fmt='o',\
#             markersize=mksz, linewidth = lw, zorder = 80, label = 'Data', \
#             color = matica99H)
#ax1.plot(ls,np.array(NeelPkVal[::])*4,'o', color='white', markersize=2, zorder = 81 )
ax1.plot(Xdec[::,0],Ydec[::,0], linewidth=lw, label = 'g=2.00',\
         color = blendClr(matica99D,clrWhite,1))
ax1.plot(Xdec[::,0],Ydec[::,2], linewidth=lw, label = 'g=0.50',\
         color = blendClr(matica99A,clrWhite,1))
ax1.plot(Xdec[::,0],Ydec[::,4], linewidth=lw, label = 'g=0.12',\
         color = blendClr(matica99I,clrWhite,1))
ax1.plot(Xdec[::,0],Ydec[::,6], linewidth=lw, label = 'g=0.05',\
         color = blendClr(matica99B,clrWhite,1))
#ax1.plot(Xdec[::,0],Ydec[::,8], linewidth=lw,\
#         color = blendClr(matica99B,clrWhite,0.8))

ax1.legend(frameon=False, loc=(0.12,-0.25), ncol=4)

ax1.set_xscale('log')
ax1.set_xlim([0.1,250])
ax1.set_xticks(np.logspace(-1,2,4), minor=False)
ax1.tick_params(direction='out', length=3, width=0.75, colors='k',
               grid_color='k', grid_alpha=0.5, which='minor')
ax1.tick_params(direction='out', length=6, width=1, colors='k',
               grid_color='k', grid_alpha=0.5, which='major')
#ax1.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
#ax1.set_xticks(np.logspace(-1,2,4), minor=False)
ax1.set_ylim([0,0.22])
ax1.set_yticks(np.linspace(0,0.2,5))


#
ax2 = fig.add_subplot(122, aspect='auto')
ls=[4,5,6,7,8]
ax2.errorbar(CABErr[::,0],CABErr[::,1],yerr=CABErr[::,2], fmt='o',\
             markersize=mksz, linewidth = lw, zorder = 80,\
             color = matica99E)

ax2.errorbar(CABErr[-1,0],CABErr[-1,1],yerr=CABErr[-1,2], fmt='o',\
             markersize=mksz, linewidth = lw, zorder = 80,\
             color = matica99D)

ax2.errorbar(CABErr[-3,0],CABErr[-3,1],yerr=CABErr[-3,2], fmt='o',\
             markersize=mksz, linewidth = lw, zorder = 80, \
             color = matica99A)

ax2.errorbar(CABErr[-5,0],CABErr[-5,1],yerr=CABErr[-5,2], fmt='o',\
             markersize=mksz, linewidth = lw, zorder = 80, \
             color = matica99I)

ax2.errorbar(CABErr[-7,0],CABErr[-7,1],yerr=CABErr[-7,2], fmt='o',\
             markersize=mksz, linewidth = lw, zorder = 80, \
             color = matica99B)

ax2.plot(CABErr[::,0],CABErr[::,1],'o', color='white', markersize=2, zorder = 81 )

ax2.plot(CABErr[::,0],CABErr[::,1], linewidth=lws,\
         color = blendClr(matica99E,clrWhite,0.5))
ax2.legend('')
ax2.set_xscale('log')
#ax2.set_xlim([3.75,8.15])
#ax2.set_xticks([4,5,6,7,8])



ax2.set_xscale('log')
ax2.set_xlim([0.005,2.50])
ax2.set_xticks(np.logspace(-2,0,3), minor=False)
ax2.tick_params(direction='out', length=3, width=0.75, colors='k',
               grid_color='k', grid_alpha=0.5, which='minor')
ax2.tick_params(direction='out', length=6, width=1, colors='k',
               grid_color='k', grid_alpha=0.5, which='major')
#ax1.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
#ax1.set_xticks(np.logspace(-1,2,4), minor=False)
ax2.set_ylim([0,0.22])
ax2.set_yticks(np.linspace(0,0.2,5))



#ax3 = fig.add_subplot(133, aspect='auto')
#ls=[4,5,6,7,8]
#ax3.errorbar(ls,AFMDWPkVal[::],yerr=AFMDWPkPerr[::], fmt='o',\
#             markersize=mksz, linewidth = lw, zorder = 80, label = 'Data', color = matica99B)
#ax3.plot(ls,AFMDWPkVal[::],'o', color='white', markersize=2, zorder = 81 )
#ax3.plot(ls,AFMDWPkThry[::],color = matica99B)
#
#ax3.set_xlim([3,9])
#ax3.set_ylim([0,0.4])
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


plt.savefig('decoherePlot.pdf')
