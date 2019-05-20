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
S2N4Global = np.genfromtxt('4SiteData/S2_Global_N4.csv', delimiter = ',')

NeelN4Thry = np.genfromtxt('4SiteData/N4_NeelOrder_Thry.csv', delimiter = ',')
NeelN4Data = np.genfromtxt('4SiteData/N4_NeelOrder_Data.csv', delimiter = ',')

AFMDWN4Thry = np.genfromtxt('4SiteData/N4_AFM_DomainWall_Thry.csv', delimiter = ',')
AFMDWN4Data = np.genfromtxt('4SiteData/N4_AFM_DomainWall_Data.csv', delimiter = ',')


S2N5Thry = np.genfromtxt('5SiteData/S2_Thry_N5.csv', delimiter = ',')
S2N5Data = np.genfromtxt('5SiteData/S2Data_N5.csv', delimiter = ',')
S2N5Global = np.genfromtxt('5SiteData/S2_Global_N5.csv', delimiter = ',')

NeelN5Thry = np.genfromtxt('5SiteData/N5_NeelOrder_Thry.csv', delimiter = ',')
NeelN5Data = np.genfromtxt('5SiteData/N5_NeelOrder_Data.csv', delimiter = ',')

AFMDWN5Thry = np.genfromtxt('5SiteData/N5_AFM_DomainWall_Thry.csv', delimiter = ',')
AFMDWN5Data = np.genfromtxt('5SiteData/N5_AFM_DomainWall_Data.csv', delimiter = ',')


S2N6Thry = np.genfromtxt('6SiteData/S2_Thry_N6.csv', delimiter = ',')
S2N6Data = np.genfromtxt('6SiteData/S2Data_N6.csv', delimiter = ',')
S2N6Global = np.genfromtxt('6SiteData/S2_Global_N6.csv', delimiter = ',')

NeelN6Thry = np.genfromtxt('6SiteData/N6_NeelOrder_Thry.csv', delimiter = ',')
NeelN6Data = np.genfromtxt('6SiteData/N6_NeelOrder_Data.csv', delimiter = ',')

AFMDWN6Thry = np.genfromtxt('5SiteData/N5_AFM_DomainWall_Thry.csv', delimiter = ',')
AFMDWN6Data = np.genfromtxt('5SiteData/N5_AFM_DomainWall_Data.csv', delimiter = ',')


S2N7Thry = np.genfromtxt('7SiteData/S2_Thry_N7.csv', delimiter = ',')
S2N7Data = np.genfromtxt('7SiteData/S2Data_N7.csv', delimiter = ',')
S2N7Global = np.genfromtxt('7SiteData/S2_Global_N7.csv', delimiter = ',')

NeelN7Thry = np.genfromtxt('7SiteData/N7_NeelOrder_Thry.csv', delimiter = ',')
NeelN7Data = np.genfromtxt('7SiteData/N7_NeelOrder_Data.csv', delimiter = ',')

AFMDWN7Thry = np.genfromtxt('7SiteData/N7_AFM_DomainWall_Thry.csv', delimiter = ',')
AFMDWN7Data = np.genfromtxt('7SiteData/N7_AFM_DomainWall_Data.csv', delimiter = ',')


S2N8Thry = np.genfromtxt('8SiteData/S2_Thry_N8.csv', delimiter = ',')
S2N8Data = np.genfromtxt('8SiteData/S2Data_N8.csv', delimiter = ',')
S2N8Global = np.genfromtxt('8SiteData/S2_Global_N8.csv', delimiter = ',')

NeelN8Thry = np.genfromtxt('8SiteData/N8_NeelOrder_Thry.csv', delimiter = ',')
NeelN8Data = np.genfromtxt('8SiteData/N8_NeelOrder_Data.csv', delimiter = ',')

AFMDWN8Thry = np.genfromtxt('8SiteData/N8_AFM_DomainWall_Thry.csv', delimiter = ',')
AFMDWN8Data = np.genfromtxt('8SiteData/N8_AFM_DomainWall_Data.csv', delimiter = ',')

#QWThryNN = np.genfromtxt('data/QSThry_Jall_BlochOsc.csv', delimiter = ',')

s2ind=2
s2thind=270
S2PkVal=[S2N4Data[s2ind,1],S2N5Data[s2ind,1],S2N6Data[s2ind,1],S2N7Data[s2ind,1],S2N8Data[s2ind,1]]
S2PkThry=[S2N4Thry[s2thind,1],S2N5Thry[s2thind,1],S2N6Thry[s2thind,1],S2N7Thry[s2thind,1],S2N8Thry[s2thind,1]]
S2PkMerr=[S2N4Data[s2ind,2],S2N5Data[s2ind,2],S2N6Data[s2ind,2],S2N7Data[s2ind,2],S2N8Data[s2ind,2]]
S2PkPerr=[S2N4Data[s2ind,3],S2N5Data[s2ind,3],S2N6Data[s2ind,3],S2N7Data[s2ind,3],S2N8Data[s2ind,3]]


Neelind=5
Neelthind=499
NeelPkVal=[NeelN4Data[Neelind,1],NeelN5Data[Neelind,1],NeelN6Data[Neelind,1],NeelN7Data[Neelind,1],NeelN8Data[Neelind,1]]
NeelPkThry=[NeelN4Thry[Neelthind,1],NeelN5Thry[Neelthind,1],NeelN6Thry[Neelthind,1],NeelN7Thry[Neelthind,1],NeelN8Thry[Neelthind,1]]
NeelPkMerr=[NeelN4Data[Neelind,2],NeelN5Data[Neelind,2],NeelN6Data[Neelind,2],NeelN7Data[Neelind,2],NeelN8Data[Neelind,2]]
NeelPkPerr=[NeelN4Data[Neelind,3],NeelN5Data[Neelind,3],NeelN6Data[Neelind,3],NeelN7Data[Neelind,3],NeelN8Data[Neelind,3]]

s2ind=-1
S2Global=np.array([S2N4Global[s2ind,1],S2N5Global[s2ind,1],S2N6Global[s2ind,1],S2N7Global[s2ind,1],S2N8Global[s2ind,1]])/1.0


Neelind=5
Neelthind=499
AFMDWPkVal=[AFMDWN4Data[Neelind,1],AFMDWN5Data[Neelind,1],AFMDWN6Data[Neelind,1],AFMDWN7Data[Neelind,1],AFMDWN8Data[Neelind,1]]
AFMDWPkThry=[AFMDWN4Thry[Neelthind,1],AFMDWN5Thry[Neelthind,1],AFMDWN6Thry[Neelthind,1],AFMDWN7Thry[Neelthind,1],AFMDWN8Thry[Neelthind,1]]
AFMDWPkMerr=[AFMDWN4Data[Neelind,2],AFMDWN5Data[Neelind,2],AFMDWN6Data[Neelind,2],AFMDWN7Data[Neelind,2],AFMDWN8Data[Neelind,2]]
AFMDWPkPerr=[AFMDWN4Data[Neelind,3],AFMDWN5Data[Neelind,3],AFMDWN6Data[Neelind,3],AFMDWN7Data[Neelind,3],AFMDWN8Data[Neelind,3]]

Er=1.24

lw=2.5
mksz=5
   
FX=6.5*2/3*3
FY=FX/4
fig = plt.figure(figsize=(FX,FY))
ax1 = fig.add_subplot(121, aspect='auto')
ls=[4,5,6,7,8]
ax1.errorbar(ls,np.array(NeelPkVal[::])*4,yerr=np.array(NeelPkPerr[::])*4, fmt='o',\
             markersize=mksz, linewidth = lw, zorder = 80, label = 'Data', color = matica99A)
ax1.plot(ls,np.array(NeelPkVal[::])*4,'o', color='white', markersize=2, zorder = 81 )
ax1.plot(ls,np.array(NeelPkThry[::])*4, color = matica99A)
ax1.set_xlim([3.75,8.15])
ax1.set_xticks([4,5,6,7,8])
ax1.set_ylim([0,1])
ax1.set_yticks(np.linspace(0,1,5))

#
ax2 = fig.add_subplot(122, aspect='auto')
ls=[4,5,6,7,8]
ax2.errorbar(ls,S2PkVal[::],yerr=S2PkPerr[::], fmt='o',\
             markersize=mksz, linewidth = lw, zorder = 80, label = 'Data', color = matica99B)
ax2.plot(ls,S2PkVal[::],'o', color='white', markersize=2, zorder = 81 )
ax2.plot(ls,S2PkThry[::],color = matica99B)

ax2.plot(ls,S2Global[::],'--',color = 'black')
ax2.plot(ls,S2Global[::]/2.0,'-',color = 'black')
ax2.set_xlim([3.75,8.15])
ax2.set_xticks([4,5,6,7,8])
ax2.set_ylim([0,1.75])
ax2.set_yticks(np.linspace(0,1.75,8))



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


plt.savefig('SingleParam.pdf')
