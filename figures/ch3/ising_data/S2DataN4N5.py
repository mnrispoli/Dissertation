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

S2N5Thry = np.genfromtxt('5SiteData/S2_Thry_N5.csv', delimiter = ',')
S2N5Data = np.genfromtxt('5SiteData/S2Data_N5.csv', delimiter = ',')

gvals = np.genfromtxt('4SiteData/g_vals_data.csv', delimiter = ',')

NeelN5Data = np.genfromtxt('5SiteData/N5_NeelOrder_Data.csv', delimiter = ',')
NeelN5Thry = np.genfromtxt('5SiteData/N5_NeelOrder_Thry.csv', delimiter = ',')

AFM1Data = np.genfromtxt('5SiteData/AFM1_Data_N5.csv', delimiter = ',')
AFM1Thry = np.genfromtxt('5SiteData/AFM1_Thry_N5.csv', delimiter = ',')

AFM2Data = np.genfromtxt('5SiteData/AFM2_Data_N5.csv', delimiter = ',')
AFM2Thry = np.genfromtxt('5SiteData/AFM2_Thry_N5.csv', delimiter = ',')

AFM3Data = np.genfromtxt('5SiteData/AFM3_Data_N5.csv', delimiter = ',')
AFM3Thry = np.genfromtxt('5SiteData/AFM3_Thry_N5.csv', delimiter = ',')

PMData = np.genfromtxt('5SiteData/PM_Data_N5.csv', delimiter = ',')
PMThry = np.genfromtxt('5SiteData/PM_Thry_N5.csv', delimiter = ',')

#S2N5Thry = np.genfromtxt('5SiteData/S2_Thry_N5.csv', delimiter = ',')
#S2N5Data = np.genfromtxt('5SiteData/S2Data_N5.csv', delimiter = ',')
#
#S2N6Thry = np.genfromtxt('6SiteData/S2_Thry_N6.csv', delimiter = ',')
#S2N6Data = np.genfromtxt('6SiteData/S2Data_N6.csv', delimiter = ',')
#
#S2N7Thry = np.genfromtxt('7SiteData/S2_Thry_N7.csv', delimiter = ',')
#S2N7Data = np.genfromtxt('7SiteData/S2Data_N7.csv', delimiter = ',')
#
#
#S2N8Thry = np.genfromtxt('8SiteData/S2_Thry_N8.csv', delimiter = ',')
#S2N8Data = np.genfromtxt('8SiteData/S2Data_N8.csv', delimiter = ',')
#QWThryNN = np.genfromtxt('data/QSThry_Jall_BlochOsc.csv', delimiter = ',')

Er=1.24

lw=2.5
lws=2
mksz=6
   
FX=6.5*2
FY=FX/3
fig = plt.figure(figsize=(FX,FY))

ax1 = fig.add_subplot(121, aspect='auto')
ax1.set_title('Put Phase Diagram and path here')
hz=np.linspace(0,1.1,100)
hx=-(2.0/4)*(hz**2)+1.0/2
ax1.plot(hx,hz)
ax1.set_xlim([0,0.6])
ax1.set_ylim([0,1.25])
ax1.set_xlabel('h_x')
ax1.set_ylabel('h_z')
# 0th order of kapitza dirac
#ax1.errorbar(S2N4Data[::,0],S2N4Data[::,1],yerr=S2N4Data[:,2], fmt='o',\
#             markersize=mksz, linewidth = lw, zorder = 40, label = 'Data', color = matica97A)
#ax1.plot(S2N4Data[::,0],S2N4Data[::,1],'o', color='white', markersize=2, zorder = 41 )
#
#ax1.plot(S2N4Thry[::,0],S2N4Thry[::,1], zorder = 4,ls='-', color = matica97A)
#
#
##ax1.plot(S2N5Thry[::,0],S2N5Thry[::,1], zorder = 5,ls='-')
#
#ax1.errorbar(S2N6Data[::,0],S2N6Data[::,1],yerr=S2N6Data[:,2], fmt='o',\
#             markersize=mksz, linewidth = lw, zorder = 60, label = 'Data', color = matica97B)
#ax1.plot(S2N6Data[::,0],S2N6Data[::,1],'o', color='white', markersize=2, zorder = 61 )
#
#ax1.plot(S2N6Thry[::,0],S2N6Thry[::,1], zorder = 6, ls='-', color = matica97B)
#
#
##ax1.plot(S2N7Thry[::,0],S2N7Thry[::,1], zorder = 7,ls='-')
#
#ax1.errorbar(S2N8Data[::,0],S2N8Data[::,1],yerr=S2N8Data[:,2], fmt='o',\
#             markersize=mksz, linewidth = lw, zorder = 80, label = 'Data', color = matica97C)
#ax1.plot(S2N8Data[::,0],S2N8Data[::,1],'o', color='white', markersize=2, zorder = 81 )
#ax1.plot(S2N8Thry[::,0],S2N8Thry[::,1], zorder = 8, ls='-', color = matica97C)
##ax1.plot(QWThryNN[::,0],QWThryNN[::,2]*1000,'-')
##ax1.plot(QWThryNN[::,0],QWThryNN[::,3]*1000,'-')

ax2 = fig.add_subplot(222, aspect='auto')
ax2.errorbar(PMData[::,0],PMData[::,1],yerr=PMData[:,2], fmt='o',\
             markersize=mksz, linewidth = lw, zorder = 40, \
             label = 'Data', color = matica99A)
ax2.plot(PMData[::,0],PMData[::,1],'o', color='white', markersize=2, zorder = 41 )

ax2.plot(PMThry[::,0],PMThry[::,1], zorder = 4,ls='-', \
         color = blendClr(matica99A,clrWhite,0.5),linewidth=lws)




ax2.errorbar(AFM1Data[::,0],AFM1Data[::,1],yerr=AFM1Data[:,2], fmt='o',\
             markersize=mksz, linewidth = lw, zorder = 40,\
             label = 'Data', color = matica99B)
ax2.plot(AFM1Data[::,0],AFM1Data[::,1],'o', \
         color='white', markersize=2, zorder = 41 )

ax2.plot(AFM1Thry[::,0],AFM1Thry[::,1], zorder = 4,ls='-', \
         color = blendClr(matica99B,clrWhite,0.5),linewidth=lws)



ax2.errorbar(AFM2Data[::,0],AFM2Data[::,1],yerr=AFM2Data[:,2], fmt='o',\
             markersize=mksz, linewidth = lw, zorder = 40, \
             label = 'Data', color = matica99G)
ax2.plot(AFM2Data[::,0],AFM2Data[::,1],'o', color='white', \
         markersize=2, zorder = 41 )

ax2.plot(AFM2Thry[::,0],AFM2Thry[::,1], zorder = 4,ls='-',\
         color = blendClr(matica99G,clrWhite,0.5),linewidth=lws)



ax2.errorbar(AFM3Data[::,0],AFM3Data[::,1],yerr=AFM3Data[:,2], fmt='o',\
             markersize=mksz, linewidth = lw, zorder = 40, \
             label = 'Data', color = matica99C)
ax2.plot(AFM3Data[::,0],AFM3Data[::,1],'o', color='white', \
         markersize=2, zorder = 41 )

ax2.plot(AFM3Thry[::,0],AFM3Thry[::,1], zorder = 4,ls='-', \
         color = blendClr(matica99C,clrWhite,0.5),linewidth=lws)





ax4 = fig.add_subplot(224, aspect='auto')
ax4.errorbar(S2N5Data[::,0],S2N5Data[::,1],yerr=S2N5Data[:,2], fmt='o',\
             markersize=mksz, linewidth = lw, zorder = 40,\
             label = 'Data', color = matica99D)
ax4.plot(S2N5Data[::,0],S2N5Data[::,1],'o', color='white',\
         markersize=2, zorder = 41 )

ax4.plot(S2N5Thry[::,0],S2N5Thry[::,1], zorder = 4,ls='-',\
         color = blendClr(matica99D,clrWhite,0.5),linewidth=lws)




ax4.errorbar(NeelN5Data[::,0],NeelN5Data[::,1]*4,yerr=NeelN5Data[:,2]*4, fmt='o',\
             markersize=mksz, linewidth = lw, zorder = 40,\
             label = 'Data', color = matica99H)
ax4.plot(NeelN5Data[::,0],NeelN5Data[::,1]*4,'o', color='white', \
         markersize=2, zorder = 41 )

ax4.plot(NeelN5Thry[::,0],NeelN5Thry[::,1]*4, zorder = 4,ls='-',\
         color = blendClr(matica99H,clrWhite,0.5),linewidth=lws)

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

#set axes and stuff here

ax2.set_yscale('linear')
ax2.set_xlim([-10,510])
ax2.set_ylim([-0.1,1.1])
ax2.set_xticks(np.linspace(0,500,9))
xtickLabs=['-13.2','-6.8','-0.5','6.0','12.3',\
           '6.0','-0.6','-6.8','-13.2']
ax2.set_xticklabels(xtickLabs)
#ax2.set_xlabel('g (E-U)/J')
ax2.set_ylabel('Pop(state)')
ax2.set_title('Spin States Populations')
ax2.grid(True, which="both", ls="-")
#ax1.legend(frameon=False, loc=(1), ncol=2)

ax4.set_yscale('linear')
ax4.set_xlim([-10,510])
ax4.set_ylim([-0.1,1.1])
ax4.set_xticks(np.linspace(0,500,9))
xtickLabs=['-13.2','-6.8','-0.5','6.0','12.3',\
           '6.0','-0.6','-6.8','-13.2']
ax4.set_xticklabels(xtickLabs)
ax4.set_xlabel('g (E-U)/J')
ax4.set_ylabel('NeelOrder or S_2')
#ax4.set_title('Spin States Populations')
ax4.grid(True, which="both", ls="-")


plt.savefig('IsingFigN5.pdf')
