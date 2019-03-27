# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 13:31:03 2019

@author: Matthew
"""
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib import patches
from scipy.optimize import curve_fit

import os
cdir=os.getcwd()
os.chdir('../..')
from plot_colors import *
os.chdir(cdir)

def blendClr(color1,color2,alphaVal):
    return np.array(color1)*alphaVal + np.array(color2)*(1-alphaVal)


def fit_func(x, a, b, c, d, e):
    return a*(1-b*np.exp(-((x-d)/np.sqrt(2*(c**2)))**2 ) - \
              b*np.exp(-((x-e)/np.sqrt(2*(c**2)))**2 ))

def lin_fit(x, m, b):
    return m*x+b

tiltVals = np.genfromtxt('data/TiltCalHz.csv', delimiter = ',')

EU_num=[]
EU_err=[]

for tv in tiltVals:
    temp=np.genfromtxt('data/'+'EUCal_'+str(tv)+'_Hz.csv', delimiter = ',')
    p=np.mean(temp,1)
    p=p[6::]
    EU_num.append(p)
    nval =  int(np.size(temp,1))
    EU_err.append(np.sqrt(p*(1-p)/nval))
#    EU.append(temp)

EU_num=np.array(EU_num)
EU_err=np.array(EU_err)

paramSave=[]
fitParamSave=[]

for eui in range(len(EU_num[1,::])):
    param, cov = curve_fit(fit_func, tiltVals, EU_num[::,eui], sigma = EU_err[::,eui], p0=[0.9,0.7,0.1,0.75,1.05])
    paramSave.append(param)
    fitParamSave.append(np.sqrt(np.diag(cov)))


pInd = 2
paramSave=np.array(paramSave)
fitParamSave=np.array(fitParamSave)

Es=(paramSave[::,4]+paramSave[::,3])/2
EsSig=np.sqrt(fitParamSave[::,3]**2 + fitParamSave[::,4]**2)/2

Us=(paramSave[::,4]-paramSave[::,3])/2
UsSig=np.sqrt(fitParamSave[::,3]**2 + fitParamSave[::,4]**2)/2


param=paramSave[pInd,::]
lw=2.5
mksz=5
   
FX=6.5*2
FY=FX/4
fig = plt.figure(figsize=(FX,FY))

ax1 = fig.add_subplot(131, aspect='auto')
ax1.errorbar(tiltVals,EU_num[:,pInd],yerr=EU_err[:,pInd], fmt='o',\
             markersize=mksz, linewidth = lw, zorder = 1, label = 'Data')
ax1.plot(tiltVals,EU_num[:,pInd],'o', color='white', markersize=2, zorder = 5 )

tvs=np.linspace(tiltVals[0]-0.1,tiltVals[-1]+0.1,1000)
ax1.plot(tvs,fit_func(tvs,param[0],param[1],param[2],param[3],param[4]),'--', \
         color='gray', zorder = -10, label = 'Fit')

ax1.set_ylim([0,1])
ax1.set_xlim([0.65,1.1])
ax1.set_xlabel('Frequency (kHz)')
ax1.set_ylabel('P(odd)')
ax1.set_title('E,U Calibration')
ax1.grid(True, which="both", ls="-")
#ax1.legend(frameon=False, loc=(1), ncol=2)

ns=len(fitParamSave[::,1])
sites=np.linspace(1,ns,ns)


paramEs, covEs = curve_fit(lin_fit, sites, Es, sigma = EsSig, p0=[-0.01, 0.89])
paramUs, covUs = curve_fit(lin_fit, sites, Us, sigma = UsSig, p0=[0, 0.135])

ax2 = fig.add_subplot(132, aspect='auto')
ax2.errorbar(sites, Es, yerr=EsSig, fmt='o', markersize=mksz, linewidth = lw, zorder = 1)
ax2.plot(sites, Es, 'o', markersize=2, color='white', zorder = 5)
ax2.plot(sites, lin_fit(sites,paramEs[0],paramEs[1]),'--', color='gray', zorder = -5)
ax2.set_ylim([0.87,0.91])
ax2.set_yticks(np.linspace(0.87,0.91,5))
ax2.set_title('On-Site Disorder')
ax2.set_xlabel('sites')
ax2.set_ylabel('Energy (kHz)')
ax2.grid(True, which="both", ls="-")

ax3 = fig.add_subplot(133, aspect='auto')
ax3.errorbar(sites, Us, yerr=UsSig, fmt='o', markersize=mksz, linewidth = lw, zorder = 1)
ax3.plot(sites, Us, 'o', markersize=2, color='white', zorder = 5)
#ax3.plot(sites, lin_fit(sites,paramUs[0],paramUs[1]),'--', zorder = -5)
ax3.set_ylim([0.125,0.145])
ax3.set_yticks(np.linspace(0.125,0.145,5))
ax3.set_title('On-Site U Variation')
ax3.set_xlabel('sites')
ax3.set_ylabel('Energy (kHz)')
ax3.grid(True, which="both", ls="-")

plt.savefig('EUCal.pdf')
