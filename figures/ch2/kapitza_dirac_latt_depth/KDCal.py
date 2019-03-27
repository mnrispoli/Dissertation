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

from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
image_file0 = get_pkg_data_filename('raw_imgs/scanAA001.fits')
image_file1 = get_pkg_data_filename('raw_imgs/scanAA004.fits')
image_file2 = get_pkg_data_filename('raw_imgs/scanAA007.fits')
image_file3 = get_pkg_data_filename('raw_imgs/scanAA010.fits')
image_data0 = fits.getdata(image_file0, ext=0)
image_data1 = fits.getdata(image_file1, ext=0)
image_data2 = fits.getdata(image_file2, ext=0)
image_data3 = fits.getdata(image_file3, ext=0)

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

kdm0dat = np.genfromtxt('data/KD_m0.csv', delimiter = ',')
kdm1dat = np.genfromtxt('data/KD_m1.csv', delimiter = ',')
kdmn1dat = np.genfromtxt('data/KD_m_n1.csv', delimiter = ',')


paramSave=[]
fitParamSave=[]

fitRan=12
Er=1.24

param0, cov0 = curve_fit(fit_func_m0, kdm0dat[0:fitRan,0], kdm0dat[0:fitRan,1], sigma = kdm0dat[0:fitRan,2], p0=[260*Er,404])
paramSave.append(param0)
fitParamSave.append(np.sqrt(np.diag(cov0)))

param1, cov1 = curve_fit(fit_func_m1, kdm1dat[0:fitRan,0], kdm1dat[0:fitRan,1], sigma = kdm1dat[0:fitRan,2], p0=[260*Er,404])
paramSave.append(param1)
fitParamSave.append(np.sqrt(np.diag(cov1)))

paramn1, covn1 = curve_fit(fit_func_m1, kdm1dat[0:fitRan,0], kdm1dat[0:fitRan,1], sigma = kdm1dat[0:fitRan,2], p0=[260*Er,404])
paramSave.append(paramn1)
fitParamSave.append(np.sqrt(np.diag(covn1)))

paramSave=np.array(paramSave)
fitParamSave=np.array(fitParamSave)

param = wtAvgParam(fitParamSave,paramSave)

lw=2.5
mksz=5
   
FX=6.5*2
FY=FX/3
fig = plt.figure(figsize=(FX,FY))

ax1 = fig.add_subplot(121, aspect='auto')

# 0th order of kapitza dirac
ax1.errorbar(kdm0dat[::,0]*1E3,kdm0dat[::,1],yerr=kdm0dat[:,2], fmt='o',\
             markersize=mksz, linewidth = lw, zorder = 1, label = 'Data')
ax1.plot(kdm0dat[::,0]*1E3,kdm0dat[::,1],'o', color='white', markersize=2, zorder = 5 )

tvs=np.linspace(0,0.03,1000)
ax1.plot(tvs*1E3,fit_func_m0(tvs,param[0],param[1]),'--', \
         color='gray', zorder = -10, label = 'Fit')


# 1st order of kapitza dirac scattering
ax1.errorbar(kdm1dat[::,0]*1E3,kdm1dat[::,1],yerr=kdm1dat[:,2], fmt='o',\
             markersize=mksz, linewidth = lw, zorder = 1, label = 'Data')
ax1.plot(kdm1dat[::,0]*1E3,kdm1dat[::,1],'o', color='white', markersize=2, zorder = 5 )

tvs=np.linspace(0,0.03,1000)
ax1.plot(tvs*1E3,fit_func_m1(tvs,param[0],param[1]),'--', \
         color='gray', zorder = -10, label = 'Fit')


# -1st order of kaptiza dirac scattering
ax1.errorbar(kdmn1dat[::,0]*1E3,kdmn1dat[::,1],yerr=kdmn1dat[:,2], fmt='o',\
             markersize=mksz, linewidth = lw, zorder = 1, label = 'Data')
ax1.plot(kdmn1dat[::,0]*1E3,kdmn1dat[::,1],'o', color='white', markersize=2, zorder = 5 )



ax1.set_ylim([0,500])
ax1.set_xlim([0,21])
ax1.set_xticks(np.linspace(0,20,6))
ax1.set_xticklabels(np.linspace(0,20,6))
ax1.set_xlabel('time (\mu s)')
ax1.set_ylabel('N (atom Number)')
ax1.set_title('Lattice Depth Calibration')
ax1.grid(True, which="both", ls="-")
#ax1.legend(frameon=False, loc=(1), ncol=2)

xinds=[400,600]
yinds=[250,850]
ax2 = fig.add_subplot(185, aspect='equal')
ax2.imshow(np.transpose(image_data0[0,xinds[0]:xinds[-1],yinds[0]:yinds[-1]]),\
           cmap=plt.cm.Greens)
ax2.axis('off')
ax2.set_title('t=1mus')

ax3 = fig.add_subplot(186, aspect='equal')
ax3.imshow(np.transpose(image_data1[0,xinds[0]:xinds[-1],yinds[0]:yinds[-1]]),\
           cmap=plt.cm.Greens)
ax3.axis('off')
ax3.set_title('t=4mus')

ax4 = fig.add_subplot(187, aspect='equal')
ax4.imshow(np.transpose(image_data2[0,xinds[0]:xinds[-1],yinds[0]:yinds[-1]]),\
           cmap=plt.cm.Greens)
ax4.axis('off')
ax4.set_title('t=7mus')

ax5 = fig.add_subplot(188, aspect='equal')
ax5.imshow(np.transpose(image_data3[0,xinds[0]:xinds[-1],yinds[0]:yinds[-1]]),\
           cmap=plt.cm.Greens)
ax5.axis('off')
ax5.set_title('t=10mus')

plt.savefig('KDCal.pdf')
