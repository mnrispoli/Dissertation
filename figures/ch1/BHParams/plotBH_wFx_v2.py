# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 08:39:49 2019

@author: Matthew
"""

import numpy as np
#from scipy.linalg import expm
import matplotlib.pyplot as plt
import matplotlib as mbpl

import os
cdir=os.getcwd()
os.chdir('../..')
from plot_colors import *
os.chdir(cdir)

mpl.rc('font', **fontImp)
mpl.rc('lines', **lineSty)

Er=1.24 #recoil energy in kHz. All energies saved in the /data/ folder are in kHz

VoN = np.genfromtxt('data/LattOut.csv', delimiter = ',')

wn0 = np.genfromtxt('data/W0.csv', dtype=complex, delimiter = ',')
wn1 = np.genfromtxt('data/W1.csv', dtype=complex, delimiter = ',')
wn2 = np.genfromtxt('data/W2.csv', dtype=complex, delimiter = ',')
wn3 = np.genfromtxt('data/W3.csv', dtype=complex, delimiter = ',')

xx=np.real(wn0[::,0])

BHE0 = np.genfromtxt('data/BHE0.csv', delimiter = ',')
BHE1 = np.genfromtxt('data/BHE1.csv', delimiter = ',')
BHE2 = np.genfromtxt('data/BHE2.csv', delimiter = ',')
BHE3 = np.genfromtxt('data/BHE3.csv', delimiter = ',')

BHJ0 = np.genfromtxt('data/BHJ0.csv', delimiter = ',')
BHJ1 = np.genfromtxt('data/BHJ1.csv', delimiter = ',')
BHJ2 = np.genfromtxt('data/BHJ2.csv', delimiter = ',')
BHJ3 = np.genfromtxt('data/BHJ3.csv', delimiter = ',')

BHU0 = np.genfromtxt('data/BHU0.csv', delimiter = ',')
BHU1 = np.genfromtxt('data/BHU1.csv', delimiter = ',')
BHU2 = np.genfromtxt('data/BHU2.csv', delimiter = ',')
BHU3 = np.genfromtxt('data/BHU3.csv', delimiter = ',')



NX=5001 # number of discrete points to evaluate 
Vp=-Vo*np.cos(xx)/2
Vext=Vp



                
FX=6.5*2
FY=FX/3 #1.618
fig = plt.figure(figsize=(FX, FY))



[XL,XH,YL,YH] = [-2.5,2.5,-4.0,16.0]
[XW,YT]=[XH-XL,YH-YL]

W=XW
H=FY
AR=XW/YT

wn0Inds=[2,4,12,36]
ax1 = fig.add_subplot(121, aspect='auto')
ax1.plot(xx/2/np.pi, Vp/Vo, ls='-', color='gray', lw=2, zorder = 0)
#ax1.plot(xx/2/np.pi,wn0[::,wn0Inds[0]]-0.25, lw=2, zorder = 4)
ax1.plot(xx/2/np.pi,wn0[::,wn0Inds[1]]-0.25, lw=2, zorder = 3, color=band0Clr, alpha=0.4, label=str(VoN[wn0Inds[1]])+' Er' )
ax1.plot(xx/2/np.pi,wn0[::,wn0Inds[2]]-0.25, lw=2, zorder = 2, color=band0Clr, alpha=0.7, label=str(VoN[wn0Inds[2]])+' Er' )
ax1.plot(xx/2/np.pi,wn0[::,wn0Inds[3]]-0.25, lw=2, zorder = 1, color=band0Clr, alpha=1, label=str(VoN[wn0Inds[3]])+' Er' )
ax1.set_xlim([XL, XH])
leg=ax1.legend(loc=(0.24,-0.2), ncol=3)
leg.get_frame().set_linewidth(0.0)
ax1.set_yticks([])
ax1.set_xlabel('x (sites)')
ax1.set_ylabel('w(x)')
ax1.set_title('Ground band Wannier function')


VPind=90

ax2 = fig.add_subplot(122, aspect='auto')
ax2.plot(xx/2/np.pi, 2.5*Vp/Vo+1, ls='-', color='gray', lw=2)
ax2.plot(xx/2/np.pi,wn0[::,VPind]/3.0, lw=2, color = band0Clr)
ax2.plot(xx/2/np.pi,wn1[::,VPind]/3.0+2*1.0/3, lw=2, color = band1Clr)
ax2.plot(xx/2/np.pi,wn2[::,VPind]/3.0+2*2.0/3, lw=2, color = band2Clr)
ax2.plot(xx/2/np.pi,wn3[::,VPind]/3.0+2*3.0/3, lw=2, color = band3Clr)
ax2.set_yticks([0,2.0/3,4.0/3,2])
ax2.set_yticklabels(['n=0','n=1','n=2','n=3'])
ax2.set_xlabel('x (sites)')
ax2.set_title('Wannier function: bands n=0->3')
ax2.set_xlim([XL, XH])
#plt.plot(xx/2/np.pi,xx*0)
plt.savefig('WannierFx_v2.pdf')