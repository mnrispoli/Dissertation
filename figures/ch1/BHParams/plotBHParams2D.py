# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 08:39:49 2019

@author: Matthew
"""

import numpy as np
#from scipy.linalg import expm
import matplotlib.pyplot as plt
import matplotlib as mpl

import os
cdir=os.getcwd()
os.chdir('../..')
from plot_colors import *
os.chdir(cdir)

mpl.rc('font', **fontImp)
mpl.rc('lines', **lineSty)

def blendClr(color1,color2,alphaVal):
    return np.array(color1)*alphaVal + np.array(color2)*(1-alphaVal)


Er=1.24 #recoil energy in kHz. All energies saved in the /data/ folder are in kHz

VoN = np.genfromtxt('data/LattOut.csv', delimiter = ',')

#wn0 = np.genfromtxt('data/W0.csv', dtype=complex, delimiter = ',')
#wn1 = np.genfromtxt('data/W1.csv', dtype=complex, delimiter = ',')
#wn2 = np.genfromtxt('data/W2.csv', dtype=complex, delimiter = ',')
#wn3 = np.genfromtxt('data/W3.csv', dtype=complex, delimiter = ',')
#
#xx=np.real(wn0[::,0])

BHE0 = np.genfromtxt('data/BHE0.csv', delimiter = ',')
BHE1 = np.genfromtxt('data/BHE1.csv', delimiter = ',')
BHE2 = np.genfromtxt('data/BHE2.csv', delimiter = ',')
BHE3 = np.genfromtxt('data/BHE3.csv', delimiter = ',')

BHJ0 = np.genfromtxt('data/BHJ0.csv', delimiter = ',')
BHJ1 = np.genfromtxt('data/BHJ1.csv', delimiter = ',')
BHJ2 = np.genfromtxt('data/BHJ2.csv', delimiter = ',')
BHJ3 = np.genfromtxt('data/BHJ3.csv', delimiter = ',')

BHU0 = np.genfromtxt('data/BHU0.csv', delimiter = ',')

U_big_Cal=132.0 #(1) error
U_axial_Cal=392.0 #(3) error probably
V_big_Cal=15.0
V_axial_Cal=15.0

U_big_calF=U_big_Cal*(1/BHU0[1,30])
U_axial_calF=U_axial_Cal*(1/BHU0[1,30])

BHU0_big_1d = BHU0[1,::]*U_big_calF
BHU0_axial_1d = BHU0[1,::]*U_axial_calF

BHU0_big = BHU0_big_1d[::]*(BHU0_big_1d[::]/BHU0_big_1d[90])
BHU0_axial = BHU0_axial_1d[::]*(BHU0_axial_1d[::]/BHU0_axial_1d[90])

# we can additionally tune the depth in the z-confining lattice which gives us 
#a range of wave function compression independent the physics lattice

ZLattF = (1.0/6.0)**(1.0/4.0) # this is the approximate scaling you get from deep lattice

BHU1 = np.genfromtxt('data/BHU1.csv', delimiter = ',')

BHU1_big = BHU1[1,::]*U_big_calF
BHU1_axial = BHU1[1,::]*U_axial_calF


BHU2 = np.genfromtxt('data/BHU2.csv', delimiter = ',')

BHU2_big = BHU2[1,::]*U_big_calF
BHU2_axial = BHU2[1,::]*U_axial_calF

BHU3 = np.genfromtxt('data/BHU3.csv', delimiter = ',')

BHU3_big = BHU3[1,::]*U_big_calF
BHU3_axial = BHU3[1,::]*U_axial_calF




#NX=5001 # number of discrete points to evaluate 
#Vp=-Vo*np.cos(xx)/2
#Vext=Vp



                
FX=6.5*2
FY=FX/3#/1.618


[XL,XH,YL,YH] = [-2.0,2.0,-4.0,16.0]
[XW,YT]=[XH-XL,YH-YL]

W=XW/2
H=FY
AR=XW/YT

fig = plt.figure(figsize=(FX, FY))

ax1 = fig.add_subplot(131, aspect='auto')
ax1.semilogy(VoN, BHJ0[1,::]*1000, linewidth=2, color = matica97A)
#ax1.semilogy(VoN,-BHJ0[1,::]/Er)
#ax1.semilogy(VoN,BHJ2[3,::]/Er)
#ax1.semilogy(VoN,-BHJ0[4,::]/Er)
#ax1.semilogy(VoN,BHJ0[5,::]/Er)
ax1.set_xlim([0,30])
ax1.set_xticks([0,15,30])
ax1.set_ylim([0.1,Er*1000/4.0])
ax1.set_ylabel('Tunneling J (Hz)')
ax1.set_xlabel('Lattice Depth (Er)')
ax1.grid(True, which="both", ls="-")

ax2 = fig.add_subplot(132, aspect='auto')
ax2.plot(VoN, BHU0_axial, linewidth = 2, color = blendClr(matica99A,clrWhite,0.6), label = 'Axial Lattice' )
ax2.plot(VoN, BHU0_axial*ZLattF, linewidth = 2, color = blendClr(matica99A,clrWhite,0.6))
ax2.fill_between(VoN, BHU0_axial*ZLattF, BHU0_axial, \
                 color = blendClr(matica99A,clrWhite,0.35), \
                 zorder = 1)

ax2.plot(VoN, BHU0_big, linewidth = 2, color = matica97D ,  label = 'Big Lattice')
ax2.plot(VoN, BHU0_big*ZLattF, linewidth = 2, color = matica97D)
ax2.fill_between(VoN, BHU0_big*ZLattF, BHU0_big, \
                 color = blendClr(matica97D,clrWhite,0.35), \
                 zorder = 1)

ax2.set_xlim([0,30])
ax2.set_xticks([0,15,30])
ax2.set_ylim([0,450])
ax2.set_ylabel('Interaction U (Hz)')
ax2.set_xlabel('Lattice Depth (Er)')
#leg=ax2.legend(loc=(2))
leg=ax2.legend(frameon=False, loc=(0.12,1), ncol=2)
#leg.get_frame().set_linewidth(0.0)
ax2.grid()


ax3 = fig.add_subplot(133, aspect='auto')
ujBigH=BHU0_big/(BHJ0[1,::]*1000)
ujBigL=BHU0_big/(BHJ0[1,::]*1000)*ZLattF
ax3.semilogy(VoN, ujBigH, linewidth=2, color = matica97E, label = 'Axial Lattice')
ax3.semilogy(VoN, ujBigL, linewidth=2, color = matica97E)
ax3.fill_between(VoN, ujBigL, ujBigH, \
                 color = blendClr(matica97E,clrWhite,0.35), \
                 zorder = 1)

ujAxH=BHU0_axial/(BHJ0[1,::]*1000)
ujAxL=BHU0_axial/(BHJ0[1,::]*1000)*ZLattF

ax3.semilogy(VoN, ujAxH, linewidth=2, color = matica97I, label = 'Big Lattice')
ax3.semilogy(VoN, ujAxL, linewidth=2, color = matica97I)

ax3.fill_between(VoN, ujAxL, ujAxH, \
                 color = blendClr(matica97E,clrWhite,0.35), \
                 alpha = 1, \
                 zorder = 1)

ax3.set_xticks([0,15,30])
ax3.set_ylabel('U/J')
ax3.set_xlabel('Lattice Depth (Er)')
ax3.set_xlim([0,30])
ax3.set_ylim([0.1,1E3])
leg=ax3.legend(frameon=False, loc=(0.12,1), ncol=2)
#leg.get_frame().set_linewidth(0.0)
ax3.grid(True, which="both", ls="-", zorder=-1000)


plt.savefig('BHParams2D.pdf')
#ax1.plot(xx/2/np.pi,wn0[::,1]/3.0)
#ax1.plot(xx/2/np.pi,wn1[::,1]/3.0+2*1.0/3)
#ax1.plot(xx/2/np.pi,wn2[::,1]/3.0+2*2.0/3)
#ax1.plot(xx/2/np.pi,wn3[::,1]/3.0+2*3.0/3)
#ax1.set_xlim([-3.5, 3.5])
#plt.plot(xx/2/np.pi,xx*0)
