#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 13:08:58 2019

@author: matthewrispoli
"""

Wn2D = np.genfromtxt('data/2DScatterWn.csv', delimiter = ',')
HO2D = np.genfromtxt('data/2DScatterHO.csv', delimiter = ',')
DE2D = np.genfromtxt('data/2DScatterDaleyEst.csv', delimiter = ',')

WnAx = np.genfromtxt('data/AxialScatterWn.csv', delimiter = ',')
HOAx = np.genfromtxt('data/AxialScatterHO.csv', delimiter = ',')
DEAx = np.genfromtxt('data/AxialScatterDaleyEst.csv', delimiter = ',')

WnBg = np.genfromtxt('data/BigScatterWn.csv', delimiter = ',')
HOBg = np.genfromtxt('data/BigScatterHO.csv', delimiter = ',')
DEBg = np.genfromtxt('data/BigScatterDaleyEst.csv', delimiter = ',')
    
Data2d = np.genfromtxt('exp_data/Sc2DOut.csv', delimiter = ',')
DataAx = np.genfromtxt('exp_data/ScAxOut.csv', delimiter = ',')
            
FX=6.5*2
FY=FX/3#/1.618

[XL,XH,YL,YH] = [0,0.25,0,4.0]
[XW,YT]=[XH-XL,YH-YL]

W=XW/2
H=FY
AR=XW/YT

fig = plt.figure(figsize=(FX, FY))

ax1 = fig.add_subplot(131, aspect='auto')
ax1.loglog(Wn2D[0,::],Wn2D[1,::],'-', label = 'wn')
ax1.loglog(HO2D[0,::],HO2D[1,::],'-', label = '\psi_HO')
ax1.loglog(DE2D[0,::],DE2D[1,::],'-', label = 'Lamb Dicke Regime')

#ax1.errorbar(Data2d[::,0],Data2d[::,1],yerr=Data2d[:,2], fmt='o',\
#             markersize=mksz, linewidth = lw, zorder = 1, label = 'Data')
#ax1.plot(Data2d[::,0],Data2d[::,1],'o', color='white', markersize=2, zorder = 5 )

ax1.grid(True, which="both", ls="-")
ax1.set_title('Scattering Rate: 2D Lattices')
ax1.set_xlabel('Lattice Depth (Er)')
ax1.set_ylabel('\Gamma_{sc}')
ax1.set_ylim([1E-4,1E-1])
ax1.set_xlim([1E-1,1E3])



ax2 = fig.add_subplot(132, aspect='auto')
ax2.loglog(WnAx[0,::],WnAx[1,::],'-', label = 'wn')
ax2.loglog(HOAx[0,::],HOAx[1,::],'-', label = '\psi_HO')

#ax2.errorbar(DataAx[::,0],DataAx[::,1],yerr=DataAx[:,2], fmt='o',\
#             markersize=mksz, linewidth = lw, zorder = 1, label = 'Data')
#ax2.plot(DataAx[::,0],DataAx[::,1],'o', color='white', markersize=2, zorder = 5 )

#ax2.loglog(DEAx[0,::],DEAx[1,::],'-')
ax2.grid(True, which="both", ls="-")
ax2.set_title('Axial Lattice')
ax2.set_xlabel('Lattice Depth (Er)')
ax2.set_ylabel('\Gamma_{sc}')
ax2.set_ylim([1E-4,1E-1])
ax2.set_xlim([1E-1,1E3])

ax3 = fig.add_subplot(133, aspect='auto')
ax3.loglog(WnBg[0,::],WnBg[1,::],'-', label = 'wn')
ax3.loglog(HOBg[0,::],HOBg[1,::],'-', label = 'psi_HO')
#ax3.loglog(DEBg[0,::],DEBg[1,::],'-')
ax3.grid(True, which="both", ls="-")
ax3.set_title('Big Lattice')
ax3.set_xlabel('Lattice Depth (Er)')
ax3.set_ylabel('\Gamma_{sc}')
ax3.set_ylim([1E-4,1E-1])
ax3.set_xlim([1E1,1E5])

plt.savefig('ScatteringRates.pdf')

Er2D=1.240
ErAx=0.255
ErBg=0.007

latt2D1=8 #ind 48
latt2D2=45 #ind 66
AxOmega=7.9
lattAx=(AxOmega/2/ErAx)**2 #ind 78
BgOmega=2
lattBg=(BgOmega/2/ErBg)**2
