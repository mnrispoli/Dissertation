# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 17:31:04 2019

@author: Matthew
"""

import numpy as np
#from scipy.linalg import expm
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import patches

import os
cdir=os.getcwd()
os.chdir('../..')
from plot_colors import *
os.chdir(cdir)

mpl.rc('font', **fontImp)
mpl.rc('lines', **lineSty)

def blendClr(color1,color2,alphaVal):
    return np.array(color1)*alphaVal + np.array(color2)*(1-alphaVal)

#VoN = np.genfromtxt('data/BSLattDepth.csv', delimiter = ',')
Vo=3
qs = np.genfromtxt('data/BSQs.csv', delimiter = ',')

BS0 = np.genfromtxt('data/BS0_'+str(Vo)+'_Er.csv', delimiter = ',')+Vo/2
BS1 = np.genfromtxt('data/BS1_'+str(Vo)+'_Er.csv', delimiter = ',')+Vo/2
BS2 = np.genfromtxt('data/BS2_'+str(Vo)+'_Er.csv', delimiter = ',')+Vo/2
BS3 = np.genfromtxt('data/BS3_'+str(Vo)+'_Er.csv', delimiter = ',')+Vo/2

bf0 = np.stack(( \
               np.genfromtxt('data/BF0_q0_'+str(Vo)+'_Er.csv', dtype=complex, delimiter = ','), \
               np.genfromtxt('data/BF0_q0p125_'+str(Vo)+'_Er.csv',dtype=complex, delimiter = ','), \
               np.genfromtxt('data/BF0_0p25_'+str(Vo)+'_Er.csv', dtype=complex, delimiter = ','), \
               np.genfromtxt('data/BF0_0p5_'+str(Vo)+'_Er.csv', dtype=complex, delimiter = ',') \
               ))

bf1 = np.stack(( \
               np.genfromtxt('data/BF1_q0_'+str(Vo)+'_Er.csv', dtype=complex, delimiter = ','), \
               np.genfromtxt('data/BF1_q0p125_'+str(Vo)+'_Er.csv',dtype=complex, delimiter = ','), \
               np.genfromtxt('data/BF1_0p25_'+str(Vo)+'_Er.csv', dtype=complex, delimiter = ','), \
               np.genfromtxt('data/BF1_0p5_'+str(Vo)+'_Er.csv', dtype=complex, delimiter = ',') \
               ))
#
#bf2 = np.stack(( \
#               np.genfromtxt('data/BF2_q0_'+str(Vo)+'_Er.csv', dtype=complex, delimiter = ','), \
#               np.genfromtxt('data/BF2_q0p125_'+str(Vo)+'_Er.csv',dtype=complex, delimiter = ','), \
#               np.genfromtxt('data/BF2_0p25_'+str(Vo)+'_Er.csv', dtype=complex, delimiter = ','), \
#               np.genfromtxt('data/BF2_0p5_'+str(Vo)+'_Er.csv', dtype=complex, delimiter = ',') \
#               ))
#
#bf3 = np.stack(( \
#               np.genfromtxt('data/BF3_q0_'+str(Vo)+'_Er.csv', dtype=complex, delimiter = ','), \
#               np.genfromtxt('data/BF3_q0p125_'+str(Vo)+'_Er.csv',dtype=complex, delimiter = ','), \
#               np.genfromtxt('data/BF3_0p25_'+str(Vo)+'_Er.csv', dtype=complex, delimiter = ','), \
#               np.genfromtxt('data/BF3_0p5_'+str(Vo)+'_Er.csv', dtype=complex, delimiter = ',') \
#               ))




NS=4.5
NX=1000
L=35
K=32*4 # pick a multiple of 2^n if you want to show periodicity of n*lattice sites nicely
Er=1.24

qs=np.linspace(-1.0/2,1.0/2-1.0/K,K)

xx=np.linspace(-NS,NS,NX)
xxout=xx.reshape((len(xx),1))

fig = plt.figure(figsize=(6.5*2,2.5*2))

#subplot of excited band
        
# bloch functions of excited band
ax1 = fig.add_subplot(231, aspect='auto')

ax1.plot(xx,4*(1-np.cos(xx*np.pi)**2)-0.5,'-',color = 'black', alpha=0.3, zorder = -10)

Vout=(1-np.cos(xx*np.pi)**2).reshape((len(xx),1))
ax1.plot(xx,np.real(bf1[0,::])/3+3,'-', color = matica99B)
ax1.plot(xx,np.imag(bf1[0,::])/3+3,'--', color = matica99B)

ax1.plot(xx,np.real(bf1[1,::])/3+2,'-', color = matica99C)
ax1.plot(xx,np.imag(bf1[1,::])/3+2,'--', color = matica99C)

ax1.plot(xx,np.real(bf1[2,::])/3+1,'-', color = matica99A)
ax1.plot(xx,np.imag(bf1[2,::])/3+1,'--', color = matica99A)

ax1.plot(xx,np.real(bf1[3,::])/3+0,'-', color = matica99D)
ax1.plot(xx,np.imag(bf1[3,::])/3+0,'--', color = matica99D)

ax1.set_title('Bloch Wavefunctions: n=1',pad=2)
ax1.set_yticks([0,1,2,3])
ax1.set_yticklabels(['q=1/2','q=1/4','q=1/8','q=0'])

ax1.set_xticks([])
ax1.set_xticklabels([])

   

# bloch functions of ground band
ax2 = fig.add_subplot(234, aspect='auto')

ax2.plot(xx,4*(1-np.cos(xx*np.pi)**2)-0.5,'-',color = 'black', alpha=0.3, zorder = -10)

Vout=(1-np.cos(xx*np.pi)**2).reshape((len(xx),1))

ax2.plot(xx,np.real(bf0[0,::])/3,'-', color = matica99B, label = 'Re[f,q]')
ax2.plot(xx,np.imag(bf0[0,::])/3,'--', color = matica99B, label = 'Im[f,q]')

ax2.plot(xx,np.real(bf0[1,::])/3+1,'-', color = matica99C, label = 'Re[f,q]')
ax2.plot(xx,np.imag(bf0[1,::])/3+1,'--', color = matica99C, label = 'Im[f,q]')

ax2.plot(xx,np.real(bf0[2,::])/3+2,'-', color = matica99A, label = 'Re[f,q]')
ax2.plot(xx,np.imag(bf0[2,::])/3+2,'--', color = matica99A, label = 'Im[f,q]')

ax2.plot(xx,np.real(bf0[3,::])/3+3,'-', color = matica99D, label = 'Re[f,q]')
ax2.plot(xx,np.imag(bf0[3,::])/3+3,'--', color = matica99D, label = 'Im[f,q]')

ax2.set_xlabel('Lattice Sites',labelpad=0.5)
ax2.set_title('Bloch Wavefunctions: n=0',pad=2)
ax2.set_yticks([0,1,2,3])
ax2.set_yticklabels(['q=0','q=1/8','q=1/4','q=1/2'])
ax2.legend(frameon=False, loc=(-0.1,-0.45), ncol=4)


# density bloch functions of ground band
ax4 = fig.add_subplot(236, aspect='auto')

ax4.plot(xx,4*(1-np.cos(xx*np.pi)**2)-0.5,'-',color = 'black', alpha=0.3, zorder = -10)

ax4.plot(xx,(abs(bf0[0,::])**2)/3,'-', color = matica99B, label = '|f,q|^2')
#ax4.plot(xx,np.imag(bf0[0,::])/3,'--', color = matica99B, label = 'Im[f,q]')

ax4.plot(xx,(abs(bf0[1,::])**2)/3+1,'-', color = matica99C, label = '|f,q|^2')
#ax4.plot(xx,np.imag(bf0[1,::])/3+1,'--', color = matica99C, label = 'Im[f,q]')

ax4.plot(xx,(abs(bf0[2,::])**2)/3+2,'-', color = matica99A, label = '|f,q|^2')
#ax4.plot(xx,np.imag(bf0[2,::])/3+2,'--', color = matica99A, label = 'Im[f,q]')

ax4.plot(xx,(abs(bf0[3,::])**2)/3+3,'-', color = matica99D, label = '|f,q|^2')
#ax4.plot(xx,np.imag(bf0[3,::])/3+3,'--', color = matica99D, label = 'Im[f,q]')

ax4.set_xlabel('Lattice Sites',labelpad=0.5)
ax4.set_title('Bloch Wavefunctions: n=0', pad=2)
ax4.set_yticks([0,1,2,3])
ax4.set_yticklabels(['q=0','q=1/8','q=1/4','q=1/2'])
ax4.legend(frameon=False, loc=(-0.05,-0.45), ncol=4)
ax4.yaxis.tick_right()

# density bloch functions of excited band
ax5 = fig.add_subplot(233, aspect='auto')

ax5.plot(xx,4*(1-np.cos(xx*np.pi)**2)-0.5,'-',color = 'black', alpha=0.3, zorder = -10)

ax5.plot(xx,(abs(bf1[0,::])**2)/3+3,'-', color = matica99B)
#ax4.plot(xx,np.imag(bf0[0,::])/3,'--', color = matica99B, label = 'Im[f,q]')

ax5.plot(xx,(abs(bf1[1,::])**2)/3+2,'-', color = matica99C)
#ax4.plot(xx,np.imag(bf0[1,::])/3+1,'--', color = matica99C, label = 'Im[f,q]')

ax5.plot(xx,(abs(bf1[2,::])**2)/3+1,'-', color = matica99A)
#ax4.plot(xx,np.imag(bf0[2,::])/3+2,'--', color = matica99A, label = 'Im[f,q]')

ax5.plot(xx,(abs(bf1[3,::])**2)/3+0,'-', color = matica99D)
#ax4.plot(xx,np.imag(bf0[3,::])/3+3,'--', color = matica99D, label = 'Im[f,q]')

#ax5.set_xlabel('Lattice Sites')
ax5.set_title('Bloch Wavefunctions: n=1',pad=2)
ax5.set_yticks([0,1,2,3])
ax5.set_yticklabels(['q=1/2', 'q=1/4', 'q=1/8', 'q=0'])
ax5.yaxis.tick_right()
ax5.legend(frameon=False, loc=(-0.05,-0.35), ncol=4)

ax5.set_xticks([])
ax5.set_xticklabels([])

# subplot of bands
ax3 = fig.add_subplot(132, aspect='auto')

#ax3.plot([qs[K/2],qs[K/2]],[-12,20], '--', \
#           color = matica99B, alpha=0.5, linewidth=2, zorder = -10)
#ax3.plot([qs[K/2+K/8+1],qs[K/2+K/8+1]],[-12,20], '-', \
#           color = matica99C, alpha=0.5, linewidth=2, zorder = -10)
#ax3.plot([qs[K/2+K/4+1],qs[K/2+K/4+1]],[-12,20], '-.', \
#           color = matica99A, alpha=0.5, linewidth=2, zorder = -10)
#ax3.plot([0.5-0.005,0.5-0.005],[-12,20], ':', \
#           color = matica99D, alpha=0.5, linewidth=2, zorder = -10)

#ax3.plot(qs,BS3, color = matica97J, linewidth=3, label = 'n=3')
#ax3.plot(qs,BS2, color = matica97B, linewidth=3, label = 'n=2')

lw=1.5

ax3.plot(qs,BS1, color = band1Clr, linewidth=lw, label = 'n=1', alpha=1)
ax3.plot(qs,BS0, color = band0Clr, linewidth=lw, label = 'n=0', alpha=1)
ax3.plot(qs,np.ones(BS0.shape)*Vo, color = 'gray',ls='--', linewidth=lw, label = 'n=0', alpha=1)


ax3.legend(frameon=False,loc=(1))

[XL,XH,YL,YH] = [-0.5, 0.5, 0.0, 6.0]
[XW,YT]=[XH-XL,YH-YL]
[W,H]=[XH-XL,YH-YL]
AR=H/W/1.618


ax3.set_xlim((XL,XH))
ax3.set_xticks(np.linspace(XL,XH,5))
ax3.set_xticklabels(np.linspace(XL,XH,5))
ax3.set_ylim([YL,YH])
ax3.set_yticks(np.linspace(YL,YH,5))
ax3.set_yticklabels(np.linspace(YL,YH,5))

NB=len(qs)
rr=0.03
rlw=1
circle1 = patches.Ellipse(xy = (0, BS0[NB/2]), width = rr, height = rr*AR, \
                          angle = 0, color = matica99B, fill=False, zorder = 10, \
                          linewidth = rlw, linestyle = '--')

circle2 = patches.Ellipse(xy = (0.125, BS0[NB/2+NB/8]), width = rr, height = rr*AR, \
                          angle = 0, color = matica99C, fill=False, zorder = 10, \
                          linewidth = rlw, linestyle = '--')

circle3 = patches.Ellipse(xy = (0.25, BS0[NB/2+NB/4]), width = rr, height = rr*AR, \
                          angle = 0, color = matica99A, fill=False, zorder = 10, \
                          linewidth = rlw, linestyle = '--')

circle4 = patches.Ellipse(xy = (0.5, BS0[NB-1]),  width = rr, height = rr*AR, \
                          angle = 0, color = matica99D, fill=False, zorder = 10, \
                          linewidth = rlw, linestyle = '--')

ax3.add_artist(circle1)
ax3.add_artist(circle2)
ax3.add_artist(circle3)
ax3.add_artist(circle4)

circle1 = patches.Ellipse(xy = (0, BS1[NB/2]), width = rr, height = rr*AR, \
                          angle = 0, color = matica99B, fill=False, zorder = 10, \
                          linewidth = rlw, linestyle = '--')

circle2 = patches.Ellipse(xy = (0.125, BS1[NB/2+NB/8]), width = rr, height = rr*AR, \
                          angle = 0, color = matica99C, fill=False, zorder = 10, \
                          linewidth = rlw, linestyle = '--')

circle3 = patches.Ellipse(xy = (0.25, BS1[NB/2+NB/4]), width = rr, height = rr*AR, \
                          angle = 0, color = matica99A, fill=False, zorder = 10, \
                          linewidth = rlw, linestyle = '--')

circle4 = patches.Ellipse(xy = (0.5, BS1[NB-1]),  width = rr, height = rr*AR, \
                          angle = 0, color = matica99D, fill=False, zorder = 10, \
                          linewidth = rlw, linestyle = '--')

ax3.add_artist(circle1)
ax3.add_artist(circle2)
ax3.add_artist(circle3)
ax3.add_artist(circle4)

ax3.set_ylabel('E_q^n (Er)', labelpad=-1)
ax3.set_xlabel('q', labelpad=0.5)
ax3.set_title('Lattice Depth: '+str(int(Vo))+' Er',pad=2)

plt.savefig('BlochFunctionsBS2_v2.pdf')
