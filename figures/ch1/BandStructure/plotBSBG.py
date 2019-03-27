#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 12:06:01 2019

@author: matthewrispoli
"""
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib import patches
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

VoN = np.genfromtxt('data/BSLattDepth_ref.csv', delimiter = ',')
qs = np.genfromtxt('data/BSQs_ref.csv', delimiter = ',')

BS0 = np.genfromtxt('data/BS0_ref.csv', delimiter = ',')
BS1 = np.genfromtxt('data/BS1_ref.csv', delimiter = ',')
BS2 = np.genfromtxt('data/BS2_ref.csv', delimiter = ',')
BS3 = np.genfromtxt('data/BS3_ref.csv', delimiter = ',')


VoN = np.genfromtxt('data/BSLattDepth.csv', delimiter = ',')
qs = np.genfromtxt('data/BSQs.csv', delimiter = ',')

BS0 = np.genfromtxt('data/BS0.csv', delimiter = ',')
BS1 = np.genfromtxt('data/BS1.csv', delimiter = ',')
BS2 = np.genfromtxt('data/BS2.csv', delimiter = ',')
BS3 = np.genfromtxt('data/BS3.csv', delimiter = ',')
BS4 = np.genfromtxt('data/BS4.csv', delimiter = ',')
BS5 = np.genfromtxt('data/BS5.csv', delimiter = ',')
BS6 = np.genfromtxt('data/BS6.csv', delimiter = ',')

#difference in energies between bands at same q
dBS01 = BS1-BS0
dBS12 = BS2-BS1
dBS23 = BS3-BS2
dBS34 = BS4-BS3
dBS45 = BS5-BS4
dBS56 = BS6-BS5

# bandwidths for the first 4 bands
BW0 = np.ones(len(VoN))
BW1 = np.ones(len(VoN))
BW2 = np.ones(len(VoN))
BW3 = np.ones(len(VoN))

# bandgaps between the first 4 bands
BG01min = np.ones(len(VoN))
BG12min = np.ones(len(VoN))
BG23min = np.ones(len(VoN))
BG34min = np.ones(len(VoN))
BG45min = np.ones(len(VoN))
BG56min = np.ones(len(VoN))

BG01max = np.ones(len(VoN))
BG12max = np.ones(len(VoN))
BG23max = np.ones(len(VoN))
BG34max = np.ones(len(VoN))
BG45max = np.ones(len(VoN))
BG56max = np.ones(len(VoN))

NB=7

Bmax = np.ones((NB,len(VoN)))
Bmin = np.ones((NB,len(VoN)))

for nn in range(np.size(BS0,0)):
    BW0[nn]=np.amax(BS0[nn,::])-np.amin(BS0[nn,::])
    BW1[nn]=np.amax(BS1[nn,::])-np.amin(BS1[nn,::])
    BW2[nn]=np.amax(BS2[nn,::])-np.amin(BS2[nn,::])
    BW3[nn]=np.amax(BS3[nn,::])-np.amin(BS3[nn,::])
    
    
    
    BG01min[nn]=np.amin(dBS01[nn,::])
    BG12min[nn]=np.amin(dBS12[nn,::])
    BG23min[nn]=np.amin(dBS23[nn,::])
    BG34min[nn]=np.amin(dBS34[nn,::])
    BG45min[nn]=np.amin(dBS45[nn,::])
    BG56min[nn]=np.amin(dBS56[nn,::])
    
    BG01max[nn]=np.amax(dBS01[nn,::])
    BG12max[nn]=np.amax(dBS12[nn,::])
    BG23max[nn]=np.amax(dBS23[nn,::])
    BG34max[nn]=np.amax(dBS34[nn,::])
    BG45max[nn]=np.amax(dBS45[nn,::])
    BG56max[nn]=np.amax(dBS56[nn,::])
    

    Bmin[0,nn]=np.amin(BS0[nn,::])
    Bmin[1,nn]=np.amin(BS1[nn,::])
    Bmin[2,nn]=np.amin(BS2[nn,::])
    Bmin[3,nn]=np.amin(BS3[nn,::])
    Bmin[4,nn]=np.amin(BS4[nn,::])
    Bmin[5,nn]=np.amin(BS5[nn,::])
    Bmin[6,nn]=np.amin(BS6[nn,::])
    
    Bmax[0,nn]=np.amax(BS0[nn,::])
    Bmax[1,nn]=np.amax(BS1[nn,::])
    Bmax[2,nn]=np.amax(BS2[nn,::])
    Bmax[3,nn]=np.amax(BS3[nn,::])
    Bmax[4,nn]=np.amax(BS4[nn,::])
    Bmax[5,nn]=np.amax(BS5[nn,::])
    Bmax[6,nn]=np.amax(BS6[nn,::])

#plt.semilogy(VoN,BW0)
#plt.semilogy(VoN,BW1)
#plt.semilogy(VoN,BW2)
#plt.semilogy(VoN,BW3)
    
#plt.loglog([1,2,3,4],[BW0[0],BW1[0],BW2[0],BW3[0]])
    
dpth1=3
dpth2=47
dpth3=80+19

fig = plt.figure(figsize=(2*1.618*4,1*4))
#fig = plt.figure(figsize=(6.5,3))

ax161 = fig.add_subplot(161, aspect='0.06')

ax161.plot(qs,BS0[dpth1,::], color = matica97A, label = 'n = 0')
ax161.plot(qs,BS1[dpth1,::], color = matica97B, label = 'n = 1')
ax161.plot(qs,BS2[dpth1,::], color = matica97C, label = 'n = 2')
ax161.plot(qs,BS3[dpth1,::], color = matica97D, label = 'n = 3')
ax161.set_ylim((-18,18))
ax161.set_xlim((-0.5,0.5))
ax161.set_yticks(np.linspace(-16,16,5))
ax161.set_ylabel('Energy (Er)')
ax161.set_xlabel('q')
ax161.set_title(str(VoN[dpth1])+' Er')
ax161.legend(frameon=False)

ax162 = fig.add_subplot(162, aspect='0.06')

ax162.plot(qs,BS0[dpth2,::], color = matica97A)
ax162.plot(qs,BS1[dpth2,::], color = matica97B)
ax162.plot(qs,BS2[dpth2,::], color = matica97C)
ax162.plot(qs,BS3[dpth2,::], color = matica97D)
ax162.set_ylim((-18,18))
ax162.set_xlim((-0.5,0.5))
ax162.set_yticks([])
ax162.set_xlabel('q')
ax162.set_title(str(VoN[dpth2])+' Er')

ax163 = fig.add_subplot(163, aspect='0.06')

ax163.plot(qs,BS0[dpth3,::], color = matica97A)
ax163.plot(qs,BS1[dpth3,::], color = matica97B)
ax163.plot(qs,BS2[dpth3,::], color = matica97C)
ax163.plot(qs,BS3[dpth3,::], color = matica97D)
ax163.set_ylim((-18,18))
ax163.set_xlim((-0.5,0.5))
ax163.set_yticks([])
ax163.set_xlabel('q')
ax163.set_title(str(VoN[dpth3])+' Er')


#ax2 = fig.add_subplot(122, aspect='2.3176')

ax2 = fig.add_subplot(122, aspect='auto')  
ax2.plot(VoN,BG01min,'-',color = matica97E, zorder = 3, label = 'n=0->n=1')
ax2.plot(VoN,BG01max,'-',color = matica97E, zorder = 3)
ax2.fill_between(VoN, BG01min, BG01max, \
                 color = blendClr(matica97E,clrWhite,0.35), \
                 zorder = 3)

ax2.plot(VoN,BG12min,'-', color = matica97F, zorder = 2, label = 'n=1->n=2' )
ax2.plot(VoN,BG12max,'-', color = matica97F, zorder = 2 )
ax2.fill_between(VoN, BG12min, BG12max, \
                 color = blendClr(matica97F,clrWhite,0.35), \
                 zorder = 2)

ax2.plot(VoN,BG23min,'-', color = matica97I, zorder = 1, label = 'n=2->n=3' )
ax2.plot(VoN,BG23max,'-', color = matica97I, zorder = 1 )
ax2.fill_between(VoN, BG23min, BG23max, \
                 color = blendClr(matica97I,clrWhite,0.35), \
                 zorder = 1)

#ax2.plot(VoN,BG34min,'-', color = matica97H, zorder = 0, label = 'n=2->n=3' )
#ax2.plot(VoN,BG34max,'-', color = matica97H, zorder = 0 )
#ax2.fill_between(VoN, BG34min, BG34max, \
#                 color = blendClr(matica97H,clrWhite,0.35), \
#                 zorder = 0)
#
#ax2.plot(VoN,BG45min,'-', color = matica97J, zorder = -1, label = 'n=2->n=3' )
#ax2.plot(VoN,BG45max,'-', color = matica97J, zorder = -1 )
#ax2.fill_between(VoN, BG45min, BG45max, \
#                 color = blendClr(matica97J,clrWhite,0.35), \
#                 zorder = -1)
#
#ax2.plot(VoN,BG56min,'-', color = matica97I, zorder = -2, label = 'n=2->n=3' )
#ax2.plot(VoN,BG56max,'-', color = matica97I, zorder = -2)
#ax2.fill_between(VoN, BG56min, BG56max, \
#                 color = blendClr(matica97I,clrWhite,0.35), \
#                 zorder = -2)
#
#ax2.set_xticks(np.linspace(0,100,11))
#ax2.set_yticks(np.linspace(0,20,6))
#ax2.set_xlim((0.25,1000))
#ax2.set_ylim((0.1,100))
#ax2.set_yscale('log')
#ax2.set_xscale('log')

ax2.set_xticks(np.linspace(0,100,11))
ax2.set_yticks(np.linspace(0,20,6))
ax2.set_xlim((0.24,100))
ax2.set_ylim((0,20))

ax2.set_xlabel('Lattice Depth (Er)')
ax2.set_ylabel('Direct Band Gap (Er)')
ax2.legend(frameon=False, loc=(0))
#plt.semilogy(VoN,BW3)
#plt.Axes.set_xticks(np.linspace(0,45,10))

kmid=32
klast=63

style="Simple,tail_width=1,head_width=5,head_length=6"
bg1 = dict(arrowstyle=style, color=matica97E)

a2 = patches.FancyArrowPatch((0,BS0[dpth3,kmid]+1), (0,BS1[dpth3,kmid]), \
                             connectionstyle="arc3,rad=0", \
                             **bg1)
a3 = patches.FancyArrowPatch((0,BS1[dpth3,kmid]-1),(0,BS0[dpth3,kmid]), \
                             connectionstyle="arc3,rad=0", \
                             **bg1)
ax163.add_patch(a2)
ax163.add_patch(a3)

bg2 = dict(arrowstyle=style, color=matica97F)

a2 = patches.FancyArrowPatch((0,BS1[dpth3,kmid]+1), (0,BS2[dpth3,kmid]), \
                             connectionstyle="arc3,rad=0", \
                             **bg2)
a3 = patches.FancyArrowPatch((0,BS2[dpth3,kmid]-1),(0,BS1[dpth3,kmid]), \
                             connectionstyle="arc3,rad=0", \
                             **bg2)
ax163.add_patch(a2)
ax163.add_patch(a3)


bg3 = dict(arrowstyle=style, color=matica97I)

a2 = patches.FancyArrowPatch((0,BS2[dpth3,kmid]+1), (0,BS3[dpth3,kmid]), \
                             connectionstyle="arc3,rad=0", \
                             **bg3)
a3 = patches.FancyArrowPatch((0,BS3[dpth3,kmid]-1),(0,BS2[dpth3,kmid]), \
                             connectionstyle="arc3,rad=0", \
                             **bg3)
ax163.add_patch(a2)
ax163.add_patch(a3)

plt.savefig('BSScale.pdf')
plt.show()

fig2 = plt.figure()
axbw = fig2.add_subplot(121, aspect='auto')  
axbw.plot(VoN,BW0)#/BW0[0])
axbw.plot(VoN,BW1)#/BW1[0])
axbw.plot(VoN,BW2)#/BW2[0])
axbw.plot(VoN,BW3)#/BW3[0])
axbw.set_xlim([0.0, 100])
axbw.set_ylim([0.001, 8])
axbw.set_xticks(np.linspace(0,100,11))
axbw.set_xscale('linear')
axbw.set_yscale('log')

axe = fig2.add_subplot(122, aspect='auto')  
axe.plot(VoN,np.mean(BS0,1))#/BW0[0])
axe.plot(VoN,np.mean(BS1,1))#/BW1[0])
axe.plot(VoN,np.mean(BS2,1))#/BW2[0])
axe.plot(VoN,np.mean(BS3,1))#/BW3[0])
axe.plot(VoN,np.mean(BS4,1))#/BW3[0])
axe.plot(VoN,np.mean(BS5,1))#/BW3[0])
axe.plot(VoN,np.mean(BS6,1))#/BW3[0])

axe.set_xlim([0.0, 1000])
axe.set_ylim([-240, 80])
axe.set_xticks(np.linspace(0,100,11))
axe.set_xscale('linear')
axe.set_yscale('linear')



axbs_alpha=0.7

fig3 = plt.figure()
axbs = fig3.add_subplot(111, aspect='auto')  

offV = -VoN/2

axbs.plot(VoN,Bmax[0,::] - offV, color = matica97A)#/BW0[0])
axbs.plot(VoN,Bmin[0,::] - offV, color = matica97A)#/BW1[0])
axbs.fill_between(VoN, Bmin[0,::] - offV, Bmax[0,::] - offV, \
                 color = blendClr(matica97A,clrWhite,0.35), \
                 alpha = axbs_alpha, \
                 zorder = 0)

axbs.plot(VoN,Bmax[1,::] - offV, color = matica97B)#/BW0[0])
axbs.plot(VoN,Bmin[1,::] - offV, color = matica97B)#/BW1[0])
axbs.fill_between(VoN, Bmin[1,::] - offV, Bmax[1,::] - offV, \
                 color = blendClr(matica97B,clrWhite,0.35), \
                 alpha = axbs_alpha, \
                 zorder = 0)

axbs.plot(VoN,Bmax[2,::] - offV, color = matica97C)#/BW0[0])
axbs.plot(VoN,Bmin[2,::] - offV, color = matica97C)#/BW1[0])
axbs.fill_between(VoN, Bmin[2,::] - offV, Bmax[2,::] - offV, \
                 color = blendClr(matica97C,clrWhite,0.35), \
                 alpha = axbs_alpha, \
                 zorder = 0)

axbs.plot(VoN,Bmax[3,::] - offV, color = matica97D)#/BW0[0])
axbs.plot(VoN,Bmin[3,::] - offV, color = matica97D)#/BW1[0])
axbs.fill_between(VoN, Bmin[3,::] - offV, Bmax[3,::] - offV, \
                 color = blendClr(matica97D,clrWhite,0.35), \
                 alpha = axbs_alpha, \
                 zorder = 0)


axbs.plot(VoN,Bmax[4,::] - offV, color = matica97E)#/BW0[0])
axbs.plot(VoN,Bmin[4,::] - offV, color = matica97E)#/BW1[0])
axbs.fill_between(VoN, Bmin[4,::] - offV, Bmax[4,::] - offV, \
                 color = blendClr(matica97E,clrWhite,0.35), \
                 alpha = axbs_alpha, \
                 zorder = 0)


axbs.plot(VoN,Bmax[5,::] - offV, color = matica97F)#/BW0[0])
axbs.plot(VoN,Bmin[5,::] - offV, color = matica97F)#/BW1[0])
axbs.fill_between(VoN, Bmin[5,::] - offV, Bmax[5,::] - offV, \
                 color = blendClr(matica97F,clrWhite,0.35), \
                 alpha = axbs_alpha, \
                 zorder = 0)


axbs.plot(VoN,Bmax[6,::] - offV, color = matica97G)#/BW0[0])
axbs.plot(VoN,Bmin[6,::] - offV, color = matica97G)#/BW1[0])
axbs.fill_between(VoN, Bmin[6,::] - offV, Bmax[6,::] - offV, \
                 color = blendClr(matica97G,clrWhite,0.35), \
                 alpha = axbs_alpha, \
                 zorder = 0)



axbs.plot(VoN,VoN/2 - offV, ls='--', color = 'gray', zorder = -1)#/BW2[0])
axbs.plot(VoN,VoN/2+VoN/2 - offV, ls='--', color = 'gray', zorder = -1)#/BW2[0])
axbs.plot(VoN,-VoN/2 - offV, ls='--', color = 'gray', zorder = -1)#/BW3[0])
axbs.fill_between(VoN, -VoN/2 - offV, VoN/2 - offV, \
                 color = 'gray', \
                 alpha = 0.05, \
                 zorder = -2)

axbs.set_xlim([0.0, 100])
axbs.set_ylim([-1, 50])
axbs.set_xticks(np.linspace(0,100,11))
axbs.set_xscale('linear')
axbs.set_yscale('linear')
