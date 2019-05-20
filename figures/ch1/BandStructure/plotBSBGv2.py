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
    
dpth1=7
dpth2=39-4
dpth3=80+3

fig = plt.figure(figsize=(2*1.618*4,1*4))
#fig = plt.figure(figsize=(6.5,3))
#
#ax161 = fig.add_subplot(161, aspect='0.06')
#
#ax161.plot(qs,BS0[dpth1,::]+VoN[dpth1]/2, color = band0Clr, label = 'n = 0')
##ax161.fill_between(qs, BS0[dpth1,::]+VoN[dpth1]/2, \
##                 BS1[dpth1,::]+VoN[dpth1]/2, \
##                 color = blendClr(matica97E,clrWhite,0.35), \
##                 zorder = -3)
#ax161.plot(qs,BS1[dpth1,::]+VoN[dpth1]/2, color = band1Clr, label = 'n = 1')
##ax161.fill_between(qs, BS1[dpth1,::]+VoN[dpth1]/2, \
##                 BS2[dpth1,::]+VoN[dpth1]/2, \
##                 color = blendClr(matica97F,clrWhite,0.35), \
##                 zorder = -3)
#ax161.plot(qs,BS2[dpth1,::]+VoN[dpth1]/2, color = band2Clr, label = 'n = 2')
##ax161.fill_between(qs, BS2[dpth1,::]+VoN[dpth1]/2, \
##                 BS3[dpth1,::]+VoN[dpth1]/2, \
##                 color = blendClr(matica97I,clrWhite,0.35), \
##                 zorder = -3)
#ax161.plot(qs,BS3[dpth1,::]+VoN[dpth1]/2, color = band3Clr, label = 'n = 3')
#ax161.plot(qs,BS4[dpth1,::]+VoN[dpth1]/2, color = band4Clr, label = 'n = 4')
#ax161.plot(qs,BS5[dpth1,::]+VoN[dpth1]/2, color = band5Clr, label = 'n = 5')
#ax161.plot(qs,np.ones(qs.shape)*VoN[dpth1],ls='--',color = 'gray')
#ax161.set_ylim((0,30))
#ax161.set_xlim((-0.5,0.5))
#ax161.set_yticks(np.linspace(0,30,5))
#ax161.set_ylabel('Energy (Er)')
#ax161.set_xlabel('q')
#ax161.set_title(str(VoN[dpth1])+' Er')
##ax161.legend(frameon=False)
#
#ax162 = fig.add_subplot(162, aspect='0.06')
#
#ax162.plot(qs,BS0[dpth2,::]+VoN[dpth2]/2, color = band0Clr)
#ax162.plot(qs,BS1[dpth2,::]+VoN[dpth2]/2, color = band1Clr)
#ax162.plot(qs,BS2[dpth2,::]+VoN[dpth2]/2, color = band2Clr)
#ax162.plot(qs,BS3[dpth2,::]+VoN[dpth2]/2, color = band3Clr)
#ax162.plot(qs,BS4[dpth2,::]+VoN[dpth2]/2, color = band4Clr)
#ax162.plot(qs,np.ones(qs.shape)*VoN[dpth2],ls='--',color = 'gray')
#ax162.set_ylim((0,30))
#ax162.set_xlim((-0.5,0.5))
#ax162.set_yticks([])
#ax162.set_xlabel('q')
#ax162.set_title(str(VoN[dpth2])+' Er')
#
#ax163 = fig.add_subplot(163, aspect='0.06')
#
#ax163.plot(qs,BS0[dpth3,::]+VoN[dpth3]/2, color = band0Clr)
#ax163.plot(qs,BS1[dpth3,::]+VoN[dpth3]/2, color = band1Clr)
#ax163.plot(qs,BS2[dpth3,::]+VoN[dpth3]/2, color = band2Clr)
#ax163.plot(qs,BS3[dpth3,::]+VoN[dpth3]/2, color = band3Clr)
#ax163.plot(qs,np.ones(qs.shape)*VoN[dpth3],ls='--',color = 'gray')
#ax163.set_ylim((0,30))
#ax163.set_xlim((-0.5,0.5))
#ax163.set_yticks([])
#ax163.set_xlabel('q')
#ax163.set_title(str(VoN[dpth3])+' Er')


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

ax2.plot(VoN,BG34min,'-', color = matica97G, zorder = 0, label = 'n=3->n=4' )
ax2.plot(VoN,BG34max,'-', color = matica97G, zorder = 0 )
ax2.fill_between(VoN, BG34min, BG34max, \
                 color = blendClr(matica97G,clrWhite,0.35), \
                 zorder = 0)

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

ax2.set_xticks(np.linspace(0,50,11))
ax2.set_yticks(np.linspace(0,16,5))
ax2.set_xlim((0.24,50))
ax2.set_ylim((0,16))

ax2.set_xlabel('Lattice Depth (Er)')
ax2.set_ylabel('Direct Band Gap (Er)')
#ax2.legend(frameon=False, loc=(0))
ax2.legend(frameon=False, loc=(0.025,-0.275), ncol=4)
#plt.semilogy(VoN,BW3)
#plt.Axes.set_xticks(np.linspace(0,45,10))




#fig2 = plt.figure()
#axbw = fig2.add_subplot(121, aspect='auto')  
#axbw.plot(VoN,BW0)#/BW0[0])
#axbw.plot(VoN,BW1)#/BW1[0])
#axbw.plot(VoN,BW2)#/BW2[0])
#axbw.plot(VoN,BW3)#/BW3[0])
#axbw.set_xlim([0.0, 100])
#axbw.set_ylim([0.001, 8])
#axbw.set_xticks(np.linspace(0,100,11))
#axbw.set_xscale('linear')
#axbw.set_yscale('log')
#
#axe = fig2.add_subplot(122, aspect='auto')  
#axe.plot(VoN,np.mean(BS0,1))#/BW0[0])
#axe.plot(VoN,np.mean(BS1,1))#/BW1[0])
#axe.plot(VoN,np.mean(BS2,1))#/BW2[0])
#axe.plot(VoN,np.mean(BS3,1))#/BW3[0])
#axe.plot(VoN,np.mean(BS4,1))#/BW3[0])
#axe.plot(VoN,np.mean(BS5,1))#/BW3[0])
#axe.plot(VoN,np.mean(BS6,1))#/BW3[0])
#
#axe.set_xlim([0.0, 1000])
#axe.set_ylim([-240, 80])
#axe.set_xticks(np.linspace(0,100,11))
#axe.set_xscale('linear')
#axe.set_yscale('linear')



axbs_alpha=0.7

#fig3 = plt.figure()
axbs = fig.add_subplot(121, aspect='auto')  

offV = -VoN/2
#matica99G,C
axbs.plot(VoN,Bmax[0,::] - offV, color = band0Clr, label = 'n=0')#/BW0[0])
axbs.plot(VoN,Bmin[0,::] - offV, color = band0Clr)#/BW1[0])
axbs.fill_between(VoN, Bmin[0,::] - offV, Bmax[0,::] - offV, \
                 color = blendClr(band0Clr,clrWhite,0.35), \
                 alpha = axbs_alpha, \
                 zorder = 0)
 #matica99D,F
axbs.plot(VoN,Bmax[1,::] - offV, color = band1Clr, label = 'n=1')#/BW0[0])
axbs.plot(VoN,Bmin[1,::] - offV, color = band1Clr)#/BW1[0])
axbs.fill_between(VoN, Bmin[1,::] - offV, Bmax[1,::] - offV, \
                 color = blendClr(band1Clr,clrWhite,0.35), \
                 alpha = axbs_alpha, \
                 zorder = 0)


axbs.plot(VoN,Bmax[2,::] - offV, color = band2Clr, label = 'n=2')#/BW0[0])
axbs.plot(VoN,Bmin[2,::] - offV, color = band2Clr)#/BW1[0])
axbs.fill_between(VoN, Bmin[2,::] - offV, Bmax[2,::] - offV, \
                 color = blendClr(band2Clr,clrWhite,0.35), \
                 alpha = axbs_alpha, \
                 zorder = 0)


axbs.plot(VoN,Bmax[3,::] - offV, color = band3Clr, label = 'n=3')#/BW0[0])
axbs.plot(VoN,Bmin[3,::] - offV, color = band3Clr)#/BW1[0])
axbs.fill_between(VoN, Bmin[3,::] - offV, Bmax[3,::] - offV, \
                 color = blendClr(band3Clr,clrWhite,0.35), \
                 alpha = axbs_alpha, \
                 zorder = 0)


axbs.plot(VoN,Bmax[4,::] - offV, color = band4Clr, label = 'n=4')#/BW0[0])
axbs.plot(VoN,Bmin[4,::] - offV, color = band4Clr)#/BW1[0])
axbs.fill_between(VoN, Bmin[4,::] - offV, Bmax[4,::] - offV, \
                 color = blendClr(band4Clr,clrWhite,0.35), \
                 alpha = axbs_alpha, \
                 zorder = 0)


axbs.plot(VoN,Bmax[5,::] - offV, color = band5Clr, label = 'n=5')#/BW0[0])
axbs.plot(VoN,Bmin[5,::] - offV, color = band5Clr)#/BW1[0])
axbs.fill_between(VoN, Bmin[5,::] - offV, Bmax[5,::] - offV, \
                 color = blendClr(band5Clr,clrWhite,0.35), \
                 alpha = axbs_alpha, \
                 zorder = 0)



axbs.plot(VoN,Bmax[6,::] - offV, color = band6Clr, label = 'n=6')#/BW0[0])
axbs.plot(VoN,Bmin[6,::] - offV, color = band6Clr)#/BW1[0])
axbs.fill_between(VoN, Bmin[6,::] - offV, Bmax[6,::] - offV, \
                 color = blendClr(band6Clr,clrWhite,0.35), \
                 alpha = axbs_alpha, \
                 zorder = 0)



axbs.plot(VoN,VoN/2 - offV, ls='--', color = 'gray', zorder = -1, label = 'Lattice')#/BW2[0])
#axbs.plot(VoN,VoN/2+VoN/2 - offV, ls='--', color = 'gray', zorder = -1)#/BW2[0])
axbs.plot(VoN,-VoN/2 - offV, ls='--', color = 'gray', zorder = -1)#/BW3[0])
axbs.fill_between(VoN, -VoN/2 - offV, VoN/2 - offV, \
                 color = 'gray', \
                 alpha = 0.1, \
                 zorder = -2)

axbs.set_xlim([0.1, 50])
axbs.set_ylim([-1, 45])
axbs.set_xticks(np.linspace(0,50,11))
axbs.set_xscale('linear')
axbs.set_yscale('linear')
axbs.set_ylabel('Energy (Er)')
axbs.set_xlabel('Lattice Depth (Er)')
#legend(frameon=False, loc=(0))
axbs.legend(frameon=False, loc=(0.175,-0.3), ncol=4)




kmid=32
klast=63
dpthInd=119+4*8

style="Simple,tail_width=1,head_width=5,head_length=6"
bg1 = dict(arrowstyle=style, color=matica97E)

a2 = patches.FancyArrowPatch((VoN[dpthInd],BS0[dpthInd,kmid]+1+VoN[dpthInd]/2.0), (VoN[dpthInd],BS1[dpthInd,kmid]+VoN[dpthInd]/2.0), \
                             connectionstyle="arc3,rad=0", \
                             **bg1)
a3 = patches.FancyArrowPatch((VoN[dpthInd],BS1[dpthInd,kmid]-1+VoN[dpthInd]/2.0),(VoN[dpthInd],BS0[dpthInd,kmid]+VoN[dpthInd]/2.0), \
                             connectionstyle="arc3,rad=0", \
                             **bg1)
axbs.add_patch(a2)
axbs.add_patch(a3)



bg1 = dict(arrowstyle=style, color=matica97E)
dpthInd=dpthInd-4
a2 = patches.FancyArrowPatch((VoN[dpthInd],BS0[dpthInd,klast]+1+VoN[dpthInd]/2.0), (VoN[dpthInd],BS1[dpthInd,klast]+VoN[dpthInd]/2.0), \
                             connectionstyle="arc3,rad=0", \
                             **bg1)
a3 = patches.FancyArrowPatch((VoN[dpthInd],BS1[dpthInd,klast]-1+VoN[dpthInd]/2.0),(VoN[dpthInd],BS0[dpthInd,klast]+VoN[dpthInd]/2.0), \
                             connectionstyle="arc3,rad=0", \
                             **bg1)
axbs.add_patch(a2)
axbs.add_patch(a3)






bg2 = dict(arrowstyle=style, color=matica97F)

a2 = patches.FancyArrowPatch((VoN[dpthInd],BS1[dpthInd,klast]+1+VoN[dpthInd]/2.0), (VoN[dpthInd],BS2[dpthInd,klast]+VoN[dpthInd]/2.0), \
                             connectionstyle="arc3,rad=0", \
                             **bg2)
a3 = patches.FancyArrowPatch((VoN[dpthInd],BS2[dpthInd,klast]-1+VoN[dpthInd]/2.0),(VoN[dpthInd],BS1[dpthInd,klast]+VoN[dpthInd]/2.0), \
                             connectionstyle="arc3,rad=0", \
                             **bg2)
axbs.add_patch(a2)
axbs.add_patch(a3)

dpthInd=dpthInd+4
bg2 = dict(arrowstyle=style, color=matica97F)

a2 = patches.FancyArrowPatch((VoN[dpthInd],BS1[dpthInd,kmid]+1+VoN[dpthInd]/2.0), (VoN[dpthInd],BS2[dpthInd,kmid]+VoN[dpthInd]/2.0), \
                             connectionstyle="arc3,rad=0", \
                             **bg2)
a3 = patches.FancyArrowPatch((VoN[dpthInd],BS2[dpthInd,kmid]-1+VoN[dpthInd]/2.0),(VoN[dpthInd],BS1[dpthInd,kmid]+VoN[dpthInd]/2.0), \
                             connectionstyle="arc3,rad=0", \
                             **bg2)
axbs.add_patch(a2)
axbs.add_patch(a3)


bg3 = dict(arrowstyle=style, color=matica97I)

a2 = patches.FancyArrowPatch((VoN[dpthInd],BS2[dpthInd,kmid]+1+VoN[dpthInd]/2), (VoN[dpthInd],BS3[dpthInd,kmid]+VoN[dpthInd]/2), \
                             connectionstyle="arc3,rad=0", \
                             **bg3)
a3 = patches.FancyArrowPatch((VoN[dpthInd],BS3[dpthInd,kmid]-1+VoN[dpthInd]/2),(VoN[dpthInd],BS2[dpthInd,kmid]+VoN[dpthInd]/2), \
                             connectionstyle="arc3,rad=0", \
                             **bg3)
axbs.add_patch(a2)
axbs.add_patch(a3)

bg3 = dict(arrowstyle=style, color=matica97I)
dpthInd=dpthInd-4
a2 = patches.FancyArrowPatch((VoN[dpthInd],BS2[dpthInd,klast]+1+VoN[dpthInd]/2), (VoN[dpthInd],BS3[dpthInd,klast]+VoN[dpthInd]/2), \
                             connectionstyle="arc3,rad=0", \
                             **bg3)
a3 = patches.FancyArrowPatch((VoN[dpthInd],BS3[dpthInd,klast]-1+VoN[dpthInd]/2),(VoN[dpthInd],BS2[dpthInd,klast]+VoN[dpthInd]/2), \
                             connectionstyle="arc3,rad=0", \
                             **bg3)
axbs.add_patch(a2)
axbs.add_patch(a3)

bg3 = dict(arrowstyle=style, color=matica97G)

a2 = patches.FancyArrowPatch((VoN[dpthInd],BS3[dpthInd,klast]+1+VoN[dpthInd]/2), (VoN[dpthInd],BS4[dpthInd,klast]+VoN[dpthInd]/2), \
                             connectionstyle="arc3,rad=0", \
                             **bg3)
a3 = patches.FancyArrowPatch((VoN[dpthInd],BS4[dpthInd,klast]-1+VoN[dpthInd]/2),(VoN[dpthInd],BS3[dpthInd,klast]+VoN[dpthInd]/2), \
                             connectionstyle="arc3,rad=0", \
                             **bg3)
axbs.add_patch(a2)
axbs.add_patch(a3)

bg3 = dict(arrowstyle=style, color=matica97G)

dpthInd=dpthInd+4
a2 = patches.FancyArrowPatch((VoN[dpthInd],BS3[dpthInd,kmid]+1+VoN[dpthInd]/2), (VoN[dpthInd],BS4[dpthInd,kmid]+VoN[dpthInd]/2), \
                             connectionstyle="arc3,rad=0", \
                             **bg3)
a3 = patches.FancyArrowPatch((VoN[dpthInd],BS4[dpthInd,kmid]-1+VoN[dpthInd]/2),(VoN[dpthInd],BS3[dpthInd,kmid]+VoN[dpthInd]/2), \
                             connectionstyle="arc3,rad=0", \
                             **bg3)
axbs.add_patch(a2)
axbs.add_patch(a3)

plt.savefig('BSScaleV2.pdf')
plt.show()