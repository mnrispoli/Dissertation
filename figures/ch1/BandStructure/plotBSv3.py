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
qs=np.array(qs)*2

BS0 = np.genfromtxt('data/BS0_ref.csv', delimiter = ',')
BS1 = np.genfromtxt('data/BS1_ref.csv', delimiter = ',')
BS2 = np.genfromtxt('data/BS2_ref.csv', delimiter = ',')
BS3 = np.genfromtxt('data/BS3_ref.csv', delimiter = ',')

#difference in energies between bands at same q
dBS01 = BS1-BS0
dBS12 = BS2-BS1
dBS23 = BS3-BS2

# bandwidths for the first 4 bands
BW0 = np.ones(len(VoN))
BW1 = np.ones(len(VoN))
BW2 = np.ones(len(VoN))
BW3 = np.ones(len(VoN))

# bandgaps between the first 4 bands
BG01min = np.ones(len(VoN))
BG12min = np.ones(len(VoN))
BG23min = np.ones(len(VoN))

BG01max = np.ones(len(VoN))
BG12max = np.ones(len(VoN))
BG23max = np.ones(len(VoN))

for nn in range(np.size(BS0,0)):
    BW0[nn]=np.amax(BS0[nn,::])-np.amin(BS0[nn,::])
    BW1[nn]=np.amax(BS1[nn,::])-np.amin(BS1[nn,::])
    BW2[nn]=np.amax(BS2[nn,::])-np.amin(BS2[nn,::])
    BW3[nn]=np.amax(BS3[nn,::])-np.amin(BS3[nn,::])
    
    BG01min[nn]=np.amin(dBS01[nn,::])
    BG12min[nn]=np.amin(dBS12[nn,::])
    BG23min[nn]=np.amin(dBS23[nn,::])
    
    BG01max[nn]=np.amax(dBS01[nn,::])
    BG12max[nn]=np.amax(dBS12[nn,::])
    BG23max[nn]=np.amax(dBS23[nn,::])

#plt.semilogy(VoN,BW0)
#plt.semilogy(VoN,BW1)
#plt.semilogy(VoN,BW2)
#plt.semilogy(VoN,BW3)
    
#plt.loglog([1,2,3,4],[BW0[0],BW1[0],BW2[0],BW3[0]])
    
dpth1=31-4*5
dpth2=31+4
dpth3=31+4*5+3
dpthext=dpth2

kmid=64
klast=127

qsext=np.hstack(\
                (qs[kmid:klast]-2, \
                 np.hstack(\
                    (\
                     np.hstack(\
                    (np.hstack((qs-1,qs)),qs+1)\
                    ),
                    qs[0:kmid]+2 )\
                 )))
                     
qsext=np.hstack(\
                (\
                qs[kmid::]-2,\
                qs-1,qs,qs+1,\
                qs[0:kmid]+2)\
                )
         


#dpthext=3
            
BSext=np.hstack(\
                (BS3[dpthext,kmid::], \
                 BS2[dpthext,0:kmid],\
                 BS1[dpthext,kmid::], \
                 BS0[dpthext,::], BS1[dpthext,0:kmid],
                 BS2[dpthext,kmid::], \
                 BS3[dpthext,0:kmid] \
                 )\
                 )-0*np.amin(BS0[dpthext,::])/2
                 
pltOffSt=VoN[dpthext]/2
                
#qsext=qsext*2

FX=6.5*2
FY=FX/2 #/1.618

fig = plt.figure(figsize=(FX, FY))


[XL,XH,YL,YH] = [-4.0,4.0,-4.0+pltOffSt,16.0+pltOffSt]
[XW,YT]=[XH-XL,YH-YL]

W=4*6.5/7
H=FY
AR=XW/YT
ax2 = fig.add_subplot(141, aspect='0.24')

listi=np.array(range(64))

ax2.plot(qsext[listi+64*0],BSext[listi+64*0]+pltOffSt, zorder = 2, color = band3Clr)
ax2.plot(qsext[listi+64*1],BSext[listi+64*1]+pltOffSt, zorder = 2, color = band2Clr)
ax2.plot(qsext[listi+64*2],BSext[listi+64*2]+pltOffSt, zorder = 2, color = band1Clr)
ax2.plot(qsext[listi+64*3],BSext[listi+64*3]+pltOffSt, zorder = 2, color = band0Clr)

ax2.plot(qsext[listi+64*4],BSext[listi+64*4]+pltOffSt, zorder = 2, color = band0Clr ,\
         label = 'n = 0')
ax2.plot(qsext[listi+64*5],BSext[listi+64*5]+pltOffSt, zorder = 2, color = band1Clr, \
         label = 'n = 1')
ax2.plot(qsext[listi+64*6],BSext[listi+64*6]+pltOffSt, zorder = 2, color = band2Clr, \
         label = 'n = 2')
ax2.plot(qsext[listi+64*7],BSext[listi+64*7]+pltOffSt, zorder = 2, color = band3Clr, \
         label = 'n = 3')

ax2.plot(qsext,4*(qsext)**2+pltOffSt,'-', color = 'gray', alpha = 0.75, zorder = 1, \
         label = 'Free Particle')



grlw=0.5
gralpha=0.2
#faux gridlines
ax2.plot([-1.5,-1.5],[YL,YH], 'k-', linewidth = grlw, alpha = gralpha , zorder = -100)
#ax2.plot([-1.5,-1.5],[-4,-3.5], 'k-', linewidth = 1, alpha = 0.5 , zorder = -100)

ax2.plot([-1,-1],[YL,YH], 'k-', linewidth = grlw, alpha = gralpha , zorder = -100)
#ax2.plot([-1,-1],[-4,-3.5], 'k-', linewidth = 1, alpha = 0.5 , zorder = -100)

ax2.plot([-0.5,-0.5],[YL,YH], 'k-', linewidth = grlw, alpha = gralpha , zorder = -100)
ax2.plot([0,0],[YL,YH], 'k-', linewidth = grlw, alpha = gralpha , zorder = -100)
ax2.plot([0.5,0.5],[YL,YH], 'k-', linewidth = grlw, alpha = gralpha , zorder = -100)
ax2.plot([1.0,1.0],[YL,YH], 'k-', linewidth = grlw, alpha = gralpha , zorder = -100)
ax2.plot([1.5,1.5],[YL,YH], 'k-', linewidth = grlw, alpha = gralpha , zorder = -100)

ax2.plot([XL,XH],[VoN[dpthext]/2+pltOffSt, VoN[dpthext]/2+pltOffSt],'--', color='gray' , zorder = -90, label='Lattice Depth')
ax2.plot([XL,XH],[-VoN[dpthext]/2+pltOffSt, -VoN[dpthext]/2+pltOffSt],'--', color='gray' , zorder = -90)

ax2.set_title('Dispersion Curve: '+str(int(VoN[dpthext]))+' Er')


ax2.set_ylim((0,20))
ax2.set_xlim((-4,4))
ax2.set_yticks(np.linspace(0,20,6))
ax2.set_xticks(np.linspace(-4,4,5))

ax2.set_xlabel('Momentum q (\hbar k))')
ax2.set_ylabel('E_k (Er)')

ax2.legend(frameon=False, loc=(1,-0.3), ncol=6)

qv = -0.5
xv = 0
duh = 0.5
style="Simple,tail_width=0.75,head_width=4,head_length=4"

bg1 = dict(arrowstyle=style, color=matica97E)

a2 = patches.FancyArrowPatch((qv,BS0[dpthext,xv]+duh+pltOffSt), (qv,BS1[dpthext,xv]+pltOffSt), \
                             connectionstyle="arc3,rad=0", \
                             **bg1)
a3 = patches.FancyArrowPatch((qv,BS1[dpthext,xv]-duh+pltOffSt),(qv,BS0[dpthext,xv]+pltOffSt), \
                             connectionstyle="arc3,rad=0", \
                             **bg1)
ax2.add_patch(a2)
ax2.add_patch(a3)

qv = 0.5
xv = 0
duh = 0.5
style="Simple,tail_width=0.75,head_width=4,head_length=4"

bg1 = dict(arrowstyle=style, color=matica97E)

a2 = patches.FancyArrowPatch((qv,BS0[dpthext,xv]+duh+pltOffSt), (qv,BS1[dpthext,xv]+pltOffSt), \
                             connectionstyle="arc3,rad=0", zorder = 10, \
                             **bg1)
a3 = patches.FancyArrowPatch((qv,BS1[dpthext,xv]-duh+pltOffSt),(qv,BS0[dpthext,xv]+pltOffSt), \
                             connectionstyle="arc3,rad=0", zorder = 10, \
                             **bg1)
ax2.add_patch(a2)
ax2.add_patch(a3)


qv = -1
xv = kmid
duh = 0.5
style="Simple,tail_width=0.75,head_width=4,head_length=4"

bg1 = dict(arrowstyle=style, color=matica97F)

a2 = patches.FancyArrowPatch((qv,BS1[dpthext,xv]+duh+pltOffSt), (qv,BS2[dpthext,xv]+pltOffSt), \
                             connectionstyle="arc3,rad=0", zorder = 10, \
                             **bg1)
a3 = patches.FancyArrowPatch((qv,BS2[dpthext,xv]-duh+pltOffSt),(qv,BS1[dpthext,xv]+pltOffSt), \
                             connectionstyle="arc3,rad=0", zorder = 10, \
                             **bg1)
ax2.add_patch(a2)
ax2.add_patch(a3)


qv = 1
xv = kmid
duh = 0.5
style="Simple,tail_width=0.75,head_width=4,head_length=4"

bg1 = dict(arrowstyle=style, color=matica97F)

a2 = patches.FancyArrowPatch((qv,BS1[dpthext,xv]+duh+pltOffSt), (qv,BS2[dpthext,xv]+pltOffSt), \
                             connectionstyle="arc3,rad=0", zorder = 10, \
                             **bg1)
a3 = patches.FancyArrowPatch((qv,BS2[dpthext,xv]-duh+pltOffSt),(qv,BS1[dpthext,xv]+pltOffSt), \
                             connectionstyle="arc3,rad=0", zorder = 10, \
                             **bg1)
ax2.add_patch(a2)
ax2.add_patch(a3)
#plt.semilogy(VoN,BW3)
#plt.Axes.set_xticks(np.linspace(0,45,10))


[XL,XH,YL,YH] = [-0.5,0.5,-4.0+pltOffSt,16.0+pltOffSt]
[XW,YT]=[XH-XL,YH-YL]

W=2*6.5/7
H=FY
AR=XW/YT*1
#
#ax122 = fig.add_subplot(122, aspect='0.06')
#
#ax122.plot(qs,BS0[dpthext,::]+pltOffSt, color = band0Clr)
#ax122.plot(qs,BS1[dpthext,::]+pltOffSt, color = band1Clr)
#ax122.plot(qs,BS2[dpthext,::]+pltOffSt, color = band2Clr)
#ax122.plot(qs,BS3[dpthext,::]+pltOffSt, color = band3Clr)
#ax122.set_ylim((YL,YH))
#ax122.set_xlim((XL,XH))
#ax122.set_yticks(np.linspace(YL,YH,6))
#ax122.set_xticks(np.linspace(XL,XH,5))
#ax122.set_xlabel('q')
#ax122.set_title('BZ:'+str(int(VoN[dpthext]))+' Er')
#
#
#ax122.plot([XL,XH],[VoN[dpthext]/2+pltOffSt, VoN[dpthext]/2+pltOffSt],'--', color='gray' , zorder = -90)
#ax122.plot([XL,XH],[-VoN[dpthext]/2+pltOffSt, -VoN[dpthext]/2+pltOffSt],'--', color='gray' , zorder = -90)
#
## left side of diagram
#xv = 3
#duh = 0.5
#style="Simple,tail_width=0.75,head_width=4,head_length=4"
#
#bg1 = dict(arrowstyle=style, color=matica97E)
#
#a2 = patches.FancyArrowPatch((qs[xv],BS0[dpthext,xv]+duh+pltOffSt), (qs[xv],BS1[dpthext,xv]+pltOffSt), \
#                             connectionstyle="arc3,rad=0", \
#                             **bg1)
#a3 = patches.FancyArrowPatch((qs[xv],BS1[dpthext,xv]-duh+pltOffSt),(qs[xv],BS0[dpthext,xv]+pltOffSt), \
#                             connectionstyle="arc3,rad=0", \
#                             **bg1)
#ax122.add_patch(a2)
#ax122.add_patch(a3)
#
#
#
## right side of diagram




#### new BZ

ax161 = fig.add_subplot(142, aspect='0.06')

ax161.plot(qs,BS0[dpth1,::]+VoN[dpth1]/2, color = band0Clr, label = 'n = 0')
#ax161.fill_between(qs, BS0[dpth1,::]+VoN[dpth1]/2, \
#                 BS1[dpth1,::]+VoN[dpth1]/2, \
#                 color = blendClr(matica97E,clrWhite,0.35), \
#                 zorder = -3)
ax161.plot(qs,BS1[dpth1,::]+VoN[dpth1]/2, color = band1Clr, label = 'n = 1')
#ax161.fill_between(qs, BS1[dpth1,::]+VoN[dpth1]/2, \
#                 BS2[dpth1,::]+VoN[dpth1]/2, \
#                 color = blendClr(matica97F,clrWhite,0.35), \
#                 zorder = -3)
ax161.plot(qs,BS2[dpth1,::]+VoN[dpth1]/2, color = band2Clr, label = 'n = 2')
#ax161.fill_between(qs, BS2[dpth1,::]+VoN[dpth1]/2, \
#                 BS3[dpth1,::]+VoN[dpth1]/2, \
#                 color = blendClr(matica97I,clrWhite,0.35), \
#                 zorder = -3)
ax161.plot(qs,BS3[dpth1,::]+VoN[dpth1]/2, color = band3Clr, label = 'n = 3')
#ax161.plot(qs,BS4[dpth1,::]+VoN[dpth1]/2, color = band4Clr, label = 'n = 4')
#ax161.plot(qs,BS5[dpth1,::]+VoN[dpth1]/2, color = band5Clr, label = 'n = 5')
ax161.plot(qs,np.ones(qs.shape)*VoN[dpth1],ls='--',color = 'gray')
ax161.set_ylim((0,20))
ax161.set_xlim((-0.5,0.5))
ax161.set_yticks([])
#ax161.set_ylabel('E_q (Er)')
ax161.set_xlabel('q')
ax161.set_title(str(VoN[dpth1])+' Er')
ax161.set_xticks(np.linspace(-0.5,0.5,5))
#ax161.legend(frameon=False)

ax162 = fig.add_subplot(143, aspect='0.06')

ax162.plot(qs,BS0[dpth2,::]+VoN[dpth2]/2, color = band0Clr)
ax162.plot(qs,BS1[dpth2,::]+VoN[dpth2]/2, color = band1Clr)
ax162.plot(qs,BS2[dpth2,::]+VoN[dpth2]/2, color = band2Clr)
ax162.plot(qs,BS3[dpth2,::]+VoN[dpth2]/2, color = band3Clr)
#ax162.plot(qs,BS4[dpth2,::]+VoN[dpth2]/2, color = band4Clr)
ax162.plot(qs,np.ones(qs.shape)*VoN[dpth2],ls='--',color = 'gray')
ax162.set_ylim((0,20))
ax162.set_xlim((-0.5,0.5))
ax162.set_yticks([])
ax162.set_xticks(np.linspace(-0.5,0.5,5))
ax162.set_xlabel('q')
ax162.set_title(str(VoN[dpth2])+' Er')

ax163 = fig.add_subplot(144, aspect='0.06')

ax163.plot(qs,BS0[dpth3,::]+VoN[dpth3]/2, color = band0Clr)
ax163.plot(qs,BS1[dpth3,::]+VoN[dpth3]/2, color = band1Clr)
ax163.plot(qs,BS2[dpth3,::]+VoN[dpth3]/2, color = band2Clr)
ax163.plot(qs,BS3[dpth3,::]+VoN[dpth3]/2, color = band3Clr)
ax163.plot(qs,np.ones(qs.shape)*VoN[dpth3],ls='--',color = 'gray')
ax163.set_ylim((0,20))
ax163.set_xlim((-0.5,0.5))
ax163.set_yticks([])
ax163.set_xticks(np.linspace(-0.5,0.5,5))
ax163.set_xlabel('q')
ax163.set_title(str(VoN[dpth3])+' Er')








xv = 128-4

bg1 = dict(arrowstyle=style, color=matica97E)

a2 = patches.FancyArrowPatch((qs[xv],BS0[dpthext,xv]+duh+pltOffSt), (qs[xv],BS1[dpthext,xv]+pltOffSt), \
                             connectionstyle="arc3,rad=0", \
                             **bg1)
a3 = patches.FancyArrowPatch((qs[xv],BS1[dpthext,xv]-duh+pltOffSt),(qs[xv],BS0[dpthext,xv]+pltOffSt), \
                             connectionstyle="arc3,rad=0", \
                             **bg1)
ax162.add_patch(a2)
ax162.add_patch(a3)

xv=4
bg1 = dict(arrowstyle=style, color=matica97E)

a2 = patches.FancyArrowPatch((qs[xv],BS0[dpthext,xv]+duh+pltOffSt), (qs[xv],BS1[dpthext,xv]+pltOffSt), \
                             connectionstyle="arc3,rad=0", \
                             **bg1)
a3 = patches.FancyArrowPatch((qs[xv],BS1[dpthext,xv]-duh+pltOffSt),(qs[xv],BS0[dpthext,xv]+pltOffSt), \
                             connectionstyle="arc3,rad=0", \
                             **bg1)
ax162.add_patch(a2)
ax162.add_patch(a3)



xv=kmid

bg2 = dict(arrowstyle=style, color=matica97F)


a2 = patches.FancyArrowPatch((qs[xv],BS1[dpthext,xv]+duh+pltOffSt), (qs[xv],BS2[dpthext,xv]+pltOffSt), \
                             connectionstyle="arc3,rad=0", \
                             **bg2)
a3 = patches.FancyArrowPatch((qs[xv],BS2[dpthext,xv]-duh+pltOffSt),(qs[xv],BS1[dpthext,xv]+pltOffSt), \
                             connectionstyle="arc3,rad=0", \
                             **bg2)
ax162.add_patch(a2)
ax162.add_patch(a3)

#ax2 = fig.add_subplot(122, aspect='2.3176')
#
#ax2 = fig.add_subplot(122, aspect='auto')  
#ax2.plot(VoN,BG01min,'-',color = matica97E, zorder = 3, label = 'n=0->n=1')
#ax2.plot(VoN,BG01max,'-',color = matica97E, zorder = 3)
#ax2.fill_between(VoN, BG01min, BG01max, \
#                 color = blendClr(matica97E,clrWhite,0.35), \
#                 zorder = 3)
#
#ax2.plot(VoN,BG12min,'-', color = matica97F, zorder = 2, label = 'n=1->n=2' )
#ax2.plot(VoN,BG12max,'-', color = matica97F, zorder = 2 )
#ax2.fill_between(VoN, BG12min, BG12max, \
#                 color = blendClr(matica97F,clrWhite,0.35), \
#                 zorder = 2)
#
#ax2.plot(VoN,BG23min,'-', color = matica97I, zorder = 1, label = 'n=2->n=3' )
#ax2.plot(VoN,BG23max,'-', color = matica97I, zorder = 1 )
#ax2.fill_between(VoN, BG23min, BG23max, \
#                 color = blendClr(matica97I,clrWhite,0.35), \
#                 zorder = 1)
#
#ax2.plot(VoN,BG34min,'-', color = matica97G, zorder = 0, label = 'n=3->n=4' )
#ax2.plot(VoN,BG34max,'-', color = matica97G, zorder = 0 )
#ax2.fill_between(VoN, BG34min, BG34max, \
#                 color = blendClr(matica97G,clrWhite,0.35), \
#                 zorder = 0)

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

#ax2.set_xticks(np.linspace(0,50,11))
#ax2.set_yticks(np.linspace(0,16,5))
#ax2.set_xlim((0.24,50))
#ax2.set_ylim((0,30))
#
#ax2.set_xlabel('Lattice Depth (Er)')
#ax2.set_ylabel('Direct Band Gap (Er)')
#ax2.legend(frameon=False, loc=(0))
#plt.semilogy(VoN,BW3)
#plt.Axes.set_xticks(np.linspace(0,45,10))

#kmid=32
#klast=63
#
#style="Simple,tail_width=1,head_width=5,head_length=6"
#bg1 = dict(arrowstyle=style, color=matica97E)
#
#a2 = patches.FancyArrowPatch((0,BS0[dpth3,kmid]+1+VoN[dpth3]/2.0), (0,BS1[dpth3,kmid]+VoN[dpth3]/2.0), \
#                             connectionstyle="arc3,rad=0", \
#                             **bg1)
#a3 = patches.FancyArrowPatch((0,BS1[dpth3,kmid]-1+VoN[dpth3]/2.0),(0,BS0[dpth3,kmid]+VoN[dpth3]/2.0), \
#                             connectionstyle="arc3,rad=0", \
#                             **bg1)
#ax163.add_patch(a2)
#ax163.add_patch(a3)
#
#bg2 = dict(arrowstyle=style, color=matica97F)
#
#a2 = patches.FancyArrowPatch((0,BS1[dpth3,kmid]+1+VoN[dpth3]/2.0), (0,BS2[dpth3,kmid]+VoN[dpth3]/2.0), \
#                             connectionstyle="arc3,rad=0", \
#                             **bg2)
#a3 = patches.FancyArrowPatch((0,BS2[dpth3,kmid]-1+VoN[dpth3]/2.0),(0,BS1[dpth3,kmid]+VoN[dpth3]/2.0), \
#                             connectionstyle="arc3,rad=0", \
#                             **bg2)
#ax163.add_patch(a2)
#ax163.add_patch(a3)
#
#
#bg3 = dict(arrowstyle=style, color=matica97I)
#
#a2 = patches.FancyArrowPatch((0,BS2[dpth3,kmid]+1+VoN[dpth3]/2), (0,BS3[dpth3,kmid]+VoN[dpth3]/2), \
#                             connectionstyle="arc3,rad=0", \
#                             **bg3)
#a3 = patches.FancyArrowPatch((0,BS3[dpth3,kmid]-1+VoN[dpth3]/2),(0,BS2[dpth3,kmid]+VoN[dpth3]/2), \
#                             connectionstyle="arc3,rad=0", \
#                             **bg3)
#ax163.add_patch(a2)
#ax163.add_patch(a3)




#
#bg3 = dict(arrowstyle=style, color=matica97I)
#
#a2 = patches.FancyArrowPatch((qs[xv],BS2[dpthext,xv]+1), (qs[xv],BS3[dpthext,xv]), \
#                             connectionstyle="arc3,rad=0", \
#                             **bg3)
#a3 = patches.FancyArrowPatch((qs[xv],BS3[dpthext,xv]-1),(qs[xv],BS2[dpthext,xv]), \
#                             connectionstyle="arc3,rad=0", \
#                             **bg3)
#ax163.add_patch(a2)
#ax163.add_patch(a3)

plt.savefig('BZ_BS_v3.pdf')