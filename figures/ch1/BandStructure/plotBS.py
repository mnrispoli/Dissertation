#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 12:06:01 2019

@author: matthewrispoli
"""
import numpy as np
import matplotlib.pyplot as plt 
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

VoN = np.genfromtxt('data/BSLattDepth_ref.csv', delimiter = ',')
qs = np.genfromtxt('data/BSQs_ref.csv', delimiter = ',')

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
    
dpth1=1
dpth2=47
dpth3=80

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
         
dpthext=31

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
                
FX=6.5
FY=FX/2 #/1.618

fig = plt.figure(figsize=(FX, FY))


[XL,XH,YL,YH] = [-2.0,2.0,-4.0,16.0]
[XW,YT]=[XH-XL,YH-YL]

W=4*6.5/7
H=FY
AR=XW/YT
ax2 = fig.add_subplot(121, aspect=str(AR))

listi=np.array(range(64))

ax2.plot(qsext[listi+64*0],BSext[listi+64*0], zorder = 2, color = matica97D)
ax2.plot(qsext[listi+64*1],BSext[listi+64*1], zorder = 2, color = matica97C)
ax2.plot(qsext[listi+64*2],BSext[listi+64*2], zorder = 2, color = matica97B)
ax2.plot(qsext[listi+64*3],BSext[listi+64*3], zorder = 2, color = matica97A)

ax2.plot(qsext[listi+64*4],BSext[listi+64*4], zorder = 2, color = matica97A ,\
         label = 'n = 0')
ax2.plot(qsext[listi+64*5],BSext[listi+64*5], zorder = 2, color = matica97B, \
         label = 'n = 1')
ax2.plot(qsext[listi+64*6],BSext[listi+64*6], zorder = 2, color = matica97C, \
         label = 'n = 2')
ax2.plot(qsext[listi+64*7],BSext[listi+64*7], zorder = 2, color = matica97D, \
         label = 'n = 3')

ax2.plot(qsext,4*(qsext)**2,'-', color = 'gray', alpha = 0.75, zorder = 1, \
         label = 'Free Particle')

ax2.legend(frameon=False, loc=(0,-0.3), ncol=5)

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

ax2.plot([XL,XH],[VoN[dpthext]/2, VoN[dpthext]/2],'--', color='gray' , zorder = -90)
ax2.plot([XL,XH],[-VoN[dpthext]/2, -VoN[dpthext]/2],'--', color='gray' , zorder = -90)

ax2.set_title('Dispersion Curve: '+str(int(VoN[dpthext]))+' Er')


ax2.set_ylim((YL,YH))
ax2.set_xlim((XL,XH))
ax2.set_yticks(np.linspace(YL,YH,6))
ax2.set_xticks(np.linspace(XL,XH,9))

ax2.set_xlabel('Momentum q (\hbar k))')
ax2.set_ylabel('E_k (Er)')


qv = -0.5
xv = 0
duh = 0.5
style="Simple,tail_width=0.75,head_width=4,head_length=4"

bg1 = dict(arrowstyle=style, color=matica97E)

a2 = patches.FancyArrowPatch((qv,BS0[dpthext,xv]+duh), (qv,BS1[dpthext,xv]), \
                             connectionstyle="arc3,rad=0", \
                             **bg1)
a3 = patches.FancyArrowPatch((qv,BS1[dpthext,xv]-duh),(qv,BS0[dpthext,xv]), \
                             connectionstyle="arc3,rad=0", \
                             **bg1)
ax2.add_patch(a2)
ax2.add_patch(a3)

qv = 0.5
xv = 0
duh = 0.5
style="Simple,tail_width=0.75,head_width=4,head_length=4"

bg1 = dict(arrowstyle=style, color=matica97E)

a2 = patches.FancyArrowPatch((qv,BS0[dpthext,xv]+duh), (qv,BS1[dpthext,xv]), \
                             connectionstyle="arc3,rad=0", zorder = 10, \
                             **bg1)
a3 = patches.FancyArrowPatch((qv,BS1[dpthext,xv]-duh),(qv,BS0[dpthext,xv]), \
                             connectionstyle="arc3,rad=0", zorder = 10, \
                             **bg1)
ax2.add_patch(a2)
ax2.add_patch(a3)


qv = -1
xv = kmid
duh = 0.5
style="Simple,tail_width=0.75,head_width=4,head_length=4"

bg1 = dict(arrowstyle=style, color=matica97F)

a2 = patches.FancyArrowPatch((qv,BS1[dpthext,xv]+duh), (qv,BS2[dpthext,xv]), \
                             connectionstyle="arc3,rad=0", zorder = 10, \
                             **bg1)
a3 = patches.FancyArrowPatch((qv,BS2[dpthext,xv]-duh),(qv,BS1[dpthext,xv]), \
                             connectionstyle="arc3,rad=0", zorder = 10, \
                             **bg1)
ax2.add_patch(a2)
ax2.add_patch(a3)


qv = 1
xv = kmid
duh = 0.5
style="Simple,tail_width=0.75,head_width=4,head_length=4"

bg1 = dict(arrowstyle=style, color=matica97F)

a2 = patches.FancyArrowPatch((qv,BS1[dpthext,xv]+duh), (qv,BS2[dpthext,xv]), \
                             connectionstyle="arc3,rad=0", zorder = 10, \
                             **bg1)
a3 = patches.FancyArrowPatch((qv,BS2[dpthext,xv]-duh),(qv,BS1[dpthext,xv]), \
                             connectionstyle="arc3,rad=0", zorder = 10, \
                             **bg1)
ax2.add_patch(a2)
ax2.add_patch(a3)
#plt.semilogy(VoN,BW3)
#plt.Axes.set_xticks(np.linspace(0,45,10))


[XL,XH,YL,YH] = [-0.5,0.5,-4.0,16.0]
[XW,YT]=[XH-XL,YH-YL]

W=2*6.5/7
H=FY
AR=XW/YT*1

ax122 = fig.add_subplot(122, aspect=str(AR))

ax122.plot(qs,BS0[dpthext,::], color = matica97A)
ax122.plot(qs,BS1[dpthext,::], color = matica97B)
ax122.plot(qs,BS2[dpthext,::], color = matica97C)
ax122.plot(qs,BS3[dpthext,::], color = matica97D)
ax122.set_ylim((YL,YH))
ax122.set_xlim((XL,XH))
ax122.set_yticks(np.linspace(YL,YH,6))
ax122.set_xticks(np.linspace(XL,XH,5))
ax122.set_xlabel('q')
ax122.set_title('BZ:'+str(int(VoN[dpthext]))+' Er')


ax122.plot([XL,XH],[VoN[dpthext]/2, VoN[dpthext]/2],'--', color='gray' , zorder = -90)
ax122.plot([XL,XH],[-VoN[dpthext]/2, -VoN[dpthext]/2],'--', color='gray' , zorder = -90)

# left side of diagram
xv = 3
duh = 0.5
style="Simple,tail_width=0.75,head_width=4,head_length=4"

bg1 = dict(arrowstyle=style, color=matica97E)

a2 = patches.FancyArrowPatch((qs[xv],BS0[dpthext,xv]+duh), (qs[xv],BS1[dpthext,xv]), \
                             connectionstyle="arc3,rad=0", \
                             **bg1)
a3 = patches.FancyArrowPatch((qs[xv],BS1[dpthext,xv]-duh),(qs[xv],BS0[dpthext,xv]), \
                             connectionstyle="arc3,rad=0", \
                             **bg1)
ax122.add_patch(a2)
ax122.add_patch(a3)
#
#bg2 = dict(arrowstyle=style, color=matica97F)
#
#
#a2 = patches.FancyArrowPatch((qs[xv],BS1[dpthext,xv]+1), (qs[xv],BS2[dpthext,xv]), \
#                             connectionstyle="arc3,rad=0", \
#                             **bg2)
#a3 = patches.FancyArrowPatch((qs[xv],BS2[dpthext,xv]-1),(qs[xv],BS1[dpthext,xv]), \
#                             connectionstyle="arc3,rad=0", \
#                             **bg2)
#ax163.add_patch(a2)
#ax163.add_patch(a3)

# right side of diagram

xv = 128-4

bg1 = dict(arrowstyle=style, color=matica97E)

a2 = patches.FancyArrowPatch((qs[xv],BS0[dpthext,xv]+duh), (qs[xv],BS1[dpthext,xv]), \
                             connectionstyle="arc3,rad=0", \
                             **bg1)
a3 = patches.FancyArrowPatch((qs[xv],BS1[dpthext,xv]-duh),(qs[xv],BS0[dpthext,xv]), \
                             connectionstyle="arc3,rad=0", \
                             **bg1)
ax122.add_patch(a2)
ax122.add_patch(a3)

xv=kmid

bg2 = dict(arrowstyle=style, color=matica97F)


a2 = patches.FancyArrowPatch((qs[xv],BS1[dpthext,xv]+duh), (qs[xv],BS2[dpthext,xv]), \
                             connectionstyle="arc3,rad=0", \
                             **bg2)
a3 = patches.FancyArrowPatch((qs[xv],BS2[dpthext,xv]-duh),(qs[xv],BS1[dpthext,xv]), \
                             connectionstyle="arc3,rad=0", \
                             **bg2)
ax122.add_patch(a2)
ax122.add_patch(a3)




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

plt.savefig('BZ_BS.pdf')