#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 15:03:02 2019

@author: matthewrispoli
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


def ujC(mu,n):
    return (mu-n+1)*(n-mu)/((n+1)*(mu-n+1)+(n-mu)*n)


                
FX=6.5*2
FY=FX/3#/1.618

[XL,XH,YL,YH] = [0,0.25,0,4.0]
[XW,YT]=[XH-XL,YH-YL]

W=XW/2
H=FY
AR=XW/YT

fig = plt.figure(figsize=(FX, FY))

ax1 = fig.add_subplot(121, aspect='auto')

axbs_alpha = 0.4

dt = 0.1
dx = 0.0075

mus1=np.linspace(0,1,199)
m1list=list(ujC(mus1,1))
m1m=m1list.index(max(m1list))

mus1b = np.linspace(0,mus1[m1m],100)
mus1t = np.linspace(1,mus1[m1m],100)

ax1.plot(ujC(mus1,1),mus1, ls='-', color = 'black')

ax1.fill_between(ujC(mus1t,1), mus1[m1m], mus1t, \
                 color = blendClr(clrWhite,clrBlack,0.8), \
                 zorder = 0)

ax1.fill_between(ujC(mus1b,1), mus1[m1m], mus1b, \
                 color = blendClr(clrWhite,clrBlack,0.8), \
                 zorder = 0)

ax1.annotate("MI (n=1)", xytext=(dx,0.5-dt), xy=(0,0) \
            )


mus2=np.linspace(1,2,199)
m2list=list(ujC(mus2,2))
m2m=m2list.index(max(m2list))

mus2b = np.linspace(1,mus2[m2m],100)
mus2t = np.linspace(2,mus2[m2m],100)

ax1.plot(ujC(mus2,2),mus2, ls='-', color = 'black')

ax1.fill_between(ujC(mus2b,2), mus2[m2m], mus2t, \
                 color = blendClr(clrWhite,clrBlack,0.8), \
                 zorder = 0)

ax1.fill_between(ujC(mus2b,2), mus2b, mus2[m2m], \
                 color = blendClr(clrWhite,clrBlack,0.8), \
                 zorder = 0)

ax1.annotate("MI (n=2)", xytext=(dx,1.5-dt), xy=(0,0) \
            )



mus3=np.linspace(2,3,100)

mus3b = np.linspace(2,2.45,100)
mus3t = np.linspace(3,2.45,100)

ax1.plot(ujC(mus3,3),mus3, ls='-', color = 'black')

ax1.fill_between(ujC(mus3b,3), mus3b, mus3t, \
                 color = 'gray', \
                 alpha = axbs_alpha, \
                 zorder = 0)

ax1.annotate("MI (n=3)", xytext=(dx,2.5-dt), xy=(0,0) \
            )

mus4=np.linspace(3,4,100)

mus4b = np.linspace(3,3.45,100)
mus4t = np.linspace(4,3.45,100)

ax1.plot(ujC(mus4,4),mus4, ls='-', color = 'black')

ax1.fill_between(ujC(mus4b,4), mus4t, mus4b, \
                 color = 'gray', \
                 alpha = axbs_alpha, \
                 zorder = 0)

#ax1.fill_between(ujC(mus4b,4), 3.49, mus4b, \
#                 color = 'gray', \
#                 alpha = axbs_alpha, \
#                 zorder = 0)

ax1.annotate("MI (n=4)", xytext=(dx,3.5-dt), xy=(0,0) \
            )

ax1.annotate("SF", xytext=(0.21,0.375), xy=(0,0) \
            )

ax1.set_xlim([XL,XH])
ax1.set_ylim([YL,YH])
ax1.set_ylabel('mu/U')
ax1.set_xlabel('zJ/U')
ax1.set_title('Mean-Field Phase Diagram')

left, bottom, width, height = [0.3175, 0.4, 0.3/3*1.5, 0.3*1.5]
ax1in = fig.add_axes([left, bottom, width, height])#, aspect='equal')

nbar=np.linspace(1,100,100)
ujc = ujC(nbar-0.5,nbar)

ax1in.loglog(nbar,ujc,'o-',\
             color=blendClr(matica97A,clrWhite,0.8),\
             zorder = -1\
             
             )
ax1in.set_ylabel('(U/J)c')
ax1in.set_xlabel('mu/U')
beta=np.polyfit(np.log10(nbar),np.log10(ujc),1)
ax1in.loglog(nbar,((nbar)**beta[0])*(10**beta[1]), 'k--')
#ax1in.loglog(range(10), color='red')

ax2 = fig.add_subplot(122, aspect='auto')
ax2.plot(0,0)
ax1.set_title('add cartoons for phases')
plt.savefig('MeanFieldDiagv2.pdf')