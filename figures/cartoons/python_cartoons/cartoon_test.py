# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 15:23:56 2019

@author: Matthew
"""

import matplotlib.pyplot as plt
from matplotlib import patches

import numpy as np
import pandas as pd

from cartoon_plot_settings import rcpars as rcpars_crtn

def setPltPars(rcparsI):
    for key_name, value_name in rcparsI.items():
        plt.rcParams[key_name]=value_name
        
def cm2inch(value):
    return value/2.54

def inch2cm(value):
    return value*2.54


cm = 1/2.54

#point offset for clarity
dx=0.05
titleFS = 10
axisFS = 8
annFS = 5

setPltPars(rcpars_crtn)




xcenter, ycenter = 0.38*cm, 0.52*cm
width, height = 1e-1*cm, 3e-1*cm
angle = -30

theta = np.deg2rad(np.arange(0.0, 360.0, 1.0))
x = 0.5 * width * np.cos(theta)
y = 0.5 * height * np.sin(theta)

rtheta = np.radians(angle)
R = np.array([
    [np.cos(rtheta), -np.sin(rtheta)],
    [np.sin(rtheta),  np.cos(rtheta)],
    ])


x, y = np.dot(R, np.array([x, y]))
x += xcenter
y += ycenter




                
                
fig = plt.figure()

ax = fig.add_subplot(211, aspect='auto')

ax.fill(x, y, alpha=1, facecolor = (1,0,0),
        edgecolor='yellow', linewidth=1, zorder=1)

e1 = patches.Ellipse((xcenter, ycenter), width, height,
                     angle=angle, linewidth=2, fill=False, zorder=2)

ax.add_patch(e1)



                
ax = fig.add_subplot(212, aspect='equal')

ax.fill(x, y, alpha=1, facecolor='#a50f15', edgecolor='#a50f15', zorder=1)
        
e2 = patches.Ellipse((xcenter, ycenter), width, height,
                     angle=angle, linewidth=2, fill=False, zorder=2)

a1 = patches.Arc((xcenter, ycenter), width/2, height/2, \
                 -30+90, 0, 180,\
                 )

e1.set_hatch('/')
e1.set_linestyle('-.')

ax.add_patch(a1)
ax.add_patch(e2)

xx=np.linspace(0,0.4,100);
yy=np.cos(xx*2*np.pi)/4+0.2
ax.plot(xx,yy)



#fig = plt.figure(figsize=(cm2inch(8.6), cm2inch(8.6/2 + 8.6/1.618)))
#
#ax1 = fig.add_subplot(1,1,1)


#
#ax1.imshow(g2data, cmap = 'bwr', vmin = -1, vmax = 1)
#ax1.axes.set_xlabel("site-j", fontsize = axisFS)
#ax1.axes.set_ylabel("site-i", fontsize = axisFS)
#
#
#
#
#ax1.axes.set_xticks((0,5.5,11))
#ax1.axes.set_xticklabels((-6,0,6))
#
#ax1.axes.set_yticks((0,5.5,11))
#ax1.axes.set_yticklabels((-6,0,6))
#
#ax1.plot([5.5,5.5],[-0.5,11.5],'k--', linewidth=0.8)
#ax1.plot([-0.5,11.5],[5.5,5.5],'k--', linewidth=0.8)
#
#ax1.axes.set_title(r"$G_2(i,j)$ Data", \
#                   fontsize=titleFS, verticalalignment='bottom')


plt.savefig('test_cartoon.pdf')
#plt.show()