# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 15:23:56 2019

@author: Matthew
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from thry_plot_settings import rcpars as rcparsThry

def setPltPars(rcparsI):
    for key_name, value_name in rcparsI.items():
        plt.rcParams[key_name]=value_name
        
def cm2inch(value):
    return value/2.54

def inch2cm(value):
    return value*2.54

#point offset for clarity
dx=0.05
titleFS = 10
axisFS = 8
annFS = 5

setPltPars(rcparsThry)

#plt.rc('text',usetex=True)

#import data section
g2data=pd.read_csv("g2_data.csv",names = np.arange(1,13,1)).values
g2thry=pd.read_csv("g2_theory.csv",names = np.arange(1,13,1)).values
g2datXtion=pd.read_csv("g2_data_xtion.csv",names = [0,1,2]).values
g2thryXtion=pd.read_csv("g2_theory_xtion.csv",names = [0,1]).values

xiBth5p5J=pd.read_csv("xi_bth_5.5J_time.csv",names = [0,1]).values
xiBth9p1J=pd.read_csv("xi_bth_9.1J_time.csv",names = [0,1]).values

xiBth5p5Jdat=pd.read_csv("W11_6s_time.csv",names = [0,1,2]).values
xiBth9p1Jdat=pd.read_csv("W18.4_6s_time.csv",names = [0,1,2]).values

xiBthFit=pd.read_csv("exp_decay.csv",names = [0,1]).values




fig = plt.figure(figsize=(cm2inch(8.6), cm2inch(8.6/2 + 8.6/1.618)))

ax1 = fig.add_subplot(2,2,1)

ax1.imshow(g2data, cmap = 'bwr', vmin = -1, vmax = 1)
ax1.axes.set_xlabel("site-j", fontsize = axisFS)
ax1.axes.set_ylabel("site-i", fontsize = axisFS)




ax1.axes.set_xticks((0,5.5,11))
ax1.axes.set_xticklabels((-6,0,6))

ax1.axes.set_yticks((0,5.5,11))
ax1.axes.set_yticklabels((-6,0,6))

ax1.plot([5.5,5.5],[-0.5,11.5],'k--', linewidth=0.8)
ax1.plot([-0.5,11.5],[5.5,5.5],'k--', linewidth=0.8)

ax1.axes.set_title(r"$G_2(i,j)$ Data", \
                   fontsize=titleFS, verticalalignment='bottom')




ax2 = fig.add_subplot(2,2,2)

pos1 = ax2.get_position()
pos2 = [pos1.x0 + 0.08, pos1.y0,  pos1.width / 1.0, pos1.height / 1.0]
ax2.axes.set_position(pos2)

ax2.imshow(g2thry, cmap = 'bwr', vmin = -1, vmax = 1)
ax2.axes.set_xlabel("site-j", fontsize = axisFS)
#ax2.axes.set_ylabel("site-j", fontsize = 12)

ax2.plot([5.5,5.5],[-0.5,11.5],'k--', linewidth=0.8)
ax2.plot([-0.5,11.5],[5.5,5.5],'k--', linewidth=0.8)

ax2.axes.set_xticks((0,5.5,11))
ax2.axes.set_xticklabels((-6,0,6))
ax2.axes.set_yticks((0,5.5,11))
ax2.axes.set_yticklabels((-6,0,6))
#ax2.axes.set_yticks(np.arange(0,12,2))
#ax2.axes.set_yticklabels(np.arange(1,13,2))

#ax2.axes.set_title(r"\TeX\ is Number " \
#                   r"$\displaystyle\sum_{n=1}^\infty\frac{-e^{i\pi}}{2^n}$!", \
#                   fontsize=20, verticalalignment='bottom')

ax2.axes.set_title(r"$G_2(i,j)$ Theory", \
                   fontsize=titleFS, verticalalignment='bottom')



ax3 = fig.add_subplot(2,2,3)
ax3.axes.set_aspect(12/0.15/1.618/1.6)

pos1 = ax3.get_position()
pos2 = [pos1.x0, pos1.y0 + 0.07,  pos1.width / 1.0, pos1.height / 1.0]
ax3.axes.set_position(pos2)

ax3_1=ax3.plot(g2thryXtion[:,0], g2thryXtion[:,1], '-',\
               alpha=0.6)

ax3_3=ax3.errorbar(g2datXtion[:,0], g2datXtion[:,1], yerr = g2datXtion[:,2], fmt='o',\
             markerfacecolor = 'w', \
             markersize = 4, \
             markeredgewidth = 1.5, \
             color=ax3_1[0].get_color(), label=r"W=9.1J")

ax3_2=ax3.plot(xiBthFit[:,0],xiBthFit[:,1],'r--', \
               label = "Fit")


ax3.axes.yaxis.labelpad = -1
ax3.axes.set_xticks((1,6.5,12))
ax3.axes.set_xticklabels((-6,0,6))
ax3.axes.set_yticks((-0.15,-0.10,-0.05,0))
ax3.axes.set_yticklabels((-0.15,-0.10,-0.05,0.00))

ax3.axes.set_xlabel("site-j", fontsize = axisFS)
ax3.axes.set_ylabel(r"$G_2($Bath Sites$)$", \
                    fontsize = axisFS)

ax3.grid()

ax3.annotate(r"$\xi$", \
             xy=(xiBthFit[-2,0]+dx,xiBthFit[-2,1]), \
             xytext=(6.75+0.5, xiBthFit[-2,1]-0.005) , \
             fontsize = annFS, \
             arrowprops = dict(arrowstyle="<->"))


ax4 = fig.add_subplot(2,2,4)
ax4.axes.set_aspect(np.log(250)/np.log(12)/1.618/1.6)

pos1 = ax4.get_position()
pos2 = [pos1.x0 + 0.08, pos1.y0 + 0.07,  pos1.width / 1.0, pos1.height / 1.0]
ax4.axes.set_position(pos2)

ax4_2=ax4.plot(xiBth9p1J[:,0], xiBth9p1J[:,1], '-', \
               alpha=0.6, zorder=0)

ax4_1=ax4.plot(xiBth5p5J[:,0], xiBth5p5J[:,1], '-', \
               alpha=0.6, zorder=10)

ax4_3=ax4.errorbar(xiBth9p1Jdat[:,0]-dx, xiBth9p1Jdat[:,1], yerr = xiBth9p1Jdat[:,2], fmt='o', \
             markerfacecolor = 'w', \
             markersize = 4, \
             markeredgewidth = 1.5, \
             color = ax4_2[0].get_color(), label="W=9.1J", \
             zorder = 20)      
                   

ax4_4=ax4.errorbar(xiBth5p5Jdat[:,0]+dx, xiBth5p5Jdat[:,1], yerr = xiBth5p5Jdat[:,2], fmt='o', \
             markerfacecolor = 'w', \
             markersize = 4, \
             markeredgewidth = 1.5, \
             color = ax4_1[0].get_color(), label="W=5.5J", \
             zorder = 30)

#ax4.legend()

ax4.set_xscale('log')
ax4.set_yscale('log')
ax4.axes.set_xlim(0.5,250)
ax4.axes.set_ylim(0.3,12)

ax4.axes.tick_params(axis='both', which='major')


#ax4.axes.set_xlabel(r"time ($2\pi Jt$)", \
#                    fontsize = 12)
ax4.axes.set_xlabel(r"time ($\tau$)", \
                    fontsize = axisFS)
ax4.axes.yaxis.labelpad = 0
ax4.axes.set_ylabel(r"$\xi_{Bath}$", \
                    fontsize = axisFS)

ax4.axes.tick_params(which='major',grid_alpha=0.8)
ax4.axes.tick_params(which='minor',grid_alpha=0.3)

ax4.grid(b=True, which='major', color='gray', linestyle='-',alpha=0.6)
ax4.grid(b=True, which='minor', color='gray', linestyle='-',alpha=0.1)
#
#ax4 = fig.add_subplot(2,2,4)
#
#ax4_2=ax4.plot(xiBth9p1J[:,0], xiBth9p1J[:,1], '-', \
#               alpha=0.6, label="Theory:W=9.1J")
#
#ax4_1=ax4.plot(xiBth5p5J[:,0], xiBth5p5J[:,1], '-', \
#               alpha=0.6, label="Theory:W=5.5J")
#
#ax4.errorbar(xiBth9p1Jdat[:,0], xiBth9p1Jdat[:,1], yerr = xiBth9p1Jdat[:,2], fmt='o', \
#                   color = ax4_2[0].get_color(), label="Data:W=9.1J")
#
#ax4.errorbar(xiBth5p5Jdat[:,0], xiBth5p5Jdat[:,1], yerr = xiBth5p5Jdat[:,2], fmt='o', \
#                   color = ax4_1[0].get_color(), label="Data:W=5.5J")
#
#ax4.legend()
#
#ax4.set_xscale('log')
#ax4.set_yscale('log')
#ax4.axes.set_xlim(0.5,250)
#ax4.axes.set_ylim(0.3,12)
#
#
#ax4.axes.set_xlabel(r"time ($2\pi Jt$)", \
#                    fontsize = 12)
#ax4.axes.set_ylabel(r"$\xi_{Bath}$", \
#                    fontsize = 12)
#
#ax4.grid()


#plt.show()
plt.savefig('fig2_py.pdf')
