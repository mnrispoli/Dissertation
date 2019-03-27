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

atomGreen = '#74c476'
intRed = '#de2d26'
thermalOrange = '#fd8d3c'
criticalPurple = '#756bb1'
mblBlue = '#2b8cbe'

style="Simple,tail_width=2,head_width=8,head_length=12"
kw = dict(arrowstyle=style, color="gray")



cm = 1/2.54

#point offset for clarity
dx=0.05
titleFS = 10
axisFS = 8
annFS = 5

dw_atom=0.15
dw_shdw=0.2

#Interaction Strength
Uint=0.2

#Ground States and interaction shifts
E1=0.3
LG1 = 0.13

E2=0.3+Uint
LG2 = 0.19

E3=0.3+3*Uint
LG3 = 0.34

#U offsets in x,y
duangle=-0.05*2*np.pi
dur=0.06
dux=dur*np.cos(duangle)
duy=dur*np.sin(duangle)

#Number of sites to be plotted
NS=3

#Annotation arrow parameters
AWL=0.1
AWH=0.05

#Atom locations parameters
atomLocs = np.linspace(0.5,NS-0.5,NS)
a_latt=(atomLocs[1]-atomLocs[0])

#Set several plot parameters
setPltPars(rcpars_crtn)


#   atom itself

def atomSolid(atom_xcenter,atom_ycenter, atom_width, atom_height, \
              atom_angle):
    #xcenter, ycenter = 0.38*cm, 0.52*cm
    
#    atom_xcenter, atom_ycenter = 1.5, 0.5
    
    #width, height = 1e-1*cm, 3e-1*cm
    
#    atom_width, atom_height = 0.15, 0.15
    
    atom_angle = -0
    
    atom_theta = np.deg2rad(np.arange(0.0, 360.0, 1.0))
    atom_x = 0.5 * atom_width * np.cos(atom_theta)
    atom_y = 0.5 * atom_height * np.sin(atom_theta)
    
    atom_rtheta = np.radians(atom_angle)
    
    atom_R = np.array([
        [np.cos(atom_rtheta), -np.sin(atom_rtheta)],
        [np.sin(atom_rtheta),  np.cos(atom_rtheta)],
        ])
    
    atom_x, atom_y = np.dot(atom_R, np.array([atom_x, atom_y]))
    atom_x += atom_xcenter
    atom_y += atom_ycenter
    
    return atom_x, atom_y

#    Shadow behind atom
def atomShadow(atom_shdw_xcenter, atom_shdw_ycenter, atom_shdw_width, \
               atom_shdw_height, atom_shdw_angle):

    atom_shdw_theta = np.deg2rad(np.arange(0.0, 360.0, 1.0))
    atom_shdw_x = 0.5 * atom_shdw_width * np.cos(atom_shdw_theta)
    atom_shdw_y = 0.5 * atom_shdw_height * np.sin(atom_shdw_theta)
    
    atom_shdw_rtheta = np.radians(atom_shdw_angle)
    
    atom_shdw_R = np.array([
        [np.cos(atom_shdw_rtheta), -np.sin(atom_shdw_rtheta)],
        [np.sin(atom_shdw_rtheta),  np.cos(atom_shdw_rtheta)],
        ])
    
    
    atom_shdw_x, atom_shdw_y = np.dot(atom_shdw_R, np.array([atom_shdw_x, atom_shdw_y]))
    atom_shdw_x += atom_shdw_xcenter
    atom_shdw_y += atom_shdw_ycenter
    
    return atom_shdw_x, atom_shdw_y




#start making plots!!!        
                
fig = plt.figure()

ax = fig.add_subplot(111, aspect='equal')

# add lattice

xx=np.linspace(atomLocs[0]-a_latt/2,atomLocs[-1]+a_latt/2,1000);

yy=(np.cos(xx*2*np.pi)+1)/2
ax.plot(xx,yy, linewidth = 2.5, color='k', zorder = 1)
ax.axis('off')


atom_shdw_x, atom_shdw_y = atomShadow(atomLocs[1]+dux,E2+duy,dw_shdw,dw_shdw,0)
ax.fill(atom_shdw_x, atom_shdw_y, alpha=1, facecolor = 'white',
        edgecolor='white', linewidth=0, zorder=10)

atom_x, atom_y = atomSolid(atomLocs[1]+dux,E2+duy,dw_atom,dw_atom,0)
ax.fill(atom_x, atom_y, alpha=0.8, facecolor = atomGreen, \
        edgecolor=atomGreen, linewidth=0, zorder=12)



atom_shdw_x, atom_shdw_y = atomShadow(atomLocs[1]-dux,E2-duy,dw_shdw,dw_shdw,0)
ax.fill(atom_shdw_x, atom_shdw_y, alpha=1, facecolor = 'white',
        edgecolor='white', linewidth=0, zorder=13)

atom_x, atom_y = atomSolid(atomLocs[1]-dux,E2-duy,dw_atom,dw_atom,0)
ax.fill(atom_x, atom_y, alpha=0.8, facecolor = atomGreen, \
        edgecolor=atomGreen, linewidth=0, zorder=14)



#e1 = patches.Ellipse((xcenter, ycenter), width, height,
#                     angle=angle, linewidth=1, fill=False, zorder=2)

#ax.add_patch(e1)



#e1.set_hatch('/')
#e1.set_linestyle('-.')


def addGSL(Elvl, Loc, LG, axIn):
# add lines for ground states
    xb=np.linspace(Loc-LG,Loc+LG,100);
    yb=np.ones(100)*Elvl;
    axIn.plot(xb,yb, linewidth = 3, alpha = 0.5, color='gray', zorder = 1)
 
addGSL(E1,atomLocs[0],LG1,ax)    
addGSL(E1,atomLocs[1],LG1,ax)    
addGSL(E1,atomLocs[2],LG1,ax)    

addGSL(E2,atomLocs[0],LG2,ax)   
addGSL(E2,atomLocs[1],LG2,ax)   
addGSL(E2,atomLocs[2],LG2,ax)   

#addGSL(E2,atomLocs[0],LG2,ax)   
#addGSL(E3,atomLocs[0],LG3,ax)   

# add annotations on top
#a2 = patches.FancyArrowPatch((atomLocs[1]-AWL,E2+AWH), (atomLocs[0]+AWL,E2+AWH), zorder = 2, \
#                             connectionstyle="arc3,rad=0.4", **kw)
#ax.add_patch(a2)

DUH=0.4

a2 = patches.FancyArrowPatch((atomLocs[1]/2+atomLocs[2]/2,E1), (atomLocs[1]/2+atomLocs[2]/2,E2), zorder = 2, \
                             connectionstyle="arc3,rad=0", **kw)
ax.add_patch(a2)

#a3 = patches.FancyArrowPatch((atomLocs[1]/2+atomLocs[2]/2,E2), (atomLocs[1]/2+atomLocs[2]/2,E1), zorder = 2, \
#                             connectionstyle="arc3,rad=0", **kw)
#ax.add_patch(a3)

#a3 = patches.FancyArrowPatch((atomLocs[1]+AWL,E2+AWH), (atomLocs[2]-AWL,E2+AWH), zorder = 2, \
#                             connectionstyle="arc3,rad=-0.4", **kw)
#ax.add_patch(a3)
#

#plt.rc('text', usetex=True)
#plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

ax.annotate(r"U", xytext=(atomLocs[1]/2+atomLocs[2]/2, E1/2+E2/2), xy=(0,0) \
            )


plt.savefig('U_n2.pdf')
#plt.show()