#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 13:23:13 2019

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

def rotx(x,y,theta):
    return x*np.cos(theta)-y*np.sin(theta)

def roty(x,y,theta):
    return y*np.cos(theta)+x*np.sin(theta)

def wz(w,zr,z):
    return w*np.sqrt(1+(z/zr)**2)

def Ef(r,z,w,l,zr):
    return (w/wz(w,zr,z))*np.exp(-(r/wz(w,zr,z))**2)*np.exp(-1j*np.pi*2.0*z/l)

def vLatt(x,y,w,l,zr,theta,Rs):
    xp1=rotx(x,y,theta)
    yp1=roty(x,y,theta)
    
    xp2=rotx(x,y,-theta)
    yp2=roty(x,y,-theta)
    
    return np.abs(\
                  Ef(xp1,yp1,w,l,zr) + Rs*Ef(xp2,yp2,w,l,zr) \
                  )**2
                  
#f = interp2d(X, Y, data, kind='cubic')
#from scipy.interpolate import interp2d
    
xi, yi = np.mgrid[-4:4:0.02,-4:4:0.02]

zi = vLatt(xi.flatten(),yi.flatten(),1,0.6,3,np.pi/4,-0.62)

plt.pcolormesh(xi,yi,zi.reshape(xi.shape), cmap = plt.cm.Greys)

plt.axis('off')
#plt.show()

plt.savefig('AxLattDensity.svg')
