# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 17:31:04 2019

@author: Matthew
"""

import numpy as np
#from scipy.linalg import expm
import matplotlib.pyplot as plt

import os
cdir=os.getcwd()
os.chdir('../..')
from plot_colors import *
os.chdir(cdir)

def HamJQ(q,VoJ,L):
    #create optical lattice single-particle hamiltonian
    lh=(L-1)/2
    
    dvec=4*np.square(q+np.linspace(-lh,lh,L))
    
    ovec=np.ones(L-1)*(-VoJ/4)
    
    Hjq=np.diag(dvec)+np.diag(ovec,1)+np.diag(ovec,-1)
    return Hjq

def FM(x,L):
    #fourier modes from bragg scattered waves
    lh=(L-1)/2
    return np.exp(-1j*x*np.linspace(-lh,lh,L))

def BFM(x,L,q,Vk):
    #add the bragg scattered  waves + phase shift of iqx plus contributions from
    #hybridized bloch waves found from diagonalized hamiltonian
    return np.dot(Vk[::],FM(x,L)*np.exp(-1j*x*q))/np.sqrt(2*np.pi*((L)/2))

def w0(BF,cp,x,K,L):
    #create wannier function from all bloch waves that depend on q
    
    #initialize output
    wout=0+0*1j
    
    #q values
    qs=np.linspace(-1.0/2,1.0/2-1.0/K,K)
    
    #add together blcoh waves
    for qi in range(K):
        wout=wout+cp[qi]*BFM(x,L,qs[qi],BF[::,qi])/np.sqrt(2*np.pi*(K-1))
    return wout


def normMe(xvec):
    #integrate vector
    return xvec/np.sqrt(np.dot(np.conj(xvec),xvec))


def calWannier(L, K, VoJ, Er, eps, NS, NX):

    #all quasi-momentums states that will be considered
    qs=np.linspace(-1.0/2,1.0/2-1.0/K,K)
    
    
    # initialize a ton of stuff for first 4 bands in lattice
    cp0=np.ones(K)+1j*np.ones(K)
    cp1=np.ones(K)+1j*np.ones(K)
    cp2=np.ones(K)+1j*np.ones(K)
    cp3=np.ones(K)+1j*np.ones(K)
    
    BS0=np.zeros(K)
    BS1=np.zeros(K)
    BS2=np.zeros(K)
    BS3=np.zeros(K)
    
    BF0=np.zeros([L,K])
    BF1=np.zeros([L,K])
    BF2=np.zeros([L,K])
    BF3=np.zeros([L,K])
    
    theta0=np.zeros(K)
    phiX0=np.ones(K)+1j*np.ones(K)
    
    theta1=np.zeros(K)
    phiX1=np.ones(K)+1j*np.ones(K)
    
    theta2=np.zeros(K)
    phiX2=np.ones(K)+1j*np.ones(K)
    
    theta3=np.zeros(K)
    phiX3=np.ones(K)+1j*np.ones(K)
    
    #test x value to get correct phases for bloch waves
    xtE=0.0001
    
    for qi in range(K):
        
        q=qs[qi]
        
        #create single-particle hamiltonian in lattice
        H=HamJQ(q,VoJ,L)       
        D, V = np.linalg.eig(H)
        
        #sort diagonalied hamiltonian
        idx = D.argsort()
        D = D[idx]
        V = V[::,idx]
        
        #band strucutre :: eigenenergies (q)
        BS0[qi]=D[0]
        BS1[qi]=D[1]
        BS2[qi]=D[2]
        BS3[qi]=D[3]
        
        #hybridized bloch waves :: eigenstates (q)
        BF0[::,qi]=V[::,0]
        BF1[::,qi]=V[::,1]
        BF2[::,qi]=V[::,2]
        BF3[::,qi]=V[::,3]
    
        #get correct phases for ground state, depends on x=0
        phiX0[qi]=BFM(xtE,L,qs[qi],BF0[:,qi])
        theta0[qi]=np.angle(phiX0[qi])
        cp0[qi]=np.exp(-1j*theta0[qi])
        
        #get correct phases for 1st excited state, depends on d/dx @ x=0
        phiX1[qi]= (BFM(xtE,L,qs[qi],BF1[:,qi]) - \
                    BFM(-xtE,L,qs[qi],BF1[:,qi]))/(2*xtE)
        theta1[qi]=np.angle(phiX1[qi])
        cp1[qi]=np.exp(-1j*theta1[qi])
        
        #get correct phases for 2nd excited state, depends on x=0 value
        phiX2[qi]= BFM(xtE,L,qs[qi],BF2[:,qi])
        theta2[qi]=np.angle(phiX2[qi])
        cp2[qi]=np.exp(-1j*theta2[qi])
        
        #get correct phases for 3rd excited state, depends on d/dx @ x=0
        phiX3[qi]= (BFM(xtE,L,qs[qi],BF3[:,qi]) - \
                    BFM(-xtE,L,qs[qi],BF3[:,qi]))/(2*xtE)
        theta3[qi]=np.angle(phiX3[qi])
        cp3[qi]=np.exp(-1j*theta3[qi])
    
    #initialize section
    xx=np.linspace(-NS,NS,NX)
    dx=xx[1]-xx[0]
    cntr=0
    wn0=np.ones(NX)+1j*np.ones(NX)
    wn1=np.ones(NX)+1j*np.ones(NX)
    wn2=np.ones(NX)+1j*np.ones(NX)
    wn3=np.ones(NX)+1J*np.ones(NX)
    
    #find x values of Wannier function with correct phases
    for x in xx:
        wn0[cntr]=w0(BF0,cp0,x*2*np.pi,K,L)
        wn1[cntr]=w0(BF1,cp1,x*2*np.pi,K,L)
        wn2[cntr]=w0(BF2,cp2,x*2*np.pi,K,L)
        wn3[cntr]=w0(BF3,cp3,x*2*np.pi,K,L)
        cntr+=1
    
    #normalize wannier Functions, I'm missing some factor that depends on L or K above
    wn0 = normMe(wn0)
    wn1 = normMe(wn1)
    wn2 = normMe(wn2)
    wn3 = normMe(wn3)
    
    return wn0, wn1, wn2, wn3, BF0

NS=2.5
NX=1000

K=35
L=35

#VoN=np.linspace(0.5,100,1)
#
#wn0=np.ones(len(VoN),K)
#wn1=np.ones(len(VoN),K)
#wn2=np.ones(len(VoN),K)
#wn3=np.ones(len(VoN),K)
#
#cntr=0

wn0, wn1, wn2, wn3, bf0= calWannier(L, K, 1, 1.24, 0, NS, NX)

xx=np.linspace(-NS,NS,NX)
xxout=xx.reshape((len(xx),1))

plt.plot(xx,4*(1-np.cos(xx*np.pi)**2)-0.5,'k-')

Vout=(1-np.cos(xx*np.pi)**2).reshape((len(xx),1))
W0out=wn0.reshape((len(xx),1))
W1out=wn1.reshape((len(xx),1))
W2out=wn2.reshape((len(xx),1))
W3out=wn3.reshape((len(xx),1))

plt.plot(xx,9*wn0+0)
plt.plot(xx,3*wn1+1)
plt.plot(xx,3*wn2+2)
plt.plot(xx,3*wn3+3)

np.savetxt("LattOut.csv",np.hstack((xxout,Vout)),delimiter=",")
np.savetxt("W0_1Er.csv",np.hstack((xxout,W0out)),delimiter=",")
np.savetxt("W1_1Er.csv",np.hstack((xxout,W1out)),delimiter=",")
np.savetxt("W2_1Er.csv",np.hstack((xxout,W2out)),delimiter=",")
np.savetxt("W3_1Er.csv",np.hstack((xxout,W3out)),delimiter=",")
