#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  2 18:23:03 2019

@author: matthewrispoli
"""

import numpy as np
#from scipy.linalg import expm
import matplotlib.pyplot as plt
#import matplotlib.image as mpimg
#import time

#make hamiltonian matrix 
def HamJQ(q,VoJ,K):
    kh=(K-1)/2
    dvec=4*np.square(q+np.linspace(-kh,kh,K))
    ovec=np.ones(K-1)*(-VoJ/4)
    Hjq=np.diag(dvec)+np.diag(ovec,1)+np.diag(ovec,-1)
    return Hjq

def FM(x,K):
    L=K-1
    return (1/np.sqrt(L))*np.exp(-1j*x*np.linspace(-L/2,L/2,K))

def w0(BF,cp,x,K):
    wout=0
    qs=np.linspace(-1/2,1/2-1/K,K)
    for ii in range(K):
        wout=wout+np.dot(cp[ii]*BF[:,ii],FM(x,K))*np.exp(-1j*x*qs[ii])
    return wout/np.sqrt(2*np.pi*(K-1))

def calJ(K,VoJ,Er,eps):
    Hjq=np.zeros([K,K])
    
    qs=np.linspace(-1/2,1/2-1/K,K)
    
    BS0=np.zeros(K)
    BS1=np.zeros(K)
    BS2=np.zeros(K)
    BS3=np.zeros(K)
    
    BF=np.zeros([K,K])
    
    for qi in range(K):
        q=qs[qi]
        H=HamJQ(q,VoJ,K)       
        D, V = np.linalg.eig(H)
        
        idx = D.argsort()
        D = D[idx]
        V = V[::,idx]
        BS0[qi]=D[0]
        BS1[qi]=D[1]
        BS2[qi]=D[2]
        BS3[qi]=D[3]
        
        BF[::,qi]=V[::,0]
    
#    plt.plot(qs,BS0)
#    plt.plot(qs,BS1)
#    plt.plot(qs,BS2)
#    plt.plot(qs,BS3)
#    plt.show()
    
    T=int(np.ceil(K/2)-4)
    cp=np.ones(K)
    for ii in range(K-1):
        if np.sign(BF[T,ii+1])==np.sign(BF[T,ii]):
            cp[ii+1]=cp[ii]
        else:
            cp[ii+1]=-cp[ii]
    #        if ii == 2:
    #            cp[ii+1]=-cp[ii]
    #        else:
    #            cmpr=np.abs((BF[T,ii+1]-BF[T,ii])/(cp[ii]*BF[T,ii]-cp[ii-1]*BF[T,ii-1]))
    #            if cmpr>0.7 and cmpr<1.3:
    #                cp[ii+1]=cp[ii]
    #            else:
    #                cp[ii+1]=-cp[ii]
    
    NX=300
    NS=4
    xx=np.linspace(-NS,NS,NX)*2*np.pi
    dx=xx[1]-xx[0]
    yy=np.ones(NX)
    ww=np.ones(NX)+1j*np.ones(NX)
    yy2=np.ones(NX)
    ww2=np.ones(NX)+1j*np.ones(NX)
    yy3=np.ones(NX)
    ww3=np.ones(NX)+1j*np.ones(NX)
    
    for xi in range(NX):
        ww[xi]=w0(BF,cp,xx[xi],K)/np.sqrt(2*np.pi)
        yy[xi]=np.square(np.abs(ww[xi]))
        ww2[xi]=w0(BF,cp,xx[xi]+(1)*2*np.pi,K)/np.sqrt(2*np.pi)
        yy2[xi]=np.square(np.abs(ww2[xi]))
        ww3[xi]=w0(BF,cp,xx[xi]+(1+1)*2*np.pi,K)/np.sqrt(2*np.pi)
        yy3[xi]=np.square(np.abs(ww3[xi]))
    
    dw=(ww[1:len(yy)]-ww[0:-1])/dx
    ddw=(dw[1:len(dy)]-dw[0:-1])/dx
    
    gr=(1+np.sqrt(5))/2-1
    hv=-VoJ*np.cos(xx[1:-1])/2
    hvd1=-VoJ*eps*(np.cos(xx[1:-1]*gr+np.random.rand()*2*np.pi))/2
    
#    hvd2=-VoJ*eps*(np.random.rand()*np.cos(2*xx[1:-1]/6))/2
#    hvd3=-VoJ*eps*(np.random.rand()*np.cos(3*xx[1:-1]/6))/2
#    hvd4=-VoJ*eps*(np.random.rand()*np.cos(4*xx[1:-1]/6))/2
#    hvd5=-VoJ*eps*(np.random.rand()*np.cos(5*xx[1:-1]/6))/2 
                             
    hp=-4*ddw
    
    J=Er*abs(sum(np.conj(ww2[1:-1])*hp) \
             +sum(np.conj(ww2[1:-1])*hvd1*ww[1:-1])\
             +sum(np.conj(ww2[1:-1])*hv*ww[1:-1]))
    W=Er*abs(sum(np.conj(ww[1:-1]*hp))\
             +sum(np.conj(ww[1:-1])*hvd1*ww[1:-1])\
             +sum(np.conj(ww[1:-1])*hv*ww[1:-1]))
    
    W=Er*abs(sum(np.conj(ww[1:-1]*hp))\
         +sum(np.conj(ww[1:-1])*hvd1*ww[1:-1])\
         +sum(np.conj(ww[1:-1])*hv*ww[1:-1]))
    tmp=np.zeros(2)
    tmp[0]=J;
    tmp[1]=W;
#             +sum(np.conj(ww2[1:-1])*hvd1*ww[1:-1])\
#             +sum(np.conj(ww2[1:-1])*hvd2*ww[1:-1])\
#             +sum(np.conj(ww2[1:-1])*hvd3*ww[1:-1])\
#             +sum(np.conj(ww2[1:-1])*hvd4*ww[1:-1])\
#             +sum(np.conj(ww2[1:-1])*hvd5*ww[1:-1]))
    Jnn=Er*abs(sum(np.conj(ww3[1:-1])*hp)+sum(np.conj(ww3[1:-1])*hv*ww[1:-1]))
    return tmp
#    plt.plot(xx/(2*np.pi),yy)
#    plt.plot(xx/(2*np.pi),yy2)
#    plt.grid()
#    plt.show()

K=35
VoJs=np.linspace(0.1,10,55)
Jout=np.zeros(len(VoJs))
epsN=np.linspace(0,0.1,55)
Er=1.24
#for VoJi in range(len(VoJs)):
#    Jout[VoJi]=calJ(K,VoJi,Er,eps)
#JoutN=np.zeros([len(epsN),100])
#for jn in range(len(epsN)):
JoutN=np.zeros([10,2])
for nn in range(10):
    JoutN[nn,::]=calJ(K,8,Er,0.08)

