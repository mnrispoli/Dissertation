# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 17:31:04 2019

@author: Matthew
"""

import numpy as np
#from scipy.linalg import expm
import matplotlib.pyplot as plt


def HamJQ(q,VoJ,L):
    lh=(L-1)/2
    
    dvec=4*np.square(q+np.linspace(-lh,lh,L))
    
    ovec=np.ones(L-1)*(-VoJ/4)
    
    Hjq=np.diag(dvec)+np.diag(ovec,1)+np.diag(ovec,-1)
    return Hjq

def FM(x,L):
    lh=(L-1)/2
    return (1/np.sqrt(2*L))*np.exp(-1j*x*np.linspace(-lh,lh,L))

def BFM(x,L,q,Vk):
    return np.dot(Vk[::],FM(x,L)*np.exp(-1j*x*q))

def w0(BF,cp,x,K,L):
    wout=0+0*1j
    
    qs=np.linspace(-1.0/2,1.0/2-1.0/K,K)
    
    for qi in range(K):
        wout=wout+cp[qi]*BFM(x,L,qs[qi],BF[::,qi])/np.sqrt(2*np.pi*(K-1))
    return wout

def normMe(xvec):
    return xvec/np.sqrt(np.dot(xvec,xvec))


def calWannier(L, K, VoJ, Er, eps, NS, NX):
    
    
    Hjq=np.zeros([L,L])
    
    qs=np.linspace(-1.0/2,1.0/2-1.0/K,K)
    
    cp0=np.ones(K)+1j*np.ones(K)
    cp1=np.ones(K)+1j*np.ones(K)
    cp2=np.ones(K)+1j*np.ones(K)
    cp3=np.ones(K)+1j*np.ones(K)
    
    BS0=np.zeros(K)
    BS1=np.zeros(K)
    BS2=np.zeros(K)
    BS3=np.zeros(K)
    
    dfBF1=np.zeros(K)
    
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
    
    xtE=0.0001
    
    for qi in range(K):
        
        q=qs[qi]
        
        H=HamJQ(q,VoJ,L)       
        D, V = np.linalg.eig(H)
        
        idx = D.argsort()
        D = D[idx]
        V = V[::,idx]
        
        BS0[qi]=D[0]
        BS1[qi]=D[1]
        BS2[qi]=D[2]
        BS3[qi]=D[3]
        
        BF0[::,qi]=V[::,0]
        BF1[::,qi]=V[::,1]
        BF2[::,qi]=V[::,2]
        BF3[::,qi]=V[::,3]
    
        phiX0[qi]=BFM(xtE,L,qs[qi],BF0[:,qi])
        theta0[qi]=np.angle(phiX0[qi])
        cp0[qi]=np.exp(-1j*theta0[qi])
        
        phiX1[qi]= (BFM(xtE,L,qs[qi],BF1[:,qi]) - \
                    BFM(-xtE,L,qs[qi],BF1[:,qi]))/(2*xtE)
        theta1[qi]=np.angle(phiX1[qi])
        cp1[qi]=np.exp(-1j*theta1[qi])
        
        phiX2[qi]= BFM(xtE,L,qs[qi],BF2[:,qi])
        theta2[qi]=np.angle(phiX2[qi])
        cp2[qi]=np.exp(-1j*theta2[qi])
        
        phiX3[qi]= (BFM(xtE,L,qs[qi],BF3[:,qi]) - \
                    BFM(-xtE,L,qs[qi],BF3[:,qi]))/(2*xtE)
        theta3[qi]=np.angle(phiX3[qi])
        cp3[qi]=np.exp(-1j*theta3[qi])
    
    xx=np.linspace(-NS,NS,NX)
    cntr=0
    wn0=np.ones(NX)+1j*np.ones(NX)
    wn1=np.ones(NX)+1j*np.ones(NX)
    wn2=np.ones(NX)+1j*np.ones(NX)
    wn3=np.ones(NX)+1J*np.ones(NX)
    
    for x in xx:
        wn0[cntr]=w0(BF0,cp0,x*2*np.pi,K,L)
        wn1[cntr]=w0(BF1,cp1,x*2*np.pi,K,L)
        wn2[cntr]=w0(BF2,cp2,x*2*np.pi,K,L)
        wn3[cntr]=w0(BF3,cp3,x*2*np.pi,K,L)
        cntr+=1
    
    wn0 = normMe(wn0)
    wn1 = normMe(wn1)
    wn2 = normMe(wn2)
    wn3 = normMe(wn3)
    
    return wn0, wn1, wn2, wn3

NS=2.5
NX=1000

wn0, wn1, wn2, wn3= calWannier(35, 55, 8, 1.24, 0, NS, NX)

xx=np.linspace(-NS,NS,NX)
plt.plot(xx,4*(1-np.cos(xx*np.pi)**2)-0.5,'k-')

plt.plot(xx,3*wn0+0)
plt.plot(xx,3*wn1+1)
plt.plot(xx,3*wn2+2)
plt.plot(xx,3*wn3+3)

#    plt.plot(qs,BS0)
#    plt.plot(qs,BS1)
#    plt.plot(qs,BS2)
#    plt.plot(qs,BS3)
#    plt.show()
#
#T=int(np.ceil(K/2)-4)
#
#cp0=np.ones(K)
#cp1=np.ones(K)
#cp2=np.ones(K)
#
#
#
#
#dfBF1[::]=BF1[T+4,::]-BF1[T+4-1,::]  
#
##for ii in range(K-1):
##    if np.sign(BF0[T,ii+1])==np.sign(BF0[T,ii]):
##        cp0[ii+1]=cp0[ii]
##    else:
##        cp0[ii+1]=-cp0[ii]
##        
##     
###        if np.sign(dfBF1[ii+1])==np.sign(dfBF1[ii]):
##    if np.sign(BF1[T+6,ii+1])==np.sign(BF1[T+6,ii]):
##        cp1[ii+1]=cp1[ii]
##    else:
##        cp1[ii+1]=-cp1[ii]
##        
##    
##        
##    if np.sign(BF2[T,ii+1])==np.sign(BF2[T,ii]):
##        cp2[ii+1]=cp2[ii]
##    else:
##        cp2[ii+1]=-cp2[ii]
#        
#
##        if ii == 2:
##            cp[ii+1]=-cp[ii]
##        else:
##            cmpr=np.abs((BF[T,ii+1]-BF[T,ii])/(cp[ii]*BF[T,ii]-cp[ii-1]*BF[T,ii-1]))
##            if cmpr>0.7 and cmpr<1.3:
##                cp[ii+1]=cp[ii]
##            else:
##                cp[ii+1]=-cp[ii]
#
#
#
#xx=np.linspace(-NS,NS,NX)*2*np.pi
#dx=xx[1]-xx[0]
#
#ww0=np.ones(NX)+1j*np.ones(NX)
#
#ww1=np.ones(NX)+1j*np.ones(NX)
#
#ww2=np.ones(NX)+1j*np.ones(NX)
#
#   
#for xi in range(NX):
#    ww0[xi]=w0(BF0,cp0,xx[xi],K)/np.sqrt(2*np.pi)
#    ww1[xi]=w0(BF1,cp1,xx[xi],K)/np.sqrt(2*np.pi)
#    ww2[xi]=w0(BF2,cp2,xx[xi],K)/np.sqrt(2*np.pi)
#        
