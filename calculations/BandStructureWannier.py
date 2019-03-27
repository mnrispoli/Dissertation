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

def w1(BF,cp,x,K):
    wout=0
    qs=np.linspace(-1/2,1/2-1/K,K)
    for ii in range(K):
        wout=wout+np.dot(cp[ii]*BF[:,ii],FM(x,K))*np.exp(-1j*x*qs[ii])
    return wout/np.sqrt(2*np.pi*(K-1))

#def w2(BF,cp,x,K):
#    wout=0
#    qs=np.linspace(-1/2,1/2-1/K,K)
#    for ii in range(K):
#        wout=wout+np.dot(cp[ii]*BF[:,ii],FM(x,K))*np.exp(-1j*x*qs[ii])
#    return wout/np.sqrt(2*np.pi*(K-1))

def calJ(K,VoJ,Er,eps,NS,NX):
    Hjq=np.zeros([K,K])
    
    qs=np.linspace(-1/2,1/2-1/K,K)
    
    BS0=np.zeros(K)
    BS1=np.zeros(K)
    BS2=np.zeros(K)
    BS3=np.zeros(K)
    
    dfBF1=np.zeros(K)
    
    BF0=np.zeros([K,K])
    BF1=np.zeros([K,K])
    BF2=np.zeros([K,K])
    
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
        
        BF0[::,qi]=V[::,0]
        BF1[::,qi]=V[::,1]
        BF2[::,qi]=V[::,2]
    
#    plt.plot(qs,BS0)
#    plt.plot(qs,BS1)
#    plt.plot(qs,BS2)
#    plt.plot(qs,BS3)
#    plt.show()
    
    T=int(np.ceil(K/2)-4)
    cp0=np.ones(K)
    cp1=np.ones(K)
    cp2=np.ones(K)
    
    dfBF1[::]=BF1[T+4+2,::]-BF1[T+4+2-1,::]  
    
    for ii in range(K-1):
        if np.sign(BF0[T,ii+1])==np.sign(BF0[T,ii]):
            cp0[ii+1]=cp0[ii]
        else:
            cp0[ii+1]=-cp0[ii]
            
         
        if np.sign(dfBF1[ii+1])==np.sign(dfBF1[ii]):
#        if np.sign(BF1[T+4,ii+1])==np.sign(BF1[T+4,ii]):
            cp1[ii+1]=cp1[ii]
        else:
            cp1[ii+1]=-cp1[ii]
            
        
            
        if np.sign(BF2[T,ii+1])==np.sign(BF2[T,ii]):
            cp2[ii+1]=cp2[ii]
        else:
            cp2[ii+1]=-cp2[ii]
            

    #        if ii == 2:
    #            cp[ii+1]=-cp[ii]
    #        else:
    #            cmpr=np.abs((BF[T,ii+1]-BF[T,ii])/(cp[ii]*BF[T,ii]-cp[ii-1]*BF[T,ii-1]))
    #            if cmpr>0.7 and cmpr<1.3:
    #                cp[ii+1]=cp[ii]
    #            else:
    #                cp[ii+1]=-cp[ii]
    

    
    xx=np.linspace(-NS,NS,NX)*2*np.pi
    dx=xx[1]-xx[0]

    ww0=np.ones(NX)+1j*np.ones(NX)

    ww1=np.ones(NX)+1j*np.ones(NX)

    ww2=np.ones(NX)+1j*np.ones(NX)

   
    for xi in range(NX):
        ww0[xi]=w0(BF0,cp0,xx[xi],K)/np.sqrt(2*np.pi)
        ww1[xi]=w0(BF1,cp1,xx[xi],K)/np.sqrt(2*np.pi)
        ww2[xi]=w0(BF2,cp2,xx[xi],K)/np.sqrt(2*np.pi)
        

    return ww0, ww1, ww2
#   

K=35
VoJs=np.linspace(0.1,10,55)
Jout=np.zeros(len(VoJs))
epsN=np.linspace(0,0.1,55)
Er=1.24
#for VoJi in range(len(VoJs)):
#    Jout[VoJi]=calJ(K,VoJi,Er,eps)
#JoutN=np.zeros([len(epsN),100])
#for jn in range(len(epsN)):
NEr=1

nRange=np.linspace(0.5,NEr,2*NEr)
nRange.reshape([len(nRange),1])

NX=700
NS=7
WOutN0=np.zeros([2*NEr,NX])
WOutN1=np.zeros([2*NEr,NX])
WOutN2=np.zeros([2*NEr,NX])
cnt=0

for nn in nRange:
    WOutN0[cnt,::], WOutN1[cnt,::], WOutN2[cnt,::]=calJ(K,8,Er,0,NS,NX)
    cnt+=1

###
Uo=132

xx=np.linspace(-NS,NS,NX)

fig = plt.figure()
plt.plot(xx,(1-np.cos(xx*2*np.pi))/10)
#plt.plot(xx,WOutN[1,::]/np.sign(WOutN[1,int(NX/2)]))
#plt.plot(xx,WOutN[3,::]/np.sign(WOutN[3,int(NX/2)]))
plt.plot(xx,(WOutN0[1,::]/np.sign(WOutN0[1,int(NX/2)])))
plt.plot(xx,(WOutN1[1,::]/np.sign(WOutN1[1,int(NX/2)])))
plt.plot(xx,(WOutN2[1,::]/np.sign(WOutN2[1,int(NX/2)])))
plt.xlim([-3,3])
#ax1.set_xlim([0,NEr])
#ax1.set_ylim([0,240])
#ax1.grid()
#
#ax2 = fig.add_subplot(122, aspect='auto')
#ax2.semilogy(nRange,JoutN[::,1]*1000)
#ax2.semilogy(nRange,JoutN[::,2]*1000)
#ax2.semilogy(nRange,JoutN[::,3]*1000)
#ax2.semilogy(nRange,JoutN[::,4]*1000)
#ax2.semilogy(nRange,JoutN[::,5]*1000)
#ax2.set_xlim([0,NEr])
#ax2.set_ylim([0,260])
#ax2.grid()

#nRangeNew=nRange.reshape((len(nRange),1))
#DataOut=np.hstack((nRangeNew,JoutN))
#np.savetxt("WannierFx.csv",DataOut, delimiter=",")



