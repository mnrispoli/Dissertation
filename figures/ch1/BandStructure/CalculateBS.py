# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 17:31:04 2019

@author: Matthew
"""

import numpy as np
#from scipy.linalg import expm
import matplotlib.pyplot as plt
import cPickle as pickle

import os
cdir=os.getcwd()
os.chdir('../..')
from plot_colors import *
os.chdir(cdir)

# set some colors for plotting

#matica99A = (0.65,0.00,0.00)
#matica99B = (0.05,0.53,0.63)
#matica99C = (0.44,0.26,0.71)
#matica99D = (0.75,0.36,0.13)
#matica99E = (0.46,0.56,0.01)
#matica99F = (0.66,0.21,0.30)
#matica99G = (0.21,0.39,0.80)
#matica99H = (0.78,0.52,0.04)
#matica99I = (0.52,0.22,0.53)
#matica99J = (0.11,0.56,0.42)
#
#matica97A = (0.37,0.51,0.71)
#matica97B = (0.88,0.61,0.14)
#matica97C = (0.56,0.69,0.19)
#matica97D = (0.92,0.39,0.21)
#matica97E = (0.53,0.47,0.70)
#matica97F = (0.77,0.43,0.10)
#matica97G = (0.36,0.62,0.78)
#matiac97H = (1.00,0.75,0.00)
#matica97I = (0.65,0.38,0.61)
#matica97J = (0.57,0.59,0.00)


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
    return (1/np.sqrt(2*L))*np.exp(-1j*x*np.linspace(-lh,lh,L))

def BFM(x,L,q,Vk):
    #add the bragg scattered  waves + phase shift of iqx plus contributions from
    #hybridized bloch waves found from diagonalized hamiltonian
    return np.dot(Vk[::],FM(x,L)*np.exp(-1j*x*q))

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

def bf(BF,cp,x,K,L):
    #create wannier function from all bloch waves that depend on q
    
    #initialize output
    bfout=np.ones(K)*(0+0*1j)
    
    #q values
    qs=np.linspace(-1.0/2,1.0/2-1.0/K,K)
    
    #add together blcoh waves
    for qi in range(K):
        bfout[qi]=cp[qi]*BFM(x,L,qs[qi],BF[::,qi])/np.sqrt(2*np.pi*(K-1))
    return bfout



def normMe(xvec):
    #normalize vector
    return xvec/np.sqrt(np.dot(xvec,xvec))

def normMeBF(xvec):
    #fix height of bloch function
    return xvec/np.amax(np.real(xvec))


def calBlochFx(L, K, VoJ, Er, eps, NS, NX):

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
    BS4=np.zeros(K)
    BS5=np.zeros(K)
    BS6=np.zeros(K)
    
    BF0=np.zeros([L,K])
    BF1=np.zeros([L,K])
    BF2=np.zeros([L,K])
    BF3=np.zeros([L,K])
    BF4=np.zeros([L,K])
    BF5=np.zeros([L,K])
    BF6=np.zeros([L,K])
    
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
        BS4[qi]=D[4]
        BS5[qi]=D[5]
        BS6[qi]=D[6]
        
        #hybridized bloch waves :: eigenstates (q)
        BF0[::,qi]=V[::,0]
        BF1[::,qi]=V[::,1]
        BF2[::,qi]=V[::,2]
        BF3[::,qi]=V[::,3]
#        BF4[::,qi]=V[::,4]
#        BF5[::,qi]=V[::,5]
#        BF6[::,qi]=V[::,6]
    
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
#    xx=np.linspace(-NS,NS,NX)
#    cntr=0
#    bf0=np.ones((NX,K))+1j*np.ones((NX,K))
#    bf1=np.ones((NX,K))+1j*np.ones((NX,K))
#    bf2=np.ones((NX,K))+1j*np.ones((NX,K))
#    bf3=np.ones((NX,K))+1J*np.ones((NX,K))
#    bf4=np.ones((NX,K))+1J*np.ones((NX,K))
#    bf5=np.ones((NX,K))+1J*np.ones((NX,K))
#    bf6=np.ones((NX,K))+1J*np.ones((NX,K))
    
    #find x values of bloch function with correct phases
#    for x in xx:
#        bf0[cntr,::]=bf(BF0,cp0,x*2*np.pi,K,L)
#        bf1[cntr,::]=bf(BF1,cp1,x*2*np.pi,K,L)
#        bf2[cntr,::]=bf(BF2,cp2,x*2*np.pi,K,L)
#        bf3[cntr,::]=bf(BF3,cp3,x*2*np.pi,K,L)
#        bf4[cntr,::]=bf(BF4,cp3,x*2*np.pi,K,L)
#        bf5[cntr,::]=bf(BF5,cp3,x*2*np.pi,K,L)
#        bf6[cntr,::]=bf(BF6,cp3,x*2*np.pi,K,L)
#        cntr+=1
    
    #normalize height of bloch functions to 1
#    for qi in range(K):
#        bf0[::,qi] = normMeBF(bf0[::,qi])
#        bf1[::,qi] = normMeBF(bf1[::,qi])
#        bf2[::,qi] = normMeBF(bf2[::,qi])
#        bf3[::,qi] = normMeBF(bf3[::,qi])
#        bf4[::,qi] = normMeBF(bf4[::,qi])
#        bf5[::,qi] = normMeBF(bf5[::,qi])
#        bf6[::,qi] = normMeBF(bf6[::,qi])
    
    return BS0, BS1, BS2, BS3, BS4, BS5, BS6, \
             BF0, BF1, BF2, BF3

NS=4.5
NX=1000
L=35
K=32*2 # pick a multiple of 2^n if you want to show periodicity of n*lattice sites nicely
Er=1.24


VoN=np.hstack((np.linspace(0.5,12,24),np.linspace(13,45,(45-13)+1)))
VoN=np.linspace(0.25,1000,4000)

qs=np.linspace(-1.0/2,1.0/2-1.0/K,K)

BS0=np.ones((len(VoN),K))
BS1=np.ones((len(VoN),K))
BS2=np.ones((len(VoN),K))
BS3=np.ones((len(VoN),K))
BS4=np.ones((len(VoN),K))
BS5=np.ones((len(VoN),K))
BS6=np.ones((len(VoN),K))

BF0=np.ones((len(VoN),L,K))*1j
BF1=np.ones((len(VoN),L,K))*1j
BF2=np.ones((len(VoN),L,K))*1j
BF3=np.ones((len(VoN),L,K))*1j


cntr=0
for nn in VoN:
    print(nn)
    BS0[cntr,::], BS1[cntr,::], BS2[cntr,::], BS3[cntr,::], BS4[cntr,::], BS5[cntr,::], BS6[cntr,::], \
        BF0[cntr,::,::], BF1[cntr,::,::], BF2[cntr,::,::], BF3[cntr,::,::] = calBlochFx(L, K, nn, Er, 0, NS, NX)
    cntr+=1
    
xx=np.linspace(-NS,NS,NX)
xxout=xx.reshape((len(xx),1))

#plt.plot(xx,4*(1-np.cos(xx*np.pi)**2)-0.5,'k-')

#Vout=(1-np.cos(xx*np.pi)**2).reshape((len(xx),1))
#plt.plot(xx,np.real(bf0[::,K/2])/3,'-', color = maticaB)
#plt.plot(xx,np.imag(bf0[::,K/2])/3,'--', color = maticaB)
#
#plt.plot(xx,np.real(bf0[::,K/2+K/8+1])/3+1,'-', color = maticaC)
#plt.plot(xx,np.imag(bf0[::,K/2+K/8+1])/3+1,'--', color = maticaC)
#
#plt.plot(xx,np.real(bf0[::,K/2+K/4+1])/3+2,'-', color = maticaA)
#plt.plot(xx,np.imag(bf0[::,K/2+K/4+1])/3+2,'--', color = maticaA)
#
#plt.plot(xx,np.real(bf0[::,0])/3+3,'-', color = maticaD)
#plt.plot(xx,np.imag(bf0[::,0])/3+3,'--', color = maticaD)
#plt.show()

def FPD(qs,l):
    #free particle dispersion curve offset by "l"
    return 4*(qs-l)**2

plt.plot(qs, BS0[1,::], '-', color = matica97A)
plt.plot(qs, FPD(qs,0), '--', color = matica97A)
#plt.plot(qs, FPD(qs,-1), '--', color = matica97A)

plt.plot(qs, BS1[1,::], '-', color = matica97B)
plt.plot(qs, FPD(qs,1), '--', color = matica97B)
plt.plot(qs, FPD(qs,-1), '--', color = matica97B)

plt.plot(qs, BS2[2,::], '-', color = matica97C)
plt.plot(qs, FPD(qs,2), '--', color = matica97C)
plt.plot(qs, FPD(qs,-2), '--', color = matica97C)

plt.plot(qs, BS3[3,::], '-', color = matica97D)

#plt.plot(qs, FPD(qs,1), '--', color = matica97A)
#plt.plot(qs, FPD(qs,-1), '--', color = matica97A)

#plt.plot(qs,BS1-np.amin(BS0),color=matica97B,'-')
#plt.plot(qs,BS1-np.amin(BS0),color=matica97B,'--')
#
#plt.plot(qs,BS2-np.amin(BS0),color=matica97C,'-')
#plt.plot(qs,BS3-np.amin(BS0),color=matica97D,'-')




#
#W0out=wn0.reshape((len(xx),1))
#W1out=wn1.reshape((len(xx),1))
#W2out=wn2.reshape((len(xx),1))
#W3out=wn3.reshape((len(xx),1))
#
#plt.plot(xx,9*wn0+0)
#plt.plot(xx,3*wn1+1)
#plt.plot(xx,3*wn2+2)
#plt.plot(xx,3*wn3+3)
#

np.savetxt("BSLattDepth.csv",VoN,delimiter=",")

np.savetxt("BS0.csv",BS0,delimiter=",")
np.savetxt("BS1.csv",BS1,delimiter=",")
np.savetxt("BS2.csv",BS2,delimiter=",")
np.savetxt("BS3.csv",BS3,delimiter=",")
np.savetxt("BS4.csv",BS4,delimiter=",")
np.savetxt("BS5.csv",BS5,delimiter=",")
np.savetxt("BS6.csv",BS6,delimiter=",")

#np.savetxt("BF0.csv",BF0,delimiter=",")
#np.savetxt("BF1.csv",BF1,delimiter=",")
#np.savetxt("BF2.csv",BF2,delimiter=",")
#np.savetxt("BF3.csv",BF3,delimiter=",")
pickle.dump( BF0, open( "BF0.pkl", "wb" ) )
pickle.dump( BF1, open( "BF1.pkl", "wb" ) )
pickle.dump( BF2, open( "BF2.pkl", "wb" ) )
pickle.dump( BF3, open( "BF3.pkl", "wb" ) )

BF0L = pickle.load( open( "BF0.pkl", "rb" ) )
np.savetxt("BSQs.csv",qs,delimiter=",")
