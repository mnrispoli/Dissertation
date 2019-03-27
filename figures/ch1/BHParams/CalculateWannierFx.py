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
    lh=(L)/2
    return np.exp(-1j*x*2*np.pi*np.linspace(-lh,lh,L))/np.sqrt(L)

def BFM(x,L,q,Vk):
    #add the bragg scattered  waves + phase shift of iqx plus contributions from
    #hybridized bloch waves found from diagonalized hamiltonian
    return np.dot(Vk[::],FM(x,L)*np.exp(-1j*x*2*np.pi*q))/np.sqrt(2*np.pi*L)

def w0(BF,cp,x,K,L):
    #create wannier function from all bloch waves that depend on q
    
    #initialize output
    wout=0+0*1j
    
    #q values
    qs=np.linspace(-1.0/2,1.0/2-1.0/K,K)
    
    #add together blcoh waves
    for qi in range(K):
        wout=wout+cp[qi]*BFM(x,L,qs[qi],BF[::,qi])
    return wout


def normMe(xvec):
    #integrate vector
    return xvec/np.sqrt(np.dot(np.conj(xvec),xvec))

def cPhaseEven(xval,L,qv,BF):
        #get correct phases for even bands @ x=xval
        phiX=BFM(xval,L,qv,BF)
        theta=np.angle(phiX)
        cp=np.exp(-1j*theta)
        return phiX, theta, cp
    
def cPhaseOdd(xval,L,qv,BF):
        #get correct phases for odd bands, depends on d/dx @ x=xval
        dx=0.001
        
        phiX=(BFM(xval+dx,L,qv,BF) - \
              BFM(xval-dx,L,qv,BF))/(2*dx)
        theta=np.angle(phiX)
        cp=np.exp(-1j*theta)
        return phiX, theta, cp
    
def D2(xvec,xx):
    NX=len(xx)
    dx=xx[1]-xx[0]
    
    d1=(xvec[1:NX]-xvec[0:(NX-1)])
    
    d2=(d1[1:(NX-1)]-d1[0:(NX-2)])/dx
    

    
    return np.append( \
                     np.append(np.array([0]), \
                     np.array([d2])), \
                     np.array([0]))

def calBHJ(psi1,psi2,Vext,xx,xinds):
    dx=xx[1]-xx[0]
    
    psi1=psi1[xinds]
    psi2=psi2[xinds]
    Vext=Vext[xinds]
    xx=xx[xinds]
    
    return -(-4*np.dot(np.conj(psi1),D2(psi2,xx))+np.dot(np.conj(psi1),Vext*psi2))*dx

def calBHU(psi,xx,xinds):
    dx=xx[1]-xx[0]
    psi=psi[xinds]
    xx=xx[xinds]
    return (np.dot(np.conj(psi)*np.conj(psi),psi*psi))*dx

def calBHE(psi,Vext,xx,xinds):
    dx=xx[1]-xx[0]
    psi=psi[xinds]
    Vext=Vext[xinds]
    xx=xx[xinds]
    return (-4*np.dot(np.conj(psi),D2(psi,xx))+np.dot(np.conj(psi),Vext*psi))*dx

def calWannier(L, K, VoJ, Er, eps, NS, NX):

    #all quasi-momentums states that will be considered
    qs=np.linspace(-1.0/2,1.0/2-1.0/K,K)
    
    
    # initialize a ton of stuff for first 4 bands in lattice   
    BS0=np.zeros(K)
    BS1=np.zeros(K)
    BS2=np.zeros(K)
    BS3=np.zeros(K)
    
    BF0=np.zeros([L,K])
    BF1=np.zeros([L,K])
    BF2=np.zeros([L,K])
    BF3=np.zeros([L,K])
    
    #things needed for making wannier fx
    #someone should just fix this to be a matrix
    cp0N0=np.ones(K)+1j*np.ones(K)
    theta0N0=np.zeros(K)
    phiX0N0=np.ones(K)+1j*np.ones(K)
    
    cp0N1=np.ones(K)+1j*np.ones(K)
    theta0N1=np.zeros(K)
    phiX0N1=np.ones(K)+1j*np.ones(K)
    
    cp0N2=np.ones(K)+1j*np.ones(K)
    theta0N2=np.zeros(K)
    phiX0N2=np.ones(K)+1j*np.ones(K)
    
    cp0N3=np.ones(K)+1j*np.ones(K)
    theta0N3=np.zeros(K)
    phiX0N3=np.ones(K)+1j*np.ones(K)
    
    cp0N4=np.ones(K)+1j*np.ones(K)
    theta0N4=np.zeros(K)
    phiX0N4=np.ones(K)+1j*np.ones(K)
    
    cp0N5=np.ones(K)+1j*np.ones(K)
    theta0N5=np.zeros(K)
    phiX0N5=np.ones(K)+1j*np.ones(K)
    
    cp1N0=np.ones(K)+1j*np.ones(K)
    theta1N0=np.zeros(K)
    phiX1N0=np.ones(K)+1j*np.ones(K)
    
    cp1N1=np.ones(K)+1j*np.ones(K)
    theta1N1=np.zeros(K)
    phiX1N1=np.ones(K)+1j*np.ones(K)
    
    #test x value to get correct phases for bloch waves
    xteN0=0
    xteN1=1+xteN0
    xteN2=2+xteN0
    xteN3=3+xteN0
    xteN4=4+xteN0
    xteN5=5+xteN0
    
    
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
    
    
    
        #get correct phases for ground state, depends on x=0 value
        phiX0N0[qi], theta0N0[qi], cp0N0[qi] = cPhaseEven(xteN0,L,qs[qi],BF0[:,qi])
        
        #get correct phases for ground state, depends on x=1 value        
        phiX0N1[qi], theta0N1[qi], cp0N1[qi] = cPhaseEven(xteN1,L,qs[qi],BF0[:,qi])
        
        #get correct phases for ground state, depends on x=2 value        
        phiX0N2[qi], theta0N2[qi], cp0N2[qi] = cPhaseEven(xteN2,L,qs[qi],BF0[:,qi])
        
        phiX0N3[qi], theta0N3[qi], cp0N3[qi] = cPhaseEven(xteN3,L,qs[qi],BF0[:,qi])
        
        phiX0N4[qi], theta0N4[qi], cp0N4[qi] = cPhaseEven(xteN4,L,qs[qi],BF0[:,qi])
        
        phiX0N5[qi], theta0N5[qi], cp0N5[qi] = cPhaseEven(xteN5,L,qs[qi],BF0[:,qi])
        
        
        #get correct phases for 1st excited state, depends on d/dx @ x=0,1
        phiX1N0[qi], theta1N0[qi], cp1N0[qi] = cPhaseOdd(xteN0,L,qs[qi],BF1[:,qi])
        
        phiX1N1[qi], theta1N1[qi], cp1N1[qi] = cPhaseOdd(xteN1,L,qs[qi],BF1[:,qi])
        
        #get correct phases for 1st excited state, depends on d/dx @ x=1
#        phiX1N1[qi]= (BFM(xtEN1,L,qs[qi],BF1[:,qi]) - \
#                    BFM(-xtEN1,L,qs[qi],BF1[:,qi]))/(2*xtEN1)
#        theta1N1[qi]=np.angle(phiX1N1[qi])
#        cp1N1[qi]=np.exp(-1j*theta1N1[qi])
        
        #get correct phases for ground state, depends on x=1 value
#        phiX0N1[qi]=BFM(xtEN1,L,qs[qi],BF0[:,qi])
#        theta0N1[qi]=np.angle(phiX0N1[qi])
#        cp0N1[qi]=np.exp(-1j*theta0N1[qi])
#        

    
    #initialize section
    xx=np.linspace(-NS,NS,NX)
    dx=xx[1]-xx[0]
    
    cntr=0
    wn0n0=np.ones(NX)+1j*np.ones(NX)
    wn0n1=np.ones(NX)+1j*np.ones(NX)
    wn0n2=np.ones(NX)+1j*np.ones(NX)
    wn0n3=np.ones(NX)+1j*np.ones(NX)
    wn0n4=np.ones(NX)+1j*np.ones(NX)
    wn0n5=np.ones(NX)+1j*np.ones(NX)
    
    wn1n0=np.ones(NX)+1j*np.ones(NX)
    wn1n1=np.ones(NX)+1J*np.ones(NX)
    
    #find x values of Wannier function with correct phases
    for x in xx:
        wn0n0[cntr]=w0(BF0,cp0N0,x,K,L)
        wn0n1[cntr]=w0(BF0,cp0N1,x,K,L)
        wn0n2[cntr]=w0(BF0,cp0N2,x,K,L)
        wn0n3[cntr]=w0(BF0,cp0N3,x,K,L)
        wn0n4[cntr]=w0(BF0,cp0N4,x,K,L)
        wn0n5[cntr]=w0(BF0,cp0N5,x,K,L)
        
        wn1n0[cntr]=w0(BF1,cp1N0,x,K,L)
        wn1n1[cntr]=w0(BF1,cp1N1,x,K,L)
        cntr+=1
    
    #normalize wannier Functions, I'm missing some factor that depends on L or K above
    wn0n0 = normMe(wn0n0)/np.sqrt(dx)
    wn0n1 = normMe(wn0n1)/np.sqrt(dx)
    wn0n2 = normMe(wn0n2)/np.sqrt(dx)
    wn0n3 = normMe(wn0n3)/np.sqrt(dx)
    wn0n4 = normMe(wn0n4)/np.sqrt(dx)
    wn0n5 = normMe(wn0n5)/np.sqrt(dx)
    
    wn1n0 = normMe(wn1n0)/np.sqrt(dx)
    wn1n1 = normMe(wn1n1)/np.sqrt(dx)
    
    return wn0n0, wn0n1, wn0n2, wn0n3, wn0n4, wn0n5, wn1n0, wn1n1 #, BF0

K=35
L=K-1
Er=1.24
NS=L/2
NX=3000

xx=np.linspace(-NS,NS,NX)
dx=xx[1]-xx[0]
xinds=[xi for xi in range(len(xx)) if (xx[xi]>=-0.5 and xx[xi]<0.5)]
xinds=[xi for xi in range(len(xx)) if (xx[xi]>=-100 and xx[xi]<100)]

VoN=np.linspace(2,8,3)
VoN=[2,4]


#
wn0n0Sv=np.ones((len(VoN),NX))*1j
wn0n1Sv=np.ones((len(VoN),NX))*1j
wn0n2Sv=np.ones((len(VoN),NX))*1j
wn0n3Sv=np.ones((len(VoN),NX))*1j
wn0n4Sv=np.ones((len(VoN),NX))*1j
wn0n5Sv=np.ones((len(VoN),NX))*1j

wn1n0Sv=np.ones((len(VoN),NX))*1j
wn1n1Sv=np.ones((len(VoN),NX))*1j

cntr=0

J1_n00=np.ones(len(VoN))*1j
J2_n00=np.ones(len(VoN))*1j
J3_n00=np.ones(len(VoN))*1j
J4_n00=np.ones(len(VoN))*1j
J5_n00=np.ones(len(VoN))*1j

J1_n11=np.ones(len(VoN))*1j

U_n00=np.ones(len(VoN))*1j
U_n11=np.ones(len(VoN))*1j

E_n00=np.ones(len(VoN))*1j
E_n11=np.ones(len(VoN))*1j


for Vo in VoN:
    print(Vo)
    
    wn0n0, wn0n1, wn0n2, wn0n3, wn0n4, wn0n5, wn1n0, wn1n1= calWannier(L, K, Vo, 1.24, 0, NS, NX)
    
    wn0n0Sv[cntr,::]=wn0n0
    wn0n1Sv[cntr,::]=wn0n1
    wn0n2Sv[cntr,::]=wn0n2
    wn0n3Sv[cntr,::]=wn0n3
    wn0n4Sv[cntr,::]=wn0n4
    wn0n5Sv[cntr,::]=wn0n5
    
    wn1n0Sv[cntr,::]=wn0n0
    wn1n1Sv[cntr,::]=wn0n1
    
    #d_wn0n1=wn0n1[1:NX]-wn0n1[0:(NX-1)]
    #d2_wn0n1=d_wn0n1[1:(NX-1)]-d_wn0n1[0:(NX-2)]
    
    #d2_wn0n1 = D2(wn0n1,xx)
    
    Vp=-Vo*np.cos(xx*2*np.pi)/2
    
    J1_n00[cntr]=(calBHJ(wn0n1,wn0n0,Vp,xx,xinds))
    print(J1_n00[cntr]*Er)
    J2_n00[cntr]=(calBHJ(wn0n0,wn0n2,Vp,xx,xinds))
    J3_n00[cntr]=np.real(calBHJ(wn0n0,wn0n3,Vp,xx,xinds))
    J4_n00[cntr]=np.real(calBHJ(wn0n0,wn0n4,Vp,xx,xinds))
    J5_n00[cntr]=np.real(calBHJ(wn0n0,wn0n5,Vp,xx,xinds))
    
    J1_n11[cntr]=np.real(calBHJ(wn1n0,wn1n1,Vp,xx,xinds))
    
    U_n00[cntr]=np.real(calBHU(wn0n0,xx,xinds))
    U_n11[cntr]=np.real(calBHU(wn1n0,xx,xinds))
    
    E_n00[cntr]=np.real(calBHE(wn0n0,Vp,xx,xinds))
    E_n11[cntr]=np.real(calBHE(wn1n0,Vp,xx,xinds))

    cntr+=1
 
    
plt.plot(xx,Vp)
plt.plot(xx,wn0n0Sv[1,::])
plt.plot(xx,wn0n1Sv[1,::])
plt.plot(xx,xx*0)
plt.xlim([-3.5,3.5])
plt.show()


plt.plot(VoN,np.real(J1_n00)*Er)
plt.plot(VoN,np.real(J2_n00)*Er)
plt.show()
#    
#    J1_n00[cntr]=np.real(-np.dot(np.conj(wn0n0),D2(wn0n1,xx))+np.dot(np.conj(wn0n0),Vp*wn0n1))
#    J1_n00[cntr]=np.real(calBHI(wn0n0,wn0n1,Vp,xx))
#    J1_n00[cntr]=
#    
#    
#    U_n00=np.dot(np.conj(wn0n0)*np.conj(wn0n0),wn0n0*wn0n0)



#print(np.real(J2_n00)/np.real(J1_n00))
#print(np.real(J3_n00)/np.real(J1_n00))

#etilt=np.linspace(-1,1,100)
#J1_n01=[]
#for et in etilt:
#    Ve=Vo*(xx)*et
#    J1_n01.append(-np.dot(np.conj(wn0n0),D2(wn1n1,xx))+np.dot(np.conj(wn0n0),(Vp+Ve)*wn1n1))


#xxout=xx.reshape((len(xx),1))
#
#plt.plot(xx,5*(1-np.cos(xx*np.pi)**2)-0.5,'k-')
#
##Vout=(1-np.cos(xx*np.pi)**2).reshape((len(xx),1))
##W0out=wn0.reshape((len(xx),1))
##W1out=wn1.reshape((len(xx),1))
##W2out=wn2.reshape((len(xx),1))
##W3out=wn3.reshape((len(xx),1))
#
#plt.plot(xx,3*wn0n0+0)
#plt.plot(xx,3*wn0n1+1)
#plt.plot(xx,3*wn0n2+2)
#plt.plot(xx,3*wn1n0+3)
#plt.plot(xx,3*wn1n1+4)
#plt.xlim([-2.5,2.5])
#plt.show()
#
#plt.plot(etilt,np.array(J1_n01)*1.24)
#np.savetxt("LattOut.csv",np.hstack((xxout,Vout)),delimiter=",")
#np.savetxt("W0_1Er.csv",np.hstack((xxout,W0out)),delimiter=",")
#np.savetxt("W1_1Er.csv",np.hstack((xxout,W1out)),delimiter=",")
#np.savetxt("W2_1Er.csv",np.hstack((xxout,W2out)),delimiter=",")
#np.savetxt("W3_1Er.csv",np.hstack((xxout,W3out)),delimiter=",")
