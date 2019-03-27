# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 17:31:04 2019

@author: Matthew
"""

import numpy as np
#from scipy.linalg import expm
import matplotlib.pyplot as plt

# set some colors for plotting

maticaA = (0.65,0.00,0.00)
maticaB = (0.05,0.53,0.63)
maticaC = (0.44,0.26,0.71)
maticaD = (0.75,0.36,0.13)
maticaE = (0.46,0.56,0.01)
maticaF = (0.66,0.21,0.30)
maticaG = (0.21,0.39,0.80)
maticaH = (0.78,0.52,0.04)
maticaI = (0.52,0.22,0.53)

clrWhite = (1,1,1)

def blendClr(color1,color2,alphaVal):
    return np.array(color1)*alphaVal + np.array(color2)*(1-alphaVal)



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
    cntr=0
    bf0=np.ones((NX,K))+1j*np.ones((NX,K))
    bf1=np.ones((NX,K))+1j*np.ones((NX,K))
    bf2=np.ones((NX,K))+1j*np.ones((NX,K))
    bf3=np.ones((NX,K))+1J*np.ones((NX,K))
    
    wn0=np.ones(NX)+1J*np.ones(NX)
    wn1=np.ones(NX)+1J*np.ones(NX)
    wn2=np.ones(NX)+1J*np.ones(NX)
    wn3=np.ones(NX)+1J*np.ones(NX)
    
    #find x values of bloch function with correct phases
    for x in xx:
        bf0[cntr,::]=bf(BF0,cp0,x*2*np.pi,K,L)
        bf1[cntr,::]=bf(BF1,cp1,x*2*np.pi,K,L)
        bf2[cntr,::]=bf(BF2,cp2,x*2*np.pi,K,L)
        bf3[cntr,::]=bf(BF3,cp3,x*2*np.pi,K,L)
        
        wn0[cntr]=w0(BF0,cp0,x*2*np.pi,K,L)
        wn1[cntr]=w0(BF1,cp1,x*2*np.pi,K,L)
        wn2[cntr]=w0(BF2,cp2,x*2*np.pi,K,L)
        wn3[cntr]=w0(BF3,cp3,x*2*np.pi,K,L)
        cntr+=1
    
    #normalize height of bloch functions to 1
    for qi in range(K):
        bf0[::,qi] = normMeBF(bf0[::,qi])
        bf1[::,qi] = normMeBF(bf1[::,qi])
        bf2[::,qi] = normMeBF(bf2[::,qi])
        bf3[::,qi] = normMeBF(bf3[::,qi])
    
    return BS0, BS1, BS2, BS3, bf0, bf1, bf2, bf3, wn0, wn1, wn2, wn3

NS=4.5
NX=1000
L=35
K=32*4 # pick a multiple of 2^n if you want to show periodicity of n*lattice sites nicely
Er=1.24
Vo=4

qs=np.linspace(-1.0/2,1.0/2-1.0/K,K)

BS0, BS1, BS2, BS3, bf0, bf1, bf2, bf3, wn0, wn1, wn2, wn3 = calBlochFx(L, K, Vo, Er, 0, NS, NX)

xx=np.linspace(-NS,NS,NX)
xxout=xx.reshape((len(xx),1))

fig = plt.figure(figsize=(1.618*6*1.2,1*6))

#subplot of excited band
        
# bloch functions of ground band
ax221 = fig.add_subplot(221, aspect='1')

ax221.plot(xx,1*(1-np.cos(xx*np.pi)**2)-0.5,'-',color = 'black', alpha=0.3, zorder = -10)
ax221.plot(xx,1*(1-np.cos(xx*np.pi)**2)+0.5,'-',color = 'black', alpha=0.3, zorder = -10)

Vout=(1-np.cos(xx*np.pi)**2).reshape((len(xx),1))
ax221.plot(xx,np.real(bf1[::,K/2])/3,'-', color = maticaB)
ax221.plot(xx,np.imag(bf1[::,K/2])/3,'--', color = maticaB)

ax221.plot(xx,np.real(bf1[::,K/2+K/8+1])/3+1,'-', color = maticaC)
ax221.plot(xx,np.imag(bf1[::,K/2+K/8+1])/3+1,'--', color = maticaC)

ax221.plot(xx,np.real(bf1[::,K/2+K/4+1])/3+2,'-', color = maticaA)
ax221.plot(xx,np.imag(bf1[::,K/2+K/4+1])/3+2,'--', color = maticaA)

ax221.plot(xx,np.real(bf1[::,0])/3+3,'-', color = maticaD)
ax221.plot(xx,np.imag(bf1[::,0])/3+3,'--', color = maticaD)
ax221.set_title('Bloch Wavefunctions: n=1')
ax221.set_yticks([0,1,2,3])
ax221.set_yticklabels(['q=0','q=1/8','q=1/4','q=1/2'])
   

# bloch functions of ground band
ax223 = fig.add_subplot(223, aspect='1')

ax223.plot(xx,4*(1-np.cos(xx*np.pi)**2)-0.5,'-',color = 'black', alpha=0.3, zorder = -10)

Vout=(1-np.cos(xx*np.pi)**2).reshape((len(xx),1))

bfplts =[]

bfplts+=ax223.plot(xx,np.real(bf0[::,K/2])/3,'-', color = maticaB, label = 'Re[f,q]')
bfplts+=ax223.plot(xx,np.imag(bf0[::,K/2])/3,'--', color = maticaB, label = 'Im[f,q]')

bfplts+=ax223.plot(xx,np.real(bf0[::,K/2+K/8])/3+1,'-', color = maticaC)
bfplts+=ax223.plot(xx,np.imag(bf0[::,K/2+K/8])/3+1,'--', color = maticaC)

bfplts+=ax223.plot(xx,np.real(bf0[::,K/2+K/4])/3+2,'-', color = maticaA)
bfplts+=ax223.plot(xx,np.imag(bf0[::,K/2+K/4])/3+2,'--', color = maticaA)

bfplts+=ax223.plot(xx,np.real(bf0[::,0])/3+3,'-', color = maticaD)
bfplts+=ax223.plot(xx,np.imag(bf0[::,0])/3+3,'--', color = maticaD)

ax223.set_xlabel('Lattice Sites')
ax223.set_title('Bloch Wavefunctions: n=0')
ax223.set_yticks([0,1,2,3])
ax223.set_yticklabels(['q=0','q=1/8','q=1/4','q=1/2'])
ax223.legend(frameon=False, loc=(0,-0.5))

# subplot of bands
ax143 = fig.add_subplot(122, aspect='auto')

ax143.plot([qs[K/2],qs[K/2]],[-12,20], '--', \
           color = maticaB, alpha=0.5, linewidth=2, zorder = -10)
ax143.plot([qs[K/2+K/8+1],qs[K/2+K/8+1]],[-12,20], '-', \
           color = maticaC, alpha=0.5, linewidth=2, zorder = -10)
ax143.plot([qs[K/2+K/4+1],qs[K/2+K/4+1]],[-12,20], '-.', \
           color = maticaA, alpha=0.5, linewidth=2, zorder = -10)
ax143.plot([0.5-0.005,0.5-0.005],[-12,20], ':', \
           color = maticaD, alpha=0.5, linewidth=2, zorder = -10)

ax143.plot(qs,BS3, color = matica97D, linewidth=3, label = 'n=3')
ax143.plot(qs,BS2, color = matica97C, linewidth=3, label = 'n=2')
ax143.plot(qs,BS1, color = matica97B, linewidth=3, label = 'n=1')
ax143.plot(qs,BS0, color = matica97A, linewidth=3, label = 'n=0')

ax143.legend(frameon=False,loc=(1.1,0.4))
ax143.set_xlim((-0.5,0.5))
ax143.set_xticks(np.linspace(-0.5,0.5,5))
ax143.set_xticklabels(np.linspace(-0.5,0.5,5))

ax143.set_ylim([-12,20])
ax143.set_yticks(np.linspace(-12,20,9))
ax143.set_yticklabels(np.linspace(-12,20,9))
ax143.set_ylabel('E_q^n (Er)')
ax143.set_xlabel('q')
ax143.set_title('Lattice Depth: '+str(int(Vo))+' Er')

#plt.savefig('BlochFunctionsBS.pdf')

#np.savetxt("BSLattDepth.csv",VoN,delimiter=",")
np.savetxt("BS0_"+str(int(Vo))+"_Er.csv",BS0,delimiter=",")
np.savetxt("BS1_"+str(int(Vo))+"_Er.csv",BS1,delimiter=",")
np.savetxt("BS2_"+str(int(Vo))+"_Er.csv",BS2,delimiter=",")
np.savetxt("BS3_"+str(int(Vo))+"_Er.csv",BS3,delimiter=",")

np.savetxt("BF0_q0_"+str(int(Vo))+"_Er.csv",bf0[::,K/2],delimiter=",")
np.savetxt("BF0_q0p125_"+str(int(Vo))+"_Er.csv",bf0[::,K/2+K/8],delimiter=",")
np.savetxt("BF0_0p25_"+str(int(Vo))+"_Er.csv",bf0[::,K/2+K/4],delimiter=",")
np.savetxt("BF0_0p5_"+str(int(Vo))+"_Er.csv",bf0[::,0],delimiter=",")

np.savetxt("BF1_q0_"+str(int(Vo))+"_Er.csv",bf1[::,K/2],delimiter=",")
np.savetxt("BF1_q0p125_"+str(int(Vo))+"_Er.csv",bf1[::,K/2+K/8],delimiter=",")
np.savetxt("BF1_0p25_"+str(int(Vo))+"_Er.csv",bf1[::,K/2+K/4],delimiter=",")
np.savetxt("BF1_0p5_"+str(int(Vo))+"_Er.csv",bf1[::,0],delimiter=",")

np.savetxt("BSQs.csv",qs,delimiter=",")


# subplot of ground band
#ax144 = fig.add_subplot(133, aspect='auto')
#
#ax144.plot(xx,4*(1-np.cos(xx*np.pi)**2)-0.5,'-',color = 'black', alpha=0.3, zorder = -10)
#ax144.plot(xx,wn0/2+0)
#ax144.plot(xx,wn1/2+1)
#ax144.plot(xx,wn2/2+2)
#ax144.plot(xx,wn3/2+3)


#
#ax231 = fig.add_subplot(4,3,4, aspect='auto')
#
#ax231.plot(xx,4*(1-np.cos(xx*np.pi)**2)-0.5,'-',color = 'black', alpha=0.3, zorder = -10)
#
#Vout=(1-np.cos(xx*np.pi)**2).reshape((len(xx),1))
#ax231.plot(xx,np.real(bf2[::,K/2])/3,'-', color = maticaB)
#ax231.plot(xx,np.imag(bf2[::,K/2])/3,'--', color = maticaB)
#
#ax231.plot(xx,np.real(bf2[::,K/2+K/8+1])/3+1,'-', color = maticaC)
#ax221.plot(xx,np.imag(bf2[::,K/2+K/8+1])/3+1,'--', color = maticaC)
#
#ax231.plot(xx,np.real(bf2[::,K/2+K/4+1])/3+2,'-', color = maticaA)
#ax231.plot(xx,np.imag(bf2[::,K/2+K/4+1])/3+2,'--', color = maticaA)
#
#ax231.plot(xx,np.real(bf2[::,0])/3+3,'-', color = maticaD)
#ax231.plot(xx,np.imag(bf2[::,0])/3+3,'--', color = maticaD)
#
#ax231 = fig.add_subplot(4,3,7, aspect='auto')
#
#ax231.plot(xx,4*(1-np.cos(xx*np.pi)**2)-0.5,'-',color = 'black', alpha=0.3, zorder = -10)
#
#Vout=(1-np.cos(xx*np.pi)**2).reshape((len(xx),1))
#ax231.plot(xx,np.real(bf1[::,K/2])/3,'-', color = maticaB)
#ax231.plot(xx,np.imag(bf1[::,K/2])/3,'--', color = maticaB)
#
#ax231.plot(xx,np.real(bf1[::,K/2+K/8+1])/3+1,'-', color = maticaC)
#ax221.plot(xx,np.imag(bf1[::,K/2+K/8+1])/3+1,'--', color = maticaC)
#
#ax231.plot(xx,np.real(bf1[::,K/2+K/4+1])/3+2,'-', color = maticaA)
#ax231.plot(xx,np.imag(bf1[::,K/2+K/4+1])/3+2,'--', color = maticaA)
#
#ax231.plot(xx,np.real(bf1[::,0])/3+3,'-', color = maticaD)
#ax231.plot(xx,np.imag(bf1[::,0])/3+3,'--', color = maticaD)


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
#np.savetxt("LattOut.csv",np.hstack((xxout,Vout)),delimiter=",")
#np.savetxt("W0_1Er.csv",np.hstack((xxout,W0out)),delimiter=",")
#np.savetxt("W1_1Er.csv",np.hstack((xxout,W1out)),delimiter=",")
#np.savetxt("W2_1Er.csv",np.hstack((xxout,W2out)),delimiter=",")
#np.savetxt("W3_1Er.csv",np.hstack((xxout,W3out)),delimiter=",")
