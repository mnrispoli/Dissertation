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

#experimental parameters (all in SI units)
gamma=2*np.pi*6.065*(1E6)
lambda_laser = 758*(1E-9)
lambda_rb_d2 = 780.241209686*(1E-9)
c = 299792458
flaser = c/lambda_laser
frbd2 = c/lambda_rb_d2
delta_rb_d2 = 2*np.pi*(flaser - frbd2)

def blueDetHeating(vn,wn,xx,Er,gamma,delta):
    dx=xx[2]-xx[1]
    return np.dot(vn,abs(wn)**2)*gamma*Er*dx/delta

def PsiHO(n,x,a):
#    x=x
    return ((np.sqrt(n)*4.0*np.pi/(a**2))**(1.0/4.0))* \
            np.exp(-((2.0*np.pi*x)**2.0)*np.sqrt(n)/2.0/(a**2)) #/ \
#            np.sqrt(2*np.pi)
            
def daleyExp(n,Er,gamma,delta):
    return gamma*Er*n/delta/np.sqrt(4.0*n)



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
    return np.exp(1j*x*np.linspace(-lh,lh,L+1))/np.sqrt(L+1)

def BFM(x,L,q,BF):
    #add the bragg scattered  waves + phase shift of iqx plus contributions from
    #hybridized bloch waves found from diagonalized hamiltonian
    return np.dot(BF[::],FM(x,L)*np.exp(1j*x*q))/np.sqrt(2*np.pi*(L+1))

def wFx(BF,cp,x,qs,K,L):
    #create wannier function from all bloch waves that depend on q
    
    #initialize output
    wout=0+0*1j
    
    #q values
#    qs=np.linspace(-1.0/2,1.0/2-1.0/K,K)
    
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
    
    d1=(xvec[1:NX]-xvec[0:(NX-1)])/dx
    
    d2=(d1[1:(NX-1)]-d1[0:(NX-2)])/dx
    

    
    return np.append( \
                     np.append(np.array([d2[0]]), \
                     np.array([d2])), \
                     np.array([d2[-1]]))

def calBHJ(psi1,psi2,Vext,xx):
    dx=xx[1]-xx[0]  

    xx=np.reshape(xx,psi1.shape)
    Vext=np.reshape(Vext,psi1.shape)
    
#    psi1=psi1[xinds]
#    psi2=psi2[xinds]
#    Vext=Vext[xinds]
#    xx=xx[xinds]
    
    J_1=-4*np.dot(np.conj(psi1),D2(psi2,xx))*dx
    J_2=np.dot(np.conj(psi1),Vext*psi2)*dx
    
    return -(J_1+J_2)

def calBHU(psi,xx):
    dx=xx[1]-xx[0]
    xx=np.reshape(xx,psi.shape)
    

#    psi=psi[xinds]
#    xx=xx[xinds]
    return (np.dot(np.conj(psi)*np.conj(psi),psi*psi))*dx

def calBHE(psi,Vext,xx):
    dx=xx[1]-xx[0]
    xx=np.reshape(xx,psi.shape)
    Vext=np.reshape(Vext,psi.shape)
    

#    psi=psi[xinds]
#    Vext=Vext[xinds]
#    xx=xx[xinds]
    E_1=-4*np.dot(np.conj(psi),D2(psi,xx))*dx
    E_2=np.dot(np.conj(psi),Vext*psi)*dx

    return (E_1+E_2)

def calWannier(L, K, VoJ, Er, eps, NS, NX, nn):

    #all quasi-momentums states that will be considered
    qs=np.linspace(-1.0/2,1.0/2-1.0/K,K)
    
    
    # initialize a ton of stuff for first 4 bands in lattice   
    BS0=np.zeros(K)
    BS1=np.zeros(K)
    BS2=np.zeros(K)
    BS3=np.zeros(K)
    
    BF0=np.zeros([L+1,K])
    BF1=np.zeros([L+1,K])
    BF2=np.zeros([L+1,K])
    BF3=np.zeros([L+1,K])
    
    #things needed for making wannier fx
    #someone should just fix this to be a matrix
    
    #nearest neighbors to calculate for first 4 bands
    NN0=nn+1
    cp0=np.ones((NN0,K))+1j*np.ones((NN0,K))
    theta0=np.zeros((NN0,K))
    phiX0=np.ones((NN0,K))+1j*np.ones((NN0,K))
    
#    NN1=6+1
    cp1=np.ones((NN0,K))+1j*np.ones((NN0,K))
    theta1=np.zeros((NN0,K))
    phiX1=np.ones((NN0,K))+1j*np.ones((NN0,K))
    
#    NN2=6+1
    cp2=np.ones((NN0,K))+1j*np.ones((NN0,K))
    theta2=np.zeros((NN0,K))
    phiX2=np.ones((NN0,K))+1j*np.ones((NN0,K))
    
#    NN3=6+1
    cp3=np.ones((NN0,K))+1j*np.ones((NN0,K))
    theta3=np.zeros((NN0,K))
    phiX3=np.ones((NN0,K))+1j*np.ones((NN0,K))
    
  
    #test x value to get correct phases for bloch waves
    xteN=[0]
    
#    nh=(NN0-1)/2
#    xteN=np.linspace(-nh,nh,NN0)
    xteN=np.linspace(0,NN0-1,NN0)
    
    for qi in range(K):
        
        q=qs[qi]
        
        #create single-particle hamiltonian in lattice
        H=HamJQ(q,VoJ,L+1)       
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
    
    
    
        #get correct phases for ground state, depends on x=xi value
        cntr=0
        for xi in xteN:
            phiX0[cntr,qi], theta0[cntr,qi], cp0[cntr,qi] = cPhaseEven(xi*2.0*np.pi,L,qs[qi],BF0[:,qi])
#            phiX1[cntr,qi], theta1[cntr,qi], cp1[cntr,qi] = cPhaseOdd(xi*2.0*np.pi,L,qs[qi],BF1[:,qi])
#            phiX2[cntr,qi], theta2[cntr,qi], cp2[cntr,qi] = cPhaseEven(xi*2.0*np.pi,L,qs[qi],BF2[:,qi])
#            phiX3[cntr,qi], theta3[cntr,qi], cp3[cntr,qi] = cPhaseOdd(xi*2.0*np.pi,L,qs[qi],BF3[:,qi])
            cntr+=1
        
    #initialize section
    xx=np.linspace(-NS,NS,NX)
    
    cntr=0
    wn0=np.ones((NN0,NX))+1j*np.ones((NN0,NX))
#    wn1=np.ones((NN0,NX))+1j*np.ones((NN0,NX))
#    wn2=np.ones((NN0,NX))+1j*np.ones((NN0,NX))
#    wn3=np.ones((NN0,NX))+1j*np.ones((NN0,NX))

    #find x values of Wannier function with correct phases
    for x in xx:
        for nn in range(NN0):
            wn0[nn,cntr]=wFx(BF0,cp0[nn,::],2*np.pi*x,qs,K,L)
#            wn1[nn,cntr]=wFx(BF1,cp1[nn,::],2*np.pi*x,qs,K,L)
#            wn2[nn,cntr]=wFx(BF2,cp2[nn,::],2*np.pi*x,qs,K,L)
#            wn3[nn,cntr]=wFx(BF3,cp3[nn,::],2*np.pi*x,qs,K,L)

        cntr+=1
    
    return wn0 #, wn1, wn2, wn3 



### main code section
    # set parametersfor evaluation of Bose-Hubbard parameters
K=35 #35 # number of quasimomentum states evaluated for
L=K-1 # "length" of lattice, really accounts for the number of considered free 
        #particle dispersion curves that couple in the prerioid brioullin zone scheme
        # linked to number of q-states you want to evaluate -> 1/L = dq
        
Er=1.24*(1E3) # lattice recoil depth in kHz for Rb lab
Er_ax = 255
Er_big = 7
NS=L/2 # number of sites to evaluate 
NX=5001 # number of discrete points to evaluate 


#initialize real-space points
xx=np.linspace(-NS,NS,NX)*2*np.pi
xx=np.reshape(xx,(NX,1))
dx=xx[1]-xx[0]

#counter
cntr=0
VoN=np.linspace(0.1,100,5)  #lattice depth to evaluate for
VoN=np.logspace(-1,5,151)
#initialize save files for wannier fx
wn0Sv=xx
wn1Sv=xx
wn2Sv=xx
wn3Sv=xx

nn=0 #number of neighbors to find for tunneling matrix elements

bhJ0s=np.zeros((nn,len(VoN)))
bhU0s=np.zeros((len(VoN),1))
bhE0s=np.zeros((len(VoN),1))

bhJ1s=np.zeros((nn,len(VoN)))
bhU1s=np.zeros((len(VoN),1))
bhE1s=np.zeros((len(VoN),1))

bhJ2s=np.zeros((nn,len(VoN)))
bhU2s=np.zeros((len(VoN),1))
bhE2s=np.zeros((len(VoN),1))

bhJ3s=np.zeros((nn,len(VoN)))
bhU3s=np.zeros((len(VoN),1))
bhE3s=np.zeros((len(VoN),1))

Gsc = np.zeros(len(VoN))
GscHO = np.zeros(len(VoN))

Gsc_ax = np.zeros(len(VoN))
GscHO_ax = np.zeros(len(VoN))

Gsc_big = np.zeros(len(VoN))
GscHO_big = np.zeros(len(VoN))

axialOff = 0.05644
bigOff = 0.00165

for Vo in VoN:

    print(Vo)
    wn0 = calWannier(L, K, Vo*1.0, 1.24, 0, NS, NX, nn)
    #, wn1, wn2, wn3
    
    Vp=-Vo*np.cos(xx)/2
    Vext=Vp
    wn0=np.transpose(wn0)
    psiho0=PsiHO(Vo*1.0,xx,2*np.pi*2)
    
    GscHO[cntr]=blueDetHeating(np.transpose(Vext+Vo/2),psiho0,xx,Er,gamma,delta_rb_d2)
    Gsc[cntr]=blueDetHeating(np.transpose(Vext+Vo/2),wn0,xx,Er,gamma,delta_rb_d2)
    
    GscHO_ax[cntr]=blueDetHeating(np.transpose(Vext+Vo/2+Vo*axialOff),psiho0,xx,Er_ax,gamma,delta_rb_d2)
    Gsc_ax[cntr]=blueDetHeating(np.transpose(Vext+Vo/2+Vo*axialOff),wn0,xx,Er_ax,gamma,delta_rb_d2)
    
    GscHO_big[cntr]=blueDetHeating(np.transpose(Vext+Vo/2+Vo*bigOff),psiho0,xx,Er_big,gamma,delta_rb_d2)
    Gsc_big[cntr]=blueDetHeating(np.transpose(Vext+Vo/2+Vo*bigOff),wn0,xx,Er_big,gamma,delta_rb_d2)
    cntr+=1        
      
            
FX=6.5*2
FY=FX/3#/1.618

[XL,XH,YL,YH] = [0,0.25,0,4.0]
[XW,YT]=[XH-XL,YH-YL]

W=XW/2
H=FY
AR=XW/YT

fig = plt.figure(figsize=(FX, FY))

ax1 = fig.add_subplot(121, aspect='auto')
ax1.plot(xx/2/np.pi,abs(wn0))   

ax1.plot(xx/2/np.pi,psiho0)   
ax1.set_xlim([-1,1])

ax2 = fig.add_subplot(122, aspect='auto')
ax2.loglog(VoN,abs(Gsc)*2.0*np.pi,'-')
ax2.loglog(VoN,abs(GscHO)*2.0*np.pi,'-')
ax2.loglog(VoN,abs(daleyExp(VoN,Er,gamma,delta_rb_d2))*2*np.pi,'-')

ax2.loglog(VoN,abs(GscHO_ax)*2.0*np.pi,'--')
ax2.loglog(VoN,abs(Gsc_ax)*2.0*np.pi,'--')
ax2.loglog(VoN,abs(daleyExp(VoN,Er_ax,gamma,delta_rb_d2))*2*np.pi,'--')

ax2.loglog(VoN,abs(GscHO_big)*2.0*np.pi,'.-')
ax2.loglog(VoN,abs(Gsc_big)*2.0*np.pi,'.-')
ax2.loglog(VoN,abs(daleyExp(VoN,Er_big,gamma,delta_rb_d2))*2*np.pi,'.-')
ax2.grid()



np.savetxt("2DScatterWn.csv",\
           np.vstack( \
                     (np.reshape(VoN,(1,len(VoN))),\
                      np.reshape(Gsc*2*np.pi,(1,len(Gsc))))\
                      ) \
           ,delimiter=",")
                     
np.savetxt("2DScatterHO.csv",\
           np.vstack( \
                     (np.reshape(VoN,(1,len(VoN))),\
                      np.reshape(GscHO*2*np.pi,(1,len(GscHO))))\
                      ) \
           ,delimiter=",")
GscDE = abs(daleyExp(VoN,Er,gamma,delta_rb_d2))                   
np.savetxt("2DScatterDaleyEst.csv",\
           np.vstack( \
                     (np.reshape(VoN,(1,len(VoN))),\
                      np.reshape(GscDE*2*np.pi,(1,len(GscDE))))\
                      ) \
           ,delimiter=",")
                     
                     
                     
                     
np.savetxt("AxialScatterWn.csv",\
           np.vstack( \
                     (np.reshape(VoN,(1,len(VoN))),\
                      np.reshape(Gsc_ax*2*np.pi,(1,len(Gsc_ax))))\
                      ) \
           ,delimiter=",")
                     
np.savetxt("AxialScatterHO.csv",\
           np.vstack( \
                     (np.reshape(VoN,(1,len(VoN))),\
                      np.reshape(GscHO_ax*2*np.pi,(1,len(GscHO_ax))))\
                      ) \
           ,delimiter=",")
                     
GscDE_ax= abs(daleyExp(VoN,Er_ax,gamma,delta_rb_d2))                    
np.savetxt("AxialScatterDaleyEst.csv",\
           np.vstack( \
                     (np.reshape(VoN,(1,len(VoN))),\
                      np.reshape(GscDE_ax*2*np.pi,(1,len(GscDE_ax))))\
                      ) \
           ,delimiter=",")
                     
                     
                     
                     
                     
np.savetxt("BigScatterWn.csv",\
           np.vstack( \
                     (np.reshape(VoN,(1,len(VoN))),\
                      np.reshape(Gsc_big*2*np.pi,(1,len(Gsc_big))))\
                      ) \
           ,delimiter=",")
                     
np.savetxt("BigScatterHO.csv",\
           np.vstack( \
                     (np.reshape(VoN,(1,len(VoN))),\
                      np.reshape(GscHO_big*2*np.pi,(1,len(GscHO_big))))\
                      ) \

           ,delimiter=",")
GscDE_big= abs(daleyExp(VoN,Er_big,gamma,delta_rb_d2))                       
np.savetxt("BigScatterDaleyEst.csv",\
           np.vstack( \
                     (np.reshape(VoN,(1,len(VoN))),\
                      np.reshape(GscDE_big*2*np.pi,(1,len(GscDE_big))))\
                      ) \
           ,delimiter=",")
                     

#    wn1=np.transpose(wn1)
#    wn2=np.transpose(wn2)
#    wn3=np.transpose(wn3)
#
#    
#    bhU0s[cntr]=np.real(calBHU(wn0[::,0],xx)*Er)
#    bhE0s[cntr]=np.real(calBHE(wn0[::,0],Vext,xx)*Er)
#    
#    bhU1s[cntr]=np.real(calBHU(wn1[::,0],xx)*Er)
#    bhE1s[cntr]=np.real(calBHE(wn1[::,0],Vext,xx)*Er)
#    
#    bhU2s[cntr]=np.real(calBHU(wn2[::,0],xx)*Er)
#    bhE2s[cntr]=np.real(calBHE(wn2[::,0],Vext,xx)*Er)
#    
#    bhU3s[cntr]=np.real(calBHU(wn3[::,0],xx)*Er)
#    bhE3s[cntr]=np.real(calBHE(wn3[::,0],Vext,xx)*Er)
#    
#    for ni in range(nn):
#        bhJ0s[ni,cntr]=(np.real(calBHJ(wn0[::,0],wn0[::,ni+1],Vext,xx))*Er)
#
#        bhJ1s[ni,cntr]=(np.real(calBHJ(wn1[::,0],wn1[::,ni+1],Vext,xx))*Er)
#        bhJ2s[ni,cntr]=(np.real(calBHJ(wn2[::,0],wn2[::,ni+1],Vext,xx))*Er)
#        bhJ3s[ni,cntr]=(np.real(calBHJ(wn3[::,0],wn3[::,ni+1],Vext,xx))*Er)
#        
#    cntr+=1
#    
#    
#    wn0Sv=np.hstack((wn0Sv,np.reshape(wn0[::,0],xx.shape)))
#    wn1Sv=np.hstack((wn1Sv,np.reshape(wn1[::,0],xx.shape)))
#    wn2Sv=np.hstack((wn2Sv,np.reshape(wn2[::,0],xx.shape)))
#    wn3Sv=np.hstack((wn3Sv,np.reshape(wn3[::,0],xx.shape)))

#plt.plot(xx/2/np.pi,2.5*Vp/Vo+1)
#plt.plot(xx/2/np.pi,wn0[::,0]/3.0)
#plt.plot(xx/2/np.pi,wn1[::,0]/3.0+2*1.0/3)
#plt.plot(xx/2/np.pi,wn2[::,0]/3.0+2*2.0/3)
#plt.plot(xx/2/np.pi,wn3[::,0]/3.0+2*3.0/3)
##plt.plot(xx/2/np.pi,xx*0)
#plt.xlim([-3.5,3.5])
#plt.show()
#
#
#np.savetxt("LattOut.csv",VoN,delimiter=",")
#np.savetxt("W0.csv",wn0Sv,delimiter=",")
#np.savetxt("W1.csv",wn1Sv,delimiter=",")
#np.savetxt("W2.csv",wn2Sv,delimiter=",")
#np.savetxt("W3.csv",wn3Sv,delimiter=",")
#
#
#np.savetxt("BHJ0.csv",\
#           np.vstack((np.reshape(VoN,(1,len(VoN))),bhJ0s)) \
#           ,delimiter=",")
#np.savetxt("BHJ1.csv",\
#           np.vstack((np.reshape(VoN,(1,len(VoN))),bhJ1s)) \
#           ,delimiter=",")
#np.savetxt("BHJ2.csv",\
#           np.vstack((np.reshape(VoN,(1,len(VoN))),bhJ2s)) \
#           ,delimiter=",")
#np.savetxt("BHJ3.csv",\
#           np.vstack((np.reshape(VoN,(1,len(VoN))),bhJ3s)) \
#           ,delimiter=",")
#
#np.savetxt("BHE0.csv",\
#           np.vstack((np.reshape(VoN,(1,len(VoN))),np.transpose(bhE0s))) \
#           ,delimiter=",")
#np.savetxt("BHE1.csv",\
#           np.vstack((np.reshape(VoN,(1,len(VoN))),np.transpose(bhE1s))) \
#           ,delimiter=",")
#np.savetxt("BHE2.csv",\
#           np.vstack((np.reshape(VoN,(1,len(VoN))),np.transpose(bhE2s))) \
#           ,delimiter=",")
#np.savetxt("BHE3.csv",\
#           np.vstack((np.reshape(VoN,(1,len(VoN))),np.transpose(bhE3s))) \
#           ,delimiter=",")
#
#np.savetxt("BHU0.csv",\
#           np.vstack((np.reshape(VoN,(1,len(VoN))),np.transpose(bhU0s))) \
#           ,delimiter=",")
#np.savetxt("BHU1.csv",\
#           np.vstack((np.reshape(VoN,(1,len(VoN))),np.transpose(bhU1s))) \
#           ,delimiter=",")
#np.savetxt("BHU2.csv",\
#           np.vstack((np.reshape(VoN,(1,len(VoN))),np.transpose(bhU2s))) \
#           ,delimiter=",")
#np.savetxt("BHU3.csv",\
#           np.vstack((np.reshape(VoN,(1,len(VoN))),np.transpose(bhU3s))) \
#           ,delimiter=",")
#
#
#plt.plot(VoN,bhJ0s[0,::])
#plt.plot(VoN,bhJ1s[0,::])
#plt.plot(VoN,bhJ2s[0,::])
#plt.plot(VoN,bhJ3s[0,::])
#plt.show()
#
#plt.plot(VoN,bhE0s)
#plt.plot(VoN,bhE1s)
#plt.plot(VoN,bhE2s)
#plt.plot(VoN,bhE3s)
#plt.show()          