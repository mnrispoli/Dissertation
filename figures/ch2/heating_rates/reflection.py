#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 13:08:58 2019

@author: matthewrispoli
"""

def rs(n1,n2,theta):
    #Fresnel coefficient for s-polarization
    #defined in field
    return abs((n1*np.cos(theta)-n2*np.sqrt(1-(n1*np.sin(theta)/n2)**2))/ \
            (n1*np.cos(theta)+n2*np.sqrt(1-(n1*np.sin(theta)/n2)**2)))

def cntrst(n1,n2,theta):
    return 4.0*rs(n1,n2,theta)/(2.0+2*rs(n1,n2,theta)**2)  

def offV(n1,n2,theta):       
    return (1.0+rs(n1,n2,theta)**2-2.0*rs(n1,n2,theta))/(4.0*rs(n1,n2,theta))

FX=6.5*2
FY=FX/3#/1.618

[XL,XH,YL,YH] = [0,0.25,0,4.0]
[XW,YT]=[XH-XL,YH-YL]

W=XW/2
H=FY
AR=XW/YT

fig = plt.figure(figsize=(FX, FY))

ax1 = fig.add_subplot(121, aspect='auto')
thetas=np.linspace(0,np.pi/2,1000)
#ax1.plot(thetas,np.abs(rs(1,1.45,thetas))**2,'-', label = 'wn')
ax1.plot(thetas,cntrst(1.0,1.45,thetas))
ax1.grid(True, which="both", ls="-")
ax1.set_title('Fresnel Reflection: S-Polarization, Intensity')
ax1.set_xlabel('Angle (rad)')
ax1.set_ylabel('Lattice Contrast')
ax1.set_xlim([0,np.pi/2])
ax1.set_ylim([-0.05,1.05])




ax2 = fig.add_subplot(122, aspect='auto')
ax2.plot(thetas,offV(1.0,1.45,thetas))
ax2.grid(True, which="both", ls="-")
ax2.set_title('Fresnel Reflection: S-Polarization, Intensity')
ax2.set_xlabel('Angle (rad)')
ax2.set_ylabel('Lattice Offset (V_o)')
ax2.set_xlim([0,np.pi/2])
ax2.set_ylim([-0.05,1.05])
#ax2.plot(thetas,offV(1.0,1.45,thetas))

#ax2.loglog(WnAx[0,::],WnAx[1,::],'-', label = 'wn')
#ax2.loglog(HOAx[0,::],HOAx[1,::],'-', label = '\psi_HO')
##ax2.loglog(DEAx[0,::],DEAx[1,::],'-')
#ax2.grid(True, which="both", ls="-")
#ax2.set_title('Axial Lattice')
#ax2.set_xlabel('Lattice Depth (Er)')
#ax2.set_ylabel('\Gamma_{sc}')
#ax2.set_ylim([1E-4,1E-1])
#ax2.set_xlim([1E-1,1E3])


plt.savefig('axLatRef.pdf')
