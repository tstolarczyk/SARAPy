# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 13:50:39 2015

@author: Stolar
"""
import numpy as np
from astropy import units as u
import matplotlib.pyplot as plt
import Atmosphere as atm

###############################################################################
# Parametrisation from M.Kole et al., Astropartcile Physics 62(2015)230-240
# Valid for standard geomagnetic conditions as thos of January 1st 2000
# Uncertainties given for slopes are been ignored, they can reach 10%
# h : atmospheric pressure in hPa
# S : Solar activity parameter : 0 at minimum (phi=250MV) and 1 at maximum (phi=1109MV), linear correlation between S and Phi
# lmdb : magnetic latitude in degrees

def NeutronFlux(Energy, height,S=0.5,latitude=40.0*u.degree):
 
# Convert to correct units and get values   
    E = (Energy.to(u.MeV)).value
    h = (height.to(u.hPa)).value
    lmdb = (latitude.to(u.degree)).value


    a =  6.e-4   + (1.85-1.35*S)*1.e-2*(1.-np.tanh(np.radians(180. -3.5*lmdb)))
    b =  1.1e-2  + (1.2 - 0.8*S)*1.e-1*(1.-np.tanh(np.radians(180. -3.5*lmdb)))
    c =  150.                     -33.*(1.-np.tanh(np.radians(180. -5.5*lmdb)))
    d = -4.0e-3  + (2.4 - 1.0*S)*1.e-2*(1.-np.tanh(np.radians(180. -4.4*lmdb)))
        
    alpha =         - 0.281*np.exp(-h/4.6)  + 0.732
    beta  =         - 0.186*np.exp(-h/13.2) + 1.308
    gamma = 0.011*h + 0.300*np.exp(-h/68.0) + 0.26
    delta =           0.660*np.exp(-h/8.5)  + 1.40

    A = (a*h + b)*np.exp(-h/c) + d
    B = A*0.9**(-alpha+beta)
    C = B*15**(-beta+gamma)
    D = C*70**(-gamma+delta)
    
    Flx = np.where((E>0.008) * (E<= 0.9), 
             (A*E**-alpha),
             np.where((E>0.9) * (E<=15), 
                      (B*E**-beta),
                      np.where((E>15) * (E<= 70), 
                               (C*E**-gamma),
                               np.where((E>70) * (E<= 1000),
                                        (D*E**-delta),-1.))))
    
     
    return Flx/(u.cm*u.cm*u.s*u.MeV)
###############################################################################
def NeutronFluxDisplay(Energy,latitude=40.0*u.degree):        

    hpressure = np.logspace(-7,3,50)*u.hPa
    print(hpressure)
    print(NeutronFlux(Energy, hpressure))
    fig=plt.figure(figsize=(8,4))
    a = fig.add_subplot(1,1,1)    
    a.set_xscale('log')
    plt.ylabel("A") #y labels
    plt.xlabel("Atmopsheric pressure") #x labels
    plt.plot(hpressure,NeutronFlux(Energy,hpressure,latitude=30.*u.degree),color='blue')
    plt.plot(hpressure,NeutronFlux(Energy,hpressure,latitude=40.*u.degree),color='red')
    plt.plot(hpressure,NeutronFlux(Energy,hpressure,latitude=50.*u.degree),color='green')
    plt.plot(hpressure,NeutronFlux(Energy,hpressure,latitude=60.*u.degree),color='black')

    alpha =                 - 0.281*np.exp(-hpressure/4.6)  + 0.732
    beta  =                 - 0.186*np.exp(-hpressure/13.2) + 1.308
    gamma = 0.011*hpressure + 0.300*np.exp(-hpressure/68.0) + 0.26
    delta =                   0.660*np.exp(-hpressure/8.5)  + 1.40
# Display alpha, beta, gamma, delta as in the paper        
    h=np.logspace(-0.5,3,25)
    fig=plt.figure(figsize=(8,4))
    plt.ylabel("alpha, beta, gamma, delta") #y labels
    plt.xlabel("Atmopsheric pressure") #x labels
    
    a1= fig.add_subplot(1,1,1)    
    a1.set_xscale('log')
    plt.plot(h,alpha,color='blue')
    plt.plot(h,beta,color='red')
    plt.plot(h,gamma,color='green')
    plt.plot(h,delta,color='black') 

##########################       
if __name__=="__main__":
    
#    E=1*u.MeV
#    NeutronFluxDisplay(E)
    
### Example of a neutron flux    
    E=1*u.MeV
    h = 25*u.hPa
    print("Neutron flux at E=",E," and h=",h," = ",NeutronFlux(E,h,S=0,latitude=62*u.degree) )

#    E=np.logspace(-2.,3.,100)
    E=np.linspace(0.1,20.,100)
    fig=plt.figure(figsize=(8,12))

    a1= fig.add_subplot(3,1,1)    
#    a1.set_xscale('log')
    a1.set_yscale('log')
    plt.ylabel("Neutronflux - Lat.=62 deg. - h=25hPa") #y labels
    plt.xlabel("Energy ") #x labels    
    plt.plot(E,NeutronFlux(E*u.MeV,25*u.hPa,S=0.,latitude=62*u.degree),color='blue')
    plt.plot(E,NeutronFlux(E*u.MeV,25*u.hPa,S=1.,latitude=62*u.degree),color='red')
 
    a2= fig.add_subplot(3,1,2)    
    a2.set_yscale('log')
    plt.ylabel("Neutronflux - Lat.=62 deg. - h=25hPa") #y labels
    plt.xlabel("Energy ") #x labels    
    plt.plot(E,NeutronFlux(E*u.MeV,25*u.hPa,S=0.,latitude=62*u.degree),color='blue')
    plt.plot(E,NeutronFlux(E*u.MeV,25*u.hPa,S=1.,latitude=62*u.degree),color='red')
   
    a3= fig.add_subplot(3,1,3)    
#    a3.set_xscale('log')
    a3.set_yscale('log')
    plt.ylabel("Neutronflux - - h=15hPa - Lat.=0,15,30,45,60,75,90 deg.") #y labels
    plt.xlabel("Energy ") #x labels    

    plt.plot(E,NeutronFlux(E*u.MeV,15*u.hPa,S=0.,latitude=0*u.degree),color='blue')
    plt.plot(E,NeutronFlux(E*u.MeV,15*u.hPa,S=0.,latitude=15*u.degree),color='blue')
    plt.plot(E,NeutronFlux(E*u.MeV,15*u.hPa,S=0.,latitude=30*u.degree),color='blue')
    plt.plot(E,NeutronFlux(E*u.MeV,15*u.hPa,S=0.,latitude=45*u.degree),color='blue')
    plt.plot(E,NeutronFlux(E*u.MeV,15*u.hPa,S=0.,latitude=60*u.degree),color='blue')
    plt.plot(E,NeutronFlux(E*u.MeV,15*u.hPa,S=0.,latitude=75*u.degree),color='blue')
    plt.plot(E,NeutronFlux(E*u.MeV,15*u.hPa,S=0.,latitude=90*u.degree),color='blue')


# Integrate neutron flux over energies and plot versus altitude
    print("=========================")
    EnMin = 0.01*u.MeV
    EnMax = 1000.*u.MeV
    print("Energy range = ", EnMin, EnMax)
  # Altitude - given by the photn flux calculation (below hmin, effect of the Earth neutron albedo
    hMax = 1000*u.hPa
    hMin = 10*u.hPa
    print("Altitude range = ", hMin, hMax)
# Integrate the neutron flux with energy 
    h0 = hMin
    flux_E = lambda Ex:NeutronFlux(Ex*u.MeV,h0,latitude=62*u.degree,S=0.8).value
    flux = scipy.integrate.quad(flux_E,EnMin.value,EnMax.value,epsrel=0.1)*(NeutronFlux(EnMin,hMin).unit*u.MeV)   # By default it will go to Joule 
    print(" - Neutron flux at h={0:8} : {1:.2e} +- {2:.2e} {3}".format(h0,flux[0].value,flux.value[1],flux.unit)    )
       
    
    
    

