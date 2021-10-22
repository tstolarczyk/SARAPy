# -*- coding: utf-8 -*-
"""
Éditeur de Spyder

Ceci est un script temporaire.
"""
import numpy as np
from astropy import units as u
import matplotlib.pyplot as plt
import scipy
from scipy import interpolate # not imported by default !
import scipy.integrate # not imported by default !
from scipy.interpolate import InterpolatedUnivariateSpline 
import Atmosphere as atm
import Spallation as nspal
from numba import jit
##########################
# Energies in GeV
Emin = 1.*u.GeV 
Emax = 1e6*u.GeV
nEbins=100
    
###############################################################################
# Energy in the file are in eV !
# Min energy =  40000 - > 40e-3
# Max energy = 2e7 eV -> 20 MeV    
###############################################################################    
def CrossSection_14C(Energy):
    En = Energy.to(u.eV)
    csdata = np.loadtxt('C14-crossSection-40keV.txt', skiprows=6,usecols=(0,1))
    E=csdata[:,0]
    CS=csdata[:,1]
# This seems to not work - maybe too many points ?
#    CS_14C=scipy.interpolate.InterpolatedUnivariateSpline(E,CS,k=order)
    CS_14C=scipy.interpolate.interp1d(E,CS)
    
    return CS_14C(En.value)*u.barn
    
###############################################################################
# Supernova spectrum
###############################################################################
class SNR:
    
    def __init__(self):

        self.ESN = 1e51*u.erg  # 10**51 erg -> GeV
        self.dist = 100*u.pc   # 100 pc - > cm
        self.EpsGam = 0.01               # Fraction of CR into gamma
        self.EpsCR = 0.1                 # Fraction of SN enéergy into CR
        self.dOmeg = 4.*scipy.pi       # Solid angle fraction
        self.Dt = 1000*u.year             # SNR duration in sec (1000 years)
        self.Emin = 0.1*u.GeV
        self.Emax = 1000.*u.GeV        
        
        self.norm = (self.EpsGam
                    *self.EpsCR
                    *self.ESN.to(u.GeV)
                    /(4*scipy.pi*self.dist.to(u.cm)**2
                      *self.dOmeg
                      *np.log(self.Emax/self.Emin))
                      /self.Dt.to(u.s))
        self.gamma = 2.                  # Spectral index
    
    # energy spectrum on Earth
    def dfdE(self,E):
        return self.norm*(E*u.GeV)**-self.gamma
        
    def __str__(self): # Return a string
        print "------------ SNR -----------------"
        print "Distance                        : ",self.dist   
        print "Total kinetic energy            : ",self.ESN
        print "  - Fraction to CR              : ",self.EpsCR  
        print "  - CR fraction to gamma        : ",self.EpsGam  
        print "  - Total energy to gamma       : ",self.EpsCR*self.EpsGam*self.ESN
        print
        print "Energy boundaries               : ",self.Emin,", ",self.Emax
        print "Spectral index                  : ",self.gamma
        print "  ->Corresponding normalisation : ",self.norm
        return "--------------------------------"
            
###############################################################################
# RC spectrum
###############################################################################
class RC:
    def __init__(self):
        self.norm= 3.01
        self.gamma= 2.68
        self.Emin = 1.0*u.GeV
        self.Emax = 1e15*u.eV
        
    def dfdE(self,E):
        return self.norm*(E**-self.gamma)/u.cm/u.cm/u.s/u.GeV
    
    def __str__(self): # Return a string
        print "------------ CR -----------------"
        print "Energy boundaries               : ",self.Emin,", ",self.Emax
        print "Spectral index                  : ",self.gamma
        print "Normalisation                   : ",self.norm
        return "--------------------------------"
        
###############################################################################
# PLot incoming spectra
###############################################################################        
def PlotSpectra():
#Generate energies    
    E = np.logspace(np.log10(Emin.value),np.log10(Emax.value),nEbins)*u.GeV # 100 linearly spaced numbers
#    print E
#    print snr.dfdE(E)
    
    fig = plt.figure(figsize=(15,4))     
    
    a1 = fig.add_subplot(1,3,1)
    a1.set_yscale('log')
    a1.set_xscale('log')
    plt.ylabel("log(dfdE) "+str(snr.dfdE(1).unit)) #y labels
    plt.xlabel("log(E)"+str(E.unit)) #x labels
    plt.plot(E.value,snr.dfdE(E).value,'r')
    
    a2 = fig.add_subplot(1,3,2)
    a2.set_yscale('log')
    a2.set_xscale('log')
    plt.ylabel("log(dfdE) "+str(rc.dfdE(1).unit)) #y labels
    plt.xlabel("log(E)"+str(E.unit)) #x labels
    plt.plot(E.value,rc.dfdE(E).value,'b')    
    
    a3 = fig.add_subplot(1,3,3)
    a3.set_yscale('log')
    a3.set_xscale('log')
    plt.ylabel("Ratio") #y labels
    plt.xlabel("log(E)"+str(E.unit)) #x labels
    plt.plot(E.value,rc.dfdE(E).value/snr.dfdE(E).value,'g')   
    
    fig=plt.figure(figsize=(8,4))
    a = fig.add_subplot(1,1,1)    
    a.set_yscale('log')
    a.set_xscale('log')
    plt.ylabel("log(dfdE) "+str(rc.dfdE(1).unit)) #y labels
    plt.xlabel("log(E)"+str(E.unit)) #x labels
    plt.plot(E.value,snr.dfdE(E).value,'r')
    plt.plot(E.value,rc.dfdE(E).value,'b')   
    
###############################################################################
# 14C production at E and h
# E is in MeV because of the cross section and NeutronFlux
# h is in meter in AirAtomDensity    
###############################################################################    
def  dn14dEdh(E,h):
    E = E.to(u.MeV)
    h = h.to(u.m)
    rate = atm.N2_Air*atm.AirAtomDensity(h) \
             * CrossSection_14C(E) \
             * nspal.NeutronFlux(E,atm.AltitudePressure(h),S=1.,latitude=62.0*u.degree) 
    return rate.to(1/(u.MeV*u.cm*u.cm*u.cm*u.s))
    
    
###############################################################################
# M A I N 
# Checks by hand
# Done for 10hPa (25.9 km) and 10 MeV
    
###############################################################################    
if __name__=="__main__":

### Incoming fluxes
#    snr=SNR()
#    print snr
#    rc=RC()
#    print rc
#    PlotSpectra() # PLot E spectra - crosschecks units

# Boundaries - limit of validity
# Energy - given by the 14C cross section data points
    EnMin = 0.1*u.MeV
    EnMax = 20.*u.MeV
  # Altitude - given by the photn flux calculation (below hmin, effect of the Earth neutron albedo
    hMin = 5000.*u.m
    hMax = 35000.*u.m
    
# Compute 14C production rate at a given altitude and neutron energy
#    h0= hMin
#    E0 = EnMax/2. 
    h0 = 25900*u.m
    E0 = 10*u.MeV
        
    print "======================================="
    print "14C production rate at :"
    print "  h= ",h0," (",atm.AltitudePressure(h0),")"
    print "  E0=",E0
    print "    Air density  = ",atm.N2_Air*atm.AirAtomDensity(h0).to(1/(u.cm)**3)
    print "    14C cross s. = ",CrossSection_14C(E0)," or ",CrossSection_14C(E0).si
    print "    n flux       = ",nspal.NeutronFlux(E0,atm.AltitudePressure(h0),S=0.,latitude=62.0*u.degree)
    print "         --->",dn14dEdh(E0,h0)
    print "======================================="    
    
   
# Integrate the 14C production with energy - use the integrant has a function - precise but can be very long
    dn14dhdE_E = lambda Ex:dn14dEdh(Ex*u.MeV,h0).value
    dn14dh = scipy.integrate.quad(dn14dhdE_E,EnMin.value,EnMax.value,epsrel=0.1)*(dn14dEdh(E0,h0).unit*u.MeV)   # By default it will go to Joule 
    print " - 14C production rate at h={0:8} : {1:.2e} +- {2:.2e} {3}".format(h0,dn14dh[0].value,dn14dh.value[1],dn14dh.unit)    
    
    # Alternative method from data points
    En = np.linspace(EnMin,EnMax,500) # Result identical if one goes to 200 to 1000, unsufficient if 120 bins
#    print "dn14dEdh =",dn14dEdh(En,h0) 
    dn14dh_1 = scipy.integrate.simps(dn14dEdh(En,h0),x=En)*(dn14dEdh(E0,h0).unit*u.MeV) # By default it will go to Joule   
    print "                    Check h={0:8} : {1:.2e}".format(h0,dn14dh_1)
#      
# Integrate the 14C production with altitude - - use the integrant has a function - precise but can be very long
    dn14dhdE_h = lambda hx:dn14dEdh(E0,hx*u.m).value
    dn14dE = scipy.integrate.quad(dn14dhdE_h,hMin.value,hMax.value,epsrel=0.1)*(dn14dEdh(E0,h0).unit*u.cm)
    print " - 14C production rate at E={0:6} : {1:.2e} +- {2:.2e} {3}".format(E0,dn14dE[0].value,dn14dE[1].value,dn14dE.unit)     

   # Alternative method 
    h = np.linspace(hMin,hMax,10)
#    print "h=",h," - dn14dEdh =",dn14dEdh(E0,h)     
    dn14dE_1 = scipy.integrate.simps(dn14dEdh(E0,h),x=h)*(dn14dEdh(E0,h0).unit*u.cm)
    print "                    Check E={0:6} : {1:.2e}".format(E0,dn14dE_1)    
#
# Plot 14C prodcution rate at all altitude versus energies - cannot be done with "quad": would be too long       
    fig = plt.figure(figsize=(8,4))
    a1= fig.add_subplot(1,1,1)    
    plt.ylabel("14C rate ("+str(dn14dE.unit)+")") #y labels
    plt.xlabel("Energy ("+str(En.unit)+")") #x labels    
    dn14dE_E = [scipy.integrate.simps(dn14dEdh(E,h),x=h) for E in En]
    plt.plot(En,dn14dE_E,"r")     

# Integrate the 14C prodcution rate with altitude and energy
# Wikipedia says the 14C production rate on Earth is around 22000 / m2/s   
    dn14dEdh_Eh = lambda x,y:dn14dEdh(x*u.MeV,y*u.m).value
    dn14 = scipy.integrate.dblquad(dn14dEdh_Eh,hMin.value,hMax.value,lambda x:EnMin.value,lambda x:EnMax.value,epsrel=0.1)*(dn14dEdh(E0,h0).unit*u.cm*u.MeV)
    print "Total 14C production rate in the atmosphere : {0:.2e} +- {1:.2e} {2}".format(dn14[0].value,dn14[1].value,dn14.unit)    
            
# Plot 14C production rate at a given altitude for all valid energy
#    fig=plt.figure(figsize=(12,8)) 
#    a1= fig.add_subplot(2,2,1)    
#    plt.ylabel("Cross section") #y labels
#    plt.xlabel("Energy ") #x labels 
#    plt.plot(En.value,CrossSection_14C(En).value,"r")  
#    
#    a3= fig.add_subplot(2,2,2)    
#    plt.ylabel("Neutron flux"+str(nspal.NeutronFlux(En,atm.AltitudePressure(h),S=0.5,latitude=40.0*u.degree).unit)) #y labels
#    plt.xlabel("Energy ") #x labels 
#    a3.set_yscale('log')
#    plt.plot(En,nspal.NeutronFlux(En,atm.AltitudePressure(h),S=0.5,latitude=40.0*u.degree),"r") 
#    
#    a4 = fig.add_subplot(2,2,3)   
#    plt.ylabel("14C rate at h="+str(h)) #y labels
#    plt.xlabel("Energy ") #x labels 
#    plt.plot(En,dn14dEdh.to(1/(u.MeV*u.m*u.m*u.m*u.s)),"r")      
    


# Produce the same result with h variation

#    dn14dEdh_Eh = atm.AirAtomDensity(h) \
#             * CrossSection_14C(En) \
#             * nspal.NeutronFlux(En,atm.AltitudePressure(h),S=0.5,latitude=40.0*u.degree)
# Recompute normalisation factors, check integrals
#    print "\n Cross checks"
#    Elin = np.logspace(np.log10(snr.Emin.value),np.log10(snr.Emax.value),nEbins)
#    Elin = Elin*u.GeV
#    print Elin
## SN total energy
#    ESNRtot = scipy.integrate.simps(Elin.value*snr.dfdE(Elin.value).value, \
#                                    Elin.value) \
#                                    *snr.dfdE(1).unit*Elin.unit
##                                    
#    print "Total SNR gamma E flux on Earth =",ESNRtot 
#    print "Total SNR gamma energy over SNR interval= ", (ESNRtot*snr.Dt.to(u.s)
#                                       *4*np.pi*snr.dist.to(u.cm)**2
#                                       *4*np.pi).to(u.erg)
##   
#                                       
#   ERCtot = scipy.integrate.simps(Elin.value*rc.dfdE(Elin.value).value, \
#                                   Elin.value) \
#                                   *rc.dfdE(1).unit*Elin.unit                                        
#   print "Total CR flux = ", ERCtot