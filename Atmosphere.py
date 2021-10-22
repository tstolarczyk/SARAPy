# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 13:45:37 2015

@author: Stolar
"""
import numpy as np
from astropy import units as u
from astropy import constants as const
import matplotlib.pyplot as plt

###############################################################################
Mair   = 0.0289644*u.kg/u.mol     # Molar mass of dry air
Rair   = 287.058*u.J/u.kg/u.K     # Ideal (universal) gas constant
Lair   = 0.0065*u.K/u.m           # Temperature lapse rate
Rgaz   = 8.31447*u.J/u.mol/u.K    #
P0     = 101.325*u.kPa            # Sea level standard atmospheric pressure
T0     = 288.15*u.K               # Sea level standard temperature
gEarth = 9.80665*u.m/u.s/u.s      # earth-surface gravitationnal acceleration

nAir   = 0.02504e27/(u.m)**3      # Air atom number density -standard conditions
###############################################################################
# Air composition in VOLUME (wikipedia - density of air)
N2_Air = 0.7808
O2_Air = 0.2095
Ar_Air = 0.009340 
CO2_Air= 0.0003978

###############################################################################
def AirDensity(P,T): # From wikipedia - density of air
    Pressure = P.to(u.Pa)
    Temperature = T.to(u.K)       
    return Pressure/Temperature/Rair
    
def AirDensityFromKK(h): # Air density parametrisation - from K. Kosack - unused
    h = h.to(u.m)    
    return 0.00125*(u.g/(u.cm)**3)*np.exp(-h/(8000*u.m))    
    
def AltitudePressure(h): # From wikipedia - density of air
    h = h.to(u.m)
    alpha = gEarth*Mair/Rgaz/Lair
    return P0*(1-Lair*h/T0)**alpha

def AltitudeTemperature(h):
    h = h.to(u.m)
    return T0 - Lair*h
    
def AirAtomDensity(h):   
    h = h.to(u.m)
    return const.N_A*AltitudePressure(h)/Rgaz/AltitudeTemperature(h)
 
###############################################################################   
if __name__ == '__main__':  
    
    P = 100000.*u.Pa
    T = 273.15*u.K
    h = 30000.*u.m

    print " ** Air density at P =",P," and T=",T," is",AirDensity(P,T).si
    print
    Ptmp = AltitudePressure(h).si
    Ttmp = AltitudeTemperature(h)
    print " ** Air physical parameters at h=",h
    print "      - Pressure    :",Ptmp," or ",Ptmp.to(u.hPa)
    print "      - Temperature :",Ttmp
    print "      - Density     :",AirDensity(Ptmp,Ttmp).to(u.g/(u.cm)**3)," Karl:",AirDensityFromKK(h).to(u.g/(u.cm)**3)
    print "      - Mol. density:",AirAtomDensity(h).si
    print "           - N2     :",N2_Air*AirAtomDensity(h).si
    altitude=np.linspace(0.,100.,100)
    altitude = altitude*u.km
    
    fig=plt.figure(figsize=(8,4))
    plt.xlabel("altitude ("+str(altitude.unit)+")") #y labels
    plt.ylabel("Atmopsheric pressure (hPa)") #x labels
    
    a1= fig.add_subplot(1,1,1)    
    plt.plot(altitude,(AltitudePressure(altitude).to(u.hPa)).value,color='blue')
       
       
    
    