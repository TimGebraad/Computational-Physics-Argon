# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 11:33:47 2016

@author: timge
"""
import math

#Conversion factors to obtain the physical quantities in S.I. Units

# Used quantities in the model
eps = 1.6*10**(-21)           #[J]=[kg*m^2/s^2]   in S.I. vs 1 in model
sigma =3.4*10**(-10)          #[m]                in S.I. vs 1 in model
m = 6.63*10**(-26)            #[kg]               in S.I. vs 1 in model
kB = 1.38*10**(-23)           #[m^2*kg/s^2*K]     in S.I. vs 1 in model
T = eps/kB                    #[K]                in S.I. vs 1 in model (consequence of previous)
t = math.sqrt(sigma**2*m/eps) #[s]                in S.I. vs 1 in model (consequence of previous)


#Conversion factors
length2nm = sigma*10**9                         #Conversion of length to nm
density2kgm = m/sigma**3                        #Conversion of density to kg/m^3
density2gcm = density2kgm*10**3/10**6           #Conversion of density to g/cm^3
temp2K = T                                      #Conversion of temperature to Kelvin
pressure2bar = m/sigma/t**2*10**(-5)            #Conversion of pressure to bar
specvol2cmg = 1/density2gcm                     #Conversion of specific volume to cm^3/g