# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 09:49:19 2016

@author: timge
"""

import pickle
import matplotlib.pyplot as plt
import numpy as np
import os

from DataAnalysis import *


def LoadPressureData (szFileName, bPlot=False):
    pkl_file = open(szFileName, 'rb')
    PresData = pickle.load(pkl_file)
    pkl_file.close()

    T = PresData['Temperature']
    Virial = PresData['Virial']
    N = PresData['Particles']
    rho = PresData['Rho']
    
    Tcor = GetCorrelation(Virial, 0.05, bPlot)
    Tcor = max(50, int(Tcor)) #Set minimum value of 50 for the correlation time
    
    VirM, VirError = AnalyseQuantity(Virial, Tcor) 
    print('Virial = ', VirM, '+/-', VirError) 

    Pres = rho*T-rho/(3*N)*VirM
    PresError = rho/(3*N)*VirError
    SpecVol = 1.0/rho
    
    print('Pressure = ', Pres, '+/-', PresError, 'at T = ', T, 'and Spec. Vol = ', SpecVol)
    return SpecVol, Pres, PresError, T
    
def PressureAnalysis (szFolder):    
    SpecVol   = {'Temp = 48 K': [],  'Temp = 84 K': [],  'Temp = 120 K': [], 'Temp = 300 K': []}
    Pres      = {'Temp = 48 K': [],  'Temp = 84 K': [],  'Temp = 120 K': [], 'Temp = 300 K': []}
    PresError = {'Temp = 48 K': [],  'Temp = 84 K': [],  'Temp = 120 K': [], 'Temp = 300 K': []}
    color     = {'Temp = 48 K': 'b', 'Temp = 84 K': 'k', 'Temp = 120 K': 'r', 'Temp = 300 K': 'm'}
    Temp      = {'Temp = 48 K': 0.4,  'Temp = 84 K': 0.7,  'Temp = 120 K': 1.0, 'Temp = 300 K': 2.5}
    Labels = []
    
    
    for subdir, dirs, files in os.walk(szFolder):   #Loop through all files of the folder
        for file in files:
            filepath = subdir + os.sep + file
            if filepath.endswith(".out"):
                print('----------------------------------------------------------->')
                print ('Loaded file:', filepath)
                SpecVoli, Presi, PresErrori, T = LoadPressureData(filepath, False)
                
                #Put the data in the appropriate Temperature regime
                for keys in Temp:
                    if Temp[keys]-0.1<T<Temp[keys]+0.1:
                        key = keys
                        break
                SpecVol[key].append(SpecVoli*specvol2cmg) 
                Pres[key].append(Presi*pressure2tentosixtNm)
                PresError[key].append(PresErrori*pressure2tentosixtNm)
                
    #Reorder the dictionaries to the specific volume
    for keyT in SpecVol:
        if not SpecVol[keyT] == []:
            order = [i[0] for i in sorted(enumerate(SpecVol[keyT]), key=lambda x:x[1])]
            SpecVol[keyT] = [ SpecVol[keyT][i] for i in order]
            Pres[keyT] = [ Pres[keyT][i] for i in order]
            PresError[keyT] = [ PresError[keyT][i] for i in order]

    #Plotting of the pressure data
    fig, ax = plt.subplots(nrows=1, ncols=1, sharex=True)
    RhoInvMax = 5.5
    for key in SpecVol:
        if not SpecVol[key] == []:   
            Labels.append(key)
#            ax.plot(np.linspace(0, RhoInvMax)*specvol2cmg, 1/np.linspace(0, RhoInvMax)*Temp[key]*pressure2tentosixtNm, color=color[key], alpha=0.2, linewidth=3, label=key)
            ax.errorbar(SpecVol[key], Pres[key], yerr=PresError[key], color=color[key], label=key)
    ax.set_title('Pressure')
    ax.set_xlabel('Specific volume (cm^3/g)')
    ax.set_ylabel('Pressure (bar)')    
    ax.legend(Labels)
#    ax.set_xlim([0, RhoInvMax]*SV2SI)
#    ax.set_ylim([-10, 25])
    return

if __name__ == '__main__':
    PressureAnalysis('Output/Pressure/')