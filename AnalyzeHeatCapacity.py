# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 10:47:37 2016

@author: timge
"""

import pickle
import matplotlib.pyplot as plt
import numpy as np
import os

from DataAnalysis import *
from SIConversion import *


def LoadCvData (szFileName, bPlot=False):
    pkl_file = open(szFileName, 'rb')
    CvData = pickle.load(pkl_file)
    pkl_file.close()

    T = CvData['Temperature']
    Ekin = CvData['Ekin']
    N = CvData['Particles']
    rho = CvData['Rho']
    
    # Determine distribution of heat capacity by data blocking
    Tcor = GetCorrelation(Ekin, 0.05, bPlot)    
    Tblock = max(int(len(Ekin)/20), int(Tcor))  #Set minimum value for the correlation time
    M = math.trunc((len(Ekin)/Tblock))
    print('Number of blocks', M, 'with length', Tblock)
    CvList = []
    for i in range(M):
        Ekmean = np.mean(Ekin[i*Tblock:Tblock-1+i*Tblock])
        dEk = Ekin[i*Tblock:Tblock-1+i*Tblock]-Ekmean
        dEkmean = np.mean(dEk*dEk)
        CvList.append(3*Ekmean**2/(2*Ekmean**2-3*N*dEkmean))    
    
    Cv = np.mean(CvList)
    CvError = np.std(CvList)
    
    if bPlot:
        fig, axs = plt.subplots(nrows=1, ncols=1, sharex=True)
        ax = axs
        ax.hist(CvList, bins=50, normed=True)
        
    print('Cv = ', Cv, '+/-', CvError, 'at T = ', T)
    return T, Cv, CvError
    
def CvAnalysis (szFolder):
    T       = []
    Cv      = []
    CvError = []
    
    for subdir, dirs, files in os.walk(szFolder): #Loop through all files of the folder
        for file in files:
            filepath = subdir + os.sep + file

            if filepath.endswith(".out"):
                print('----------------------------------------------------------->')
                print ('Loaded file:', filepath)
                Ti, Cvi, CvErrori = LoadCvData(filepath)
                T.append(Ti*temp2K)
                Cv.append(Cvi)
                CvError.append(CvErrori)
                
    # Plotting of the heat capacity data
    print(temp2K)
    Tmax = 500
    fig, ax = plt.subplots(nrows=1, ncols=1, sharex=True)
    ax.errorbar(T, Cv, yerr=CvError, color='black',  fmt='D')
    ax.plot([0, Tmax], [1.5, 1.5], 'gold') #Expectation for gas phase (3/2)
    ax.plot([0, Tmax], [3.0, 3.0], 'indigo') #Expectation for liquid phase (3)
    ax.set_title('Heat Capacity')
    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel('Heat Capacity (kB)')
    ax.set_xlim([0, Tmax])
    ax.set_ylim([0, 6])

if __name__ == '__main__':
    CvAnalysis('Output/CV/')