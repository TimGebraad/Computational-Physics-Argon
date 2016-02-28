# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 14:36:28 2016

@author: timge
"""
import math 
import numpy as np
import matplotlib.pyplot as plt

##Determination of the correlation function and the corresponding correlation time
def GetCorrelation (data, frac=0.2, bPlot=False):
    
    # Calculate correlation function
    kmax=int(frac*len(data))
    ck = np.zeros((kmax))
    Am = np.mean(data)
    for k in range (kmax):
        D = []
        for n in range (len(data)-k):
            D.append((data[n]-Am)*(data[n+k]-Am)) 
        ck[k]= np.mean(D)
        del D

    # Calculate correlation time
    Tcor = np.sum(ck)/ck[0]
    print ('Average', Am, ', Standard deviation', math.sqrt(np.var(data)), ', Data points', len(data))
    
    if bPlot:
        f, (ax1, ax2) = plt.subplots(2, 1)
        ax1.plot(data)
        ax2.plot(ck)
        
    return Tcor
    
    
##Determination of the mean value and error of the variable by the use of data blocking
def AnalyseQuantity (data, Tcor):
    Tblock = round(Tcor)
    M = math.trunc((len(data)/Tblock))
    print('Number of blocks', M, 'with length', Tblock)
    
    A = []
    for i in range(M):
        A.append(np.mean(data[i*Tblock:Tblock-1+i*Tblock]))
        
    Mp = np.mean(A)     
    error = np.std(A)

    return Mp, error
 
if __name__ == '__main__':
   Data = np.loadtxt('Output/Ekin - Equilibrium.out')
   Tcor = GetCorrelation(Data, 0.2, True)
   AnalyseQuantity(Data,Tcor)
