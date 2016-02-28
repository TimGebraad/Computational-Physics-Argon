# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 16:38:48 2016

@author: timge
"""
import  pickle
import matplotlib.pyplot as plt
import numpy as np

from DataAnalysis import *
from SIConversion import *


def LoadStructureFactor (szFileName):
    pkl_file = open(szFileName, 'rb')
    GrDict = pickle.load(pkl_file)
    pkl_file.close()

    r = GrDict['Range']
    T = GrDict['Temperature']
    Rho = GrDict['Rho']
    N = GrDict['Particles']
    Gr = GrDict['Structure']
    
    # Average the data in each point and calculate the standard deviation of the distribution
    GrErr = np.std(Gr, axis=1)
    Gr = np.mean(Gr, axis=1)
    
    print('Structure factor at T = ', T, 'and Density = ', Rho)
    
    # Convert into a correlation function
    Gr = Gr*2/Rho/(N-1)*N
    GrErr = GrErr*2/Rho/(N-1)*N
    return r, Gr, GrErr, Rho, T
    
def StructureAnalysis (szFolder):
    fig, ax = plt.subplots(nrows=1, ncols=1, sharex=True)
    Labels = []
    for subdir, dirs, files in os.walk(szFolder):     #Loop through all files of the folder
        for file in files:
            filepath = subdir + os.sep + file

            if filepath.endswith(".out"):
                print('----------------------------------------------------------->')
                print ('Loaded file:', filepath)
                r, Gr, GrErr, Rho, T = LoadStructureFactor(filepath)
                plt.errorbar(r*length2nm, Gr, yerr=GrErr, fmt='o')
                Labels.append(('Rho', int(Rho*density2gcm*100)/100, ', T', int(temp2K*T)))
        
    ax.set_title('Structure Factor')
    ax.set_ylabel('g(r)')
    ax.set_xlabel('r [nm]')
    ax.set_ylim([0, 11])
    ax.legend(Labels)
    return

    
if __name__ == '__main__':
    StructureAnalysis('Output/GR/')