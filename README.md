ArgonGas.py
Contains the model for the argon gas structured in a class. The class consists of methods that can be divided roughly into three parts:
The initialisation in which the particles are positioned into a fcc structure and in which they obtain an initial velocity according to a Maxwell distribution.
The part where the system is iterated through time including the Verlet algorithm, the calculation of the force and the storage of data in arrays. This part also consists of animating the particles in the box and the different energies in time to be able to view the system reach equilibrium. 
The last part exports the data of the different energies and the data for the heat capacity, pressure and structure factor which can later be analyzed and displayed by other files.

SIConversion.py
This is a simple tool that generates the conversion factors between the units used in our model and S.I. units which can be used by the other modules to display the generated data in S.I. units.

DataAnalysis.py
This module has two functions that take care of the data analysis. The first one GetCorrelation takes as input a data string and calculates the correlation time based upon its autocorrelation. The other function AnalyseQuantity takes the data and the correlation time as input and estimates the mean value of the data including the error using data blocking, where the data blocks are the correlation time.

AnalyzeHeatCapacity.py
This module analyzes the heat capacity. It has two functions, the first one reads a file containing a dictionary of variables concerning the heat capacity and calculates the heat capacity itself based on data blocking. For this the heat capacity is calculated for each data block and from this distribution the heat capacity and the error are obtained.
The second functions takes a folder as input and runs over all files in that folder containing the heat capacity data, is let's the first function analyze the data and thereafter plots all the data in the folder in a single plot.

AnalyzePressure.py
This module analyzes the pressure. It has two functions, the first one reads a file containing a dictionary of variables concerning the pressure and calculates the pressure based on data blocking.
The second functions takes a folder as input and runs over all files in that folder containing the pressure data, is let's the first function analyze the data and thereafter plots all the data in the folder in a single plot.

AnalyzeStructureFactor.py
This module analyzes the structure factor. It has two functions, the first one reads a file containing a dictionary of variables concerning the pressure and calculates the structure factor including errors based on the distribution of the data.
The second functions takes a folder as input and runs over all files in that folder containing the structure factor data, is let's the first function analyze the data and thereafter plots all the data in the folder in a single plot.
