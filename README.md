# Computational-Physics-Argon

import numpy as np
import math
import time

from mpl_toolkits.mplot3d import Axes3D
from matplotlib.pyplot import cm
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d.art3d import juggle_axes
from matplotlib.widgets import Button
import pickle

from numba import jit

class ArgonGas:
        ##--Initialization ----------------------------------------------------------------------
    def __init__(self, N, rho, T, Ensemble):    #Initialize with number of particles N, density Rho and initial Temperature, for a type of Ensemble
        self.N = N
        self.rho = rho
        self.L = (self.N/self.rho)**(1/3.0)
        print(self.L)
        self.T0 = T
        self.T = T
        
        self.Ensemble = Ensemble
        print(self.Ensemble)
        self.dt = 0.004
        self.bins = 200        #bin size for plotting of g(r)
        self.Lmax = self.L/2
        
        self.Ekin = []
        self.Vtot = []
        self.Etot = []  
        self.time = []
        
        self.iStart = 0
        
        self.Virial = []
        
        self.Gr = np.zeros((self.bins,0))
        
        self.Tequil = 0
        
        self.InitializePosVel()
        
        self.fig = plt.figure()
        self.axBox = self.fig.add_subplot(121,projection = '3d')
        self.axEn  = self.fig.add_subplot(122)
        
                
        self.axBExp = plt.axes([0.85, 0.01, 0.1, 0.075])        
        self.bExp = Button(self.axBExp, 'Export')        
        self.bExp.on_clicked(self.Export)
               
        return
    
    def InitializePosVel(self):  
        Nc = int((self.N/4)**(1/3))   #Calculate number of unit cells in each direction
        print(Nc)
        
        #Initialize the initial positions based of fcc stacking        
        self.r = np.zeros((3, self.N) )            
        n=0
        for i in range(Nc):
            for j in range(Nc):
                for k in range(Nc):
                    self.r[0, n] = i*self.L/Nc
                    self.r[1, n] = j*self.L/Nc
                    self.r[2, n] = k*self.L/Nc
                    n +=1
                    self.r[0, n] = i*self.L/Nc 
                    self.r[1, n] = j*self.L/Nc + 0.5*self.L/Nc
                    self.r[2, n] = k*self.L/Nc + 0.5*self.L/Nc
                    n +=1
                    self.r[0, n] = i*self.L/Nc + 0.5*self.L/Nc
                    self.r[1, n] = j*self.L/Nc 
                    self.r[2, n] = k*self.L/Nc + 0.5*self.L/Nc
                    n +=1
                    self.r[0, n] = i*self.L/Nc + 0.5*self.L/Nc
                    self.r[1, n] = j*self.L/Nc + 0.5*self.L/Nc
                    self.r[2, n] = k*self.L/Nc
                    n +=1
        self.r += self.L/Nc/4   #Center the particles (so they aren't at the boundary)
#        print('Initial positions')
#        print(self.r)
                    
        #Initialize the velocities based on a Maxwell distribution and random direction
        sigma = np.sqrt(self.T)
        self.v = np.random.normal(0, sigma, (3, self.N))
        for i in range(3):          #Set net velocity in each direction to 0
            self.v[i, :] -= sum(self.v[0, :])/self.N
#        print('Initial velocities')
#        print (self.v)
        return

    def PlotIC (self):
        self.fig2, self.ax = plt.subplots()
        self.ax2 = self.fig2.add_subplot(111, projection='3d')
        plt.quiver(self.r[0,:], self.r[1,:], self.r[2,:], self.v[0,:], self.v[1,:], self.v[2,:], length=self.L/10)
        self.scat2 = self.ax2.scatter(self.r[0,:], self.r[1,:], self.r[2,:], c='b', cmap='jet')
        return
               
        ##--Animating the particle box -------------------------------------------------------------
    def Animate(self):
        self.ani = animation.FuncAnimation(self.fig, self.Update, interval=1, init_func=self.SetupPlot, blit=True)
                                           
    def SetupPlot(self):
        c = ['b', 'r', 'g', 'y']
        self.scat = self.axBox.scatter(self.r[0,:], self.r[1,:], self.r[2,:],c=c, s=10, animated=True)
        self.axBox.set_xlim3d(0, self.L)
        self.axBox.set_ylim3d(0, self.L)
        self.axBox.set_zlim3d(0, self.L)
        
        self.axEn.set_ylabel('Energy')
        self.axEn.set_xlabel('Time')
        self.axEn.legend(['Kinetic', 'Potential', 'Total'])
        self.plEkin, = self.axEn.plot(self.time, self.Ekin, 'b-', label="Kinetic")
        self.plVtot, = self.axEn.plot(self.time, self.Vtot, 'r-', label="Potential")
        self.plEtot, = self.axEn.plot(self.time, self.Etot, 'm-', label="Total")
        self.plEkin.axes.set_xlim(0, 1)
        self.plEkin.axes.set_ylim(-10, 10)                       
        
        self.fig.canvas.mpl_connect('button_press_event', self.OnClick)
        
        return self.scat, self.plEkin, self.plVtot, self.plEtot
    
    def OnClick(self, event):
        if event.inaxes==self.axEn:
            self.iStart = int(event.xdata/self.dt)
            print('Starting measurement at', event.xdata, 's')
        return
    
    def Update(self, iter):
        next(self.DoTimeSteps())
        self.scat._offsets3d = ( np.ma.ravel(self.r[0,:]) , np.ma.ravel(self.r[1,:]) , np.ma.ravel(self.r[2,:]) )
                    
        self.time.append(iter*self.dt)
        if iter%25==0:              #Update the energy plot
            self.plEkin.set_data(self.time, self.Ekin)
            self.plVtot.set_data(self.time, self.Vtot)
            self.plEtot.set_data(self.time, self.Etot)  
            
            self.plEkin.axes.set_xlim(0, self.time[-1])
            self.plVtot.axes.set_xlim(0, self.time[-1])
            self.plEtot.axes.set_xlim(0, self.time[-1])
            
        
        if (iter-10)%100==0:        #Rescale y axis
            miny = min([min(self.Ekin), min(self.Vtot), min(self.Etot)])
            maxy = max([max(self.Ekin), max(self.Vtot), max(self.Etot)])
            self.plEkin.axes.set_ylim(miny, maxy)
            self.plVtot.axes.set_ylim(miny, maxy)
            self.plEtot.axes.set_ylim(miny, maxy)             
            
            r2, histo = self.GetGr()
            GrNew = np.zeros((self.bins,1))
            for i in range(len(GrNew)):
                GrNew[i,0]=histo[i]
            self.Gr = np.concatenate((self.Gr,GrNew), axis=1)
        
        plt.draw()
        return self.scat, self.plEkin, self.plVtot, self.plEtot
        
    def DoTimeSteps(self):
        F , V, Vir= GetFExtern(self.r, self.N, self.L)
        while True:                      
            F = self.DoTimeStep(F)
            yield self.r
            
    def DoTimeStep(self, F):     
        #Verlet algorithm
        self.v += F/2*self.dt
        self.r += self.v*self.dt 
        self.r = np.remainder(self.r, self.L)   #Return the particle into the box
#        t1 = time.clock()
        F, V, Vir = GetFExtern(self.r, self.N, self.L)
#        print(time.clock()-t1)
        self.v += F/2*self.dt
                             
        if self.Ensemble == 'Canonical':
            factor = math.sqrt(3*self.N*self.T0/sum(sum(self.v*self.v)))
            print('Factor', factor)
            self.v *= factor
            
        self.Ekin.append(0.5*sum(sum(self.v*self.v)))
        self.Vtot.append(V)
        self.Etot.append(self.Ekin[-1] + V)
        self.T = 2/3*self.Ekin[-1]/self.N
        print('Temperature', self.T, self.T0)
        self.Virial.append(Vir)
        print('Virial', self.Virial[-1])
#        print(self.Ekin[-1], self.Vtot[-1], self.Etot[-1])
#        print(np.sum(self.v, axis=1))
        return F
                
        ##--Calculcation of the Distances----------------------------------------------------
    def GetDistances(self, i):
        rvector = np.transpose(self.r[:,i]-np.transpose(self.r))
        rvector -= np.rint(rvector/self.L)*self.L        
        r2 = sum(rvector*rvector)
        return rvector, r2
       
        ##-- Perform a measurement----------------------------------------------------------------------
    def Export (self, event):
        print('Start export')
        np.savetxt('Output/Ekin.out', self.Ekin[self.iStart:-1], delimiter=',')
        np.savetxt('Output/Vtot.out', self.Vtot[self.iStart:-1], delimiter=',')
        np.savetxt('Output/Etot.out', self.Etot[self.iStart:-1], delimiter=',')
        print(self.Gr)
        
        #Export data for heat capacity
        CvDict = {'Temperature': self.T, 'Particles': self.N, 'Rho': self.rho, 'Ekin':self.Ekin[self.iStart:-1] }        
        output = open('Output/Cv.out', 'wb')
        pickle.dump(CvDict, output)
        output.close()
        
        #Export data for virial
        VirDict = {'Temperature': self.T, 'Particles': self.N, 'Rho': self.rho, 'Virial': self.Virial[self.iStart:-1]}
        output = open('Output/Virial.out', 'wb')
        pickle.dump(VirDict, output)
        output.close()
        
        #Export data for Structure factor
        r, hist = self.GetGr()
        GrDict = {'Temperature': self.T, 'Particles': self.N, 'Rho': self.rho, 'Range': r[0:-1], 'Structure': self.Gr}
        output = open('Output/Gr.out', 'wb')
        pickle.dump(GrDict, output)
        output.close()
    
        

        print('Export complete')
        return        
        
        ##-- Calculation of 'thermodynamic' Identities--------------------------------------------------------
    def GetGr(self):
        r2tot = np.zeros(self.N*(self.N-1))
        for j in range(self.N):
            rvector, r2 = self.GetDistances(j)
            r2 = np.delete(r2,j)
            r2tot[j*(self.N-1):(j+1)*(self.N-1)]= np.sqrt(r2)
        histo,r2 = np.histogram(r2tot,self.bins,(0,self.L/2))
        histo = histo/r2[0:-1]**2/(4*math.pi*self.L/2/self.bins)  #Normalize by volume of shell
        histo[0] = 0
        return r2, histo
                
    


@jit
def GetFExtern(r, N, L):
    F = np.zeros((3, N))
    V = 0
    Virial = 0
    for i in range(N):
        for j in range(i):
            dx = r[0,i]-r[0,j]
            dx = dx - np.rint(dx/L)*L
            dy = r[1,i]-r[1,j]
            dy = dy - np.rint(dy/L)*L
            dz = r[2,i]-r[2,j]
            dz = dz - np.rint(dz/L)*L
            dr2 = dx*dx+dy*dy+dz*dz
            dr2 = 1/dr2
            dr6 = dr2*dr2*dr2
                    
            V += dr6*(dr6-1)
            Fx = dr6*(2*dr6-1)*dr2*dx
            Fy = dr6*(2*dr6-1)*dr2*dy
            Fz = dr6*(2*dr6-1)*dr2*dz
            
            F[0,i] += Fx
            F[1,i] += Fy
            F[2,i] += Fz            
            F[0,j] -= Fx
            F[1,j] -= Fy
            F[2,j] -= Fz
            Virial += dr2*dr6*(24-48*dr6)*(dx*dx+dy*dy+dz*dz)
    V = V*4
    F = F*24
    return F, V, Virial

  
if __name__ == '__main__':
   MyArgonGas = ArgonGas(500, 1.4, 1.0, 'Canonical')
#   MyArgonGas.PlotIC()
   MyArgonGas.Animate()
   plt.show()
