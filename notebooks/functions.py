import numpy as np
import matplotlib.pyplot as plt
from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets as widgets

class Schroedinger_Equation:
    def __init__(self, hbar=1, m=1):
        self.hbar = hbar
        self.m = m
        
        
    def initialize_grid(self, potential_type, tmin, tmax, Nt, xmin, xmax, Nx, k0=1, d=1):
        self.t_grid = np.linspace(tmin, tmax, Nt) #
        self.dt = self.t_grid[1] - self.t_grid[0]
        
        self.dx=(xmax-xmin)/(Nx)
        self.x_grid = np.arange(xmin,xmax,self.dx)  # Understand boundaries!!!
        
        self.dk = 2*np.pi/(xmax-xmin)
        self.k_grid = np.arange(-self.dk*Nx/2,self.dk*Nx/2,self.dk)
        
        self.Nt = Nt # Anzahl der Zeitschritte
        self.Nx = Nx
        self.Nk = Nx
        self.xmin = xmin
        self.xmax = xmax
        self.kmin = self.k_grid[0]
        self.kmax = self.k_grid[-1]
        
        self.fourier_setup()
        
        self.psi_x = np.zeros((self.Nt,self.Nx),dtype=complex)
        self.psi_k = np.zeros((self.Nt,self.Nx),dtype=complex)
        
        
        if potential_type == 'Glauber_Zustand':
            p=0.5  # Make variable
            q=0   
            self.psi_x[0] = 0.5 * np.exp(1j*self.x_grid*p/self.hbar) * np.exp(-0.5*(self.x_grid-q)**2)  # No Pi factor?
            self.psi_k[0] = self.discrete_fourier_trafo(self.psi_x[0])
        else:
            self.psi_k[0] = np.exp( -(self.k_grid-k0)**2 * d**2 ) # Vary d !!!
            self.psi_x[0] = self.inv_discrete_fourier_trafo(self.psi_k[0])
        
        
    def fourier_setup(self):
        self.discrete_fourier_trafo_setup = np.zeros((self.Nk,self.Nx),dtype=complex)
        for i in range(0,self.Nk):
            for j in range(0,self.Nx):
                self.discrete_fourier_trafo_setup[i,j] = self.dx / np.sqrt(2*np.pi) * np.exp(-1j*self.k_grid[i]*self.x_grid[j])
        
        self.inv_discrete_fourier_trafo_setup = np.zeros((self.Nx,self.Nk),dtype=complex)
        for j in range(0,self.Nx):
            for i in range(0,self.Nk):
                self.inv_discrete_fourier_trafo_setup[j,i] = self.dk / np.sqrt(2*np.pi) * np.exp(1j*self.k_grid[i]*self.x_grid[j])
        
        
    def discrete_fourier_trafo(self,psi_x):
        psi_k = self.discrete_fourier_trafo_setup @ psi_x
        return psi_k
    
    def inv_discrete_fourier_trafo(self,psi_k):
        psi_x = self.inv_discrete_fourier_trafo_setup @ psi_k
        return psi_x
    
    def calculate_rho(self,psi):
        rho = np.abs(psi)
        return rho
    
    def split_operator_method(self, potential_type = 'Freies_Teilchen', v=0):
        k_propagate_half = np.exp( -1j * (self.hbar/(2.*self.m)) * self.k_grid**2 * self.dt / 2.)
        k_propagate = k_propagate_half * k_propagate_half #freies Teilchen
        
        if potential_type == 'Freies_Teilchen':
            for t in range(1, self.Nt):
                self.psi_k[t] = self.psi_k[t-1] * k_propagate
                self.psi_x[t] = self.inv_discrete_fourier_trafo( self.psi_k[t] )
                
        elif potential_type == 'Barriere' or potential_type=='Glauber_Zustand':
            self.init_v(potential_type, v)
            x_propagate = np.exp(-1j*self.v_vector * self.dt / self.hbar)
            for t in range(1, self.Nt): #tqdm( range(1, self.Nt) ):
                psi_k_temporary = self.psi_k[t-1] * k_propagate_half
                psi_x_temporary = self.inv_discrete_fourier_trafo(psi_k_temporary) * x_propagate
                self.psi_k[t] = self.discrete_fourier_trafo(psi_x_temporary) * k_propagate_half
                self.psi_x[t] = self.inv_discrete_fourier_trafo(self.psi_k[t])
                
                
    def init_v(self, potential_type, v0, x0=40, b=2): #Gauss Potential
        if potential_type == 'Barriere':
            self.v_vector = v0*np.exp(-0.5*((self.x_grid-x0)/b)**2)
        elif potential_type == 'Glauber_Zustand':
            self.v_vector = 0.5 * self.x_grid**2
            
    
    def plot_psi(self, t, potential_type):
        fig, axs = plt.subplots(3)
        fig.set_size_inches(10, 12)
        #fig.suptitle('Freies Teilchen') # Variational title
        
        if potential_type=='Freies_Teilchen':
            psi_x_pos = 0
            psi_k_pos = 1
            rho_pos = 2
            
        elif potential_type=='Barriere' or potential_type=='Glauber_Zustand':
            psi_x_pos = 1
            psi_k_pos = 2
            rho_pos = 0
        
        
        self.line1, = axs[psi_x_pos].plot(self.x_grid, np.real(self.psi_x[t]), 'red', label='Re($\psi(x,t)$)')
        self.line2, = axs[psi_x_pos].plot(self.x_grid, np.imag(self.psi_x[t]), 'green', label='Im($\psi(x,t)$)')
        
        self.line3, = axs[psi_k_pos].plot(self.k_grid, np.real(self.psi_k[t]), 'red', label='Re($\~\psi(k,t)$)')
        self.line4, = axs[psi_k_pos].plot(self.k_grid, np.imag(self.psi_k[t]), 'green', label='Im($\~\psi(k,t)$)')
        
        rho_x = self.calculate_rho(self.psi_x)
        self.line5, = axs[rho_pos].plot(self.x_grid, rho_x[t], 'black', label='$|\psi(x,t)|^2$')
        
        
        if potential_type=='Barriere' or potential_type=='Glauber_Zustand':
            axs[rho_pos].plot(self.x_grid, self.v_vector, 'blue', label='$V(x)$')
        
        
        axs[psi_x_pos].set_xlim([self.xmin, self.xmax]) #Limits not ends of vectors!!!
        axs[psi_k_pos].set_xlim([self.kmin, self.kmax]) #Variable k-Limits
        axs[rho_pos].set_xlim([self.xmin, self.xmax])
        
        axs[psi_x_pos].set_ylim([-1, 1])
        axs[psi_k_pos].set_ylim([-1, 1])
        axs[rho_pos].set_ylim([0, 1])
        
        plt.rcParams['font.size'] = '20'
        
        axs[0].legend(loc='upper right')
        axs[1].legend(loc='upper right')
        axs[2].legend(loc='upper right')
        
        def change_time(time=0):
            self.line1.set_ydata(np.real(self.psi_x[time]))
            self.line2.set_ydata(np.imag(self.psi_x[time]))
            self.line3.set_ydata(np.real(self.psi_k[time]))
            self.line4.set_ydata(np.imag(self.psi_k[time]))
            self.line5.set_ydata(rho_x[time])
        
        interact(change_time, time = widgets.IntSlider(min=0, max=self.Nt-1, step=1, value=0))
        
        
    def solve(self, potential_type='Freies_Teilchen' ,v0=0, tmin=None, tmax=None, Nt=None, xmin=None, xmax=None, Nx=None, k0=None, d=None):
        # Check potential_type Input!!!
        
        if tmin==None or tmax==None or Nt==None:
            if potential_type=='Freies_Teilchen':
                tmin = 0
                tmax = 4
                Nt = 40
            elif potential_type=='Barriere':
                tmin = 0
                tmax = 200
                Nt = 400
            elif potential_type=='Glauber_Zustand':
                tmin = 0
                tmax = 20
                Nt = 100
        else:
            assert tmin < tmax, "Min<Max required"
            assert Nt>0, 'Nt>0 required'

        if xmin==None or xmax==None or Nx==None:
            if potential_type=='Freies_Teilchen':
                xmin = -30
                xmax = 30
                Nx = 2**10
            elif potential_type=='Barriere':
                xmin = -100
                xmax = 100
                Nx = 2**10
            elif potential_type=='Glauber_Zustand':
                xmin = -6
                xmax = 6
                Nx = 2**10
        else:
            assert xmin < xmax, "Min<Max required"
            assert Nx>0, 'Nt>0 required'
            
        if k0==None:
            if potential_type=='Freies_Teilchen':
                k0 = 1
            elif potential_type=='Barriere':
                k0 = 2
        if d==None:
            if potential_type=='Freies_Teilchen':
                d = 1
            elif potential_type=='Barriere':
                d = 2
        
        
        self.initialize_grid(potential_type, tmin, tmax, Nt, xmin, xmax, Nx, k0, d)
        
        self.split_operator_method(potential_type, v0)
        
        self.plot_psi(0, potential_type)
        
