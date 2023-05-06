import numpy as np
import matplotlib.pyplot as plt
from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets as widgets
from IPython.display import clear_output
import plotly.graph_objects as go
from plotly.subplots import make_subplots

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
        rho = np.abs(psi)**2
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
        rho_x = self.calculate_rho(self.psi_x)
        
        fig = make_subplots(rows=3, cols=1,
                            subplot_titles=['Wellenfunktion im Ortsraum', 'Wellenfunktion im Impulsraum',
                                            'Aufenthaltswahrscheinlichkeit im Ortsraum'],
                            )
        
        for time in range(self.Nt):
            fig.add_trace(go.Scatter(x=self.x_grid, y=np.real(self.psi_x[time]), visible=False, 
                                     name=r'$\text{Re}\left(\psi\left(x,t\right)\right)$'), 1, 1)
            fig.add_trace(go.Scatter(x=self.x_grid, y=np.imag(self.psi_x[time]), visible=False, 
                                     name=r'$\text{Im}\left(\psi\left(x,t\right)\right)$'), 1, 1)
    
            fig.add_trace(go.Scatter(x=self.k_grid, y=np.real(self.psi_k[time]), visible=False, 
                                     name=r'$\text{Re}\left(\~\psi\left(k,t\right)\right)$'), 2, 1)
            fig.add_trace(go.Scatter(x=self.k_grid, y=np.imag(self.psi_k[time]), visible=False, 
                                     name=r'$\text{Im}\left(\~\psi\left(k,t\right)\right)$'), 2, 1)
            
            
            fig.add_trace(go.Scatter(x=self.x_grid, y=rho_x[time], visible=False, 
                                     name=r'$|\psi(x,t)|^2$'), 3, 1)
    
        fig.data[0].visible = True
        fig.data[1].visible = True
        fig.data[2].visible = True
        fig.data[3].visible = True
        fig.data[4].visible = True

        if potential_type=='Barriere' or potential_type=='Glauber_Zustand': #Plot potential
            fig.add_trace(go.Scatter(x=self.x_grid, y=self.v_vector, visible=True, 
                                     name=r'$V(x)$'), 3, 1)      
        
        steps = []
        for i in range(0, 5*self.Nt, 5):
            if potential_type=='Barriere' or potential_type=='Glauber_Zustand': # Potential always visible
                step = dict(method="update", args=[{"visible": [False] * (5*self.Nt+1)}], label=str(i//5))
                step["args"][0]["visible"][i:i+5] = [True,True,True,True,True]
                step["args"][0]["visible"][-1] = True 
            else:
                step = dict(method="update", args=[{"visible": [False] * 5*self.Nt}], label=str(i//5))
                step["args"][0]["visible"][i:i+5] = [True,True,True,True,True]
            steps.append(step)
            
        slider = [dict(active=0, steps=steps, currentvalue={"prefix": "t="})]

        fig.update_layout(sliders=slider, width=1000, height=800)
        
        if potential_type=='Freies_Teilchen':
            fig.update_xaxes(range=[-5, 15])
            fig.update_yaxes(range=[-1, 1])
            
        elif potential_type=='Barriere':
            fig.update_xaxes(range=[-70, 70], row=1, col=1)
            fig.update_xaxes(range=[-10, 10], row=2, col=1)
            fig.update_xaxes(range=[-70, 70], row=3, col=1)
            fig.update_yaxes(range=[-1, 1])
            
        elif potential_type=='Glauber_Zustand':
            fig.update_xaxes(range=[-5, 5])
            fig.update_yaxes(range=[-1, 1])
        
        fig.update_yaxes(range=[0, 1], row=3, col=1)
        
        fig.show()
  
        
    def solve(self, potential_type='Freies_Teilchen' ,v0=0, tmin=None, tmax=None, Nt=None, xmin=None, xmax=None, Nx=None, k0=None, d=None):
        # Check potential_type Input!!!
        
        if tmin==None or tmax==None or Nt==None:
            if potential_type=='Freies_Teilchen':
                tmin = 0
                tmax = 5
                Nt = 50
            elif potential_type=='Barriere':
                tmin = 0
                tmax = 100
                Nt = 100
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
        

#######################################################################        

        
class leiter_operatoren():
    def __init__(self, hbar=1, m=1, x_min=-10, x_max=10, Nx=2**9):
        self.hbar = hbar
        self.m = m
    
        self.x_min = x_min
        self.x_max = x_max
        self.Nx = Nx # Less grid points increases stability after several a/a^dagger operations
    
        self.dx = (x_max-x_min)/(Nx)#Nx-1???
        self.x_grid = np.arange(self.x_min,self.x_max,self.dx)
    
    def init_psi(self, x0=1):
        self.psi = np.zeros(self.Nx,dtype=complex)
        self.psi[:] = np.exp( -(self.x_grid/x0)**2 / 2 ) * np.sqrt(1/np.pi*x0) # not 1./np.sqrt(x0*np.sqrt(np.pi))?
    
    def diff_psi(self, psi):
        d_psi = np.zeros(self.Nx,dtype=complex)
        d_psi[0] = (psi[1] - psi[0]) / self.dx
        d_psi[-1] = (psi[-1] - psi[-2]) / self.dx
        for i in range(1,self.Nx-1):
            d_psi[i] = (psi[i+1] - psi[i-1]) / (2 * self.dx)
        return d_psi
    
    def apply_a(self):
        d_psi = self.diff_psi(self.psi)
        self.psi = (d_psi + self.x_grid * self.psi) / np.sqrt(2)
    
    def apply_a_dagger(self):
        d_psi = self.diff_psi(self.psi)
        self.psi = (-d_psi + self.x_grid * self.psi) / np.sqrt(2)
    
    def calc_rho(self,psi):
        return np.abs(psi)**2
    
    def plot_psi(self):
        rho = self.calc_rho(self.psi)
        
        fig, axs = plt.subplots(2)
        axs[0].plot(self.x_grid, np.real(self.psi[:]), "red", label='Re[$\psi(x,t)$]')
        axs[0].plot(self.x_grid, np.imag(self.psi[:]), "green", label='Im[$\psi(x,t)$]')

        axs[1].plot(self.x_grid, rho[:], "black", label='$|\psi(x,t)|^2$')
    
        axs[0].legend(loc='upper right')
        axs[1].legend(loc='upper right')
        fig.set_size_inches(10, 9)
    
        axs[0].set_xlim([self.x_min, self.x_max])
        axs[1].set_xlim([self.x_min, self.x_max])

        plt.rcParams['font.size'] = '20'
        plt.show()
    
    def display(self):
        self.init_psi()
        
        a_button = widgets.Button(description='Absteiger')
        a_dagger_button = widgets.Button(description='Aufsteiger')
        reset_button = widgets.Button(description='Reset')
        
        def a_button_function(_):
            with out:
                clear_output()
                self.apply_a()
                self.plot_psi()
                
        a_button.on_click(a_button_function)
        
        def a_dagger_button_function(_):
            with out:
                clear_output()
                self.apply_a_dagger()
                self.plot_psi()
                
        a_dagger_button.on_click(a_dagger_button_function)
        
        def reset_button_function(_):
            with out:
                clear_output()
                self.init_psi()
                self.plot_psi()
        
        reset_button.on_click(reset_button_function)
        
        out = widgets.Output()
        reset_button_function('')
        
        display(widgets.HBox([a_button,a_dagger_button,reset_button] ),out)
        