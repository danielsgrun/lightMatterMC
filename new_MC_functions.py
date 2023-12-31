# New Monte-Carlo suite of codes #
# Includes definitions of ODT potential, RK4 time-evol.,
# AC-Stark shift, recoil from photon abs./re-emission etc. #
# Daniel S. Grun, Innsbruck 2023

import numpy as np
from tqdm import tqdm
from joblib import Parallel, delayed 
from photon_functions import *

amu = 1.66e-27 # Atomic ma87ss unit --> kg
m = 166 * amu # Erbium166 isotope mass in kg
# m = 87 * amu # Rb
kB = 1.38e-23 # Boltzmann constant, m^2 kg s^-2 K^-1
e0 = 8.854e-12 # electric permittivity, C^2 N^-1 m^-2
c = 299792458 # speed of light, m/s
g = 9.81 # gravitational acceleration, m/s^2
hbar = 1.05457e-34 # reduced Planck's constant

conversion = 0.1482e-24 * 1.113e-16 # conversion from a.u. to S.I.
alpha_GS = 430 * conversion # polarizability
lambd_trap = 488e-9 # laser wavelength, m

lambd583 = 583e-9
lambd841 = 841e-9

Gamma583 = 2*np.pi * 180e3
Gamma841 = 2*np.pi * 8e3

def generateInitialCond(P,T,w0,n_samples=int(1e3),alpha=alpha_GS,lambd=lambd_trap):
    u0 = P*alpha / (np.pi*c*e0*w0**2)
    zR = np.pi*w0**2/lambd

    T = 0.1 * u0/kB

    omega_perp = np.sqrt(4*u0/(m*w0**2)) # radial trap frequency, Hz
    omega_par = np.sqrt(2*u0/(m*zR**2)) # longitudinal trap frequency, Hz  
    
    dx_par = np.sqrt(kB*T/(m*omega_par**2))
    dx_perp = np.sqrt(kB*T/(m*omega_perp**2))
    dv = np.sqrt(kB*T/m)
      
    vz0 = np.random.normal(loc=0, scale=dv, size=n_samples)
    vy0 = np.random.normal(loc=0, scale=dv, size=n_samples)
    vx0 = np.random.normal(loc=0, scale=dv, size=n_samples)  
      
    [x0, y0, z0] = np.array([np.random.normal(loc=0, scale=1*dx_perp, size=n_samples),
                             np.random.normal(loc=0, scale=1*dx_perp, size=n_samples),
                             np.random.normal(loc=0, scale=1*dx_par, size=n_samples)])
    
    return np.c_[x0,y0,z0, vx0,vy0,vz0]



def trapPotential(r,P,w0,alpha=alpha_GS, lambd=lambd_trap):
    x,y,z = r
    U0 = P*alpha / (np.pi * e0 * c * w0**2)
    zR = np.pi * w0**2 / lambd
    U = -U0/(1+(z/zR)**2) * np.exp(-2/(w0**2) * (x**2 + y**2) / (1+(z/zR)**2))
    return U


def trapPotDerivs(r,U,w0,lambd=lambd_trap):
    x,y,z = r
    zR = np.pi*w0**2/lambd
    z_term = np.sqrt(1+(z/zR)**2)
    
    a_x = 4/m*x/(w0**2*z_term**2)*U
    a_y = 4/m*y/(w0**2*z_term**2)*U
    a_z = 2/m*z/(zR**4*w0**2*z_term**4)*(zR**2*(w0**2-2*(x**2+y**2)) + w0**2*z**2)*U
    
    return [a_x, a_y, a_z]

def kineticEnergy(v):
    return m/2 * sum(abs(v)**2)


def trapEvol(rvVector, P, w0):
    r = rvVector[:3]
    v = rvVector[3:]
    DrDt = [0]*6
    U = trapPotential(r,P,w0)

    ax,ay,az = trapPotDerivs(r,U,w0)
    
    for i in range(3):
        DrDt[i] = v[i]
    
    DrDt[3] = ax; DrDt[4] = ay; DrDt[5] = az - g
    
    return np.array(DrDt)


def odeRK45_solver(rvFunc, y0, dt, P, w0):
    k1 = dt * rvFunc(y0, P, w0)
    k2 = dt * rvFunc(y0+0.5*k1, P, w0)
    k3 = dt * rvFunc(y0+0.5*k2, P, w0)
    k4 = dt * rvFunc(y0+k3, P, w0)
    y = y0 + (k1 + 2*k2 + 2*k3 + k4)/6
        
    return np.array(y)

def getAcStarkShift(rvVector, P, w0, alpha_E):
    r = rvVector[:3]
    U_g = trapPotential(r, P, w0)
    U_e = trapPotential(r, P, w0, alpha=alpha_E)
    acStarkShift = -(U_e-U_g)/hbar
    return acStarkShift

def recoilVel(absProj, lambd, case):
    
    absProj = np.array(absProj)
    absProj /= np.sqrt(sum(abs(absProj)**2))
    
    spontDir = 1-2*np.random.random((3))
    spontDir /= np.sqrt(sum(abs(spontDir)**2))
    
    recoilDir = absProj*(case=='absorption')+spontDir*(case=='emission')
    
    k = 2*np.pi / lambd
    phRecoil = hbar*k/m * recoilDir

    return phRecoil



def solveMC(P,T,w0,times,n_samples=int(1e3), initCond='auto'):
    
    if initCond == 'auto':
        initialConds = generateInitialCond(P,T,w0)
    else:
        initialConds = initCond
    
    def prepareSolver(sample_index):
        i = sample_index
        sol = odeRK45_solver(trapEvol, initialConds[i], times, P, w0)
        return sol

    sols = Parallel(n_jobs=4)(delayed(prepareSolver)(i) 
                              for i in tqdm(range(n_samples)))
    sols = np.array(sols)
        
    return sols



if __name__ == "__main__":
    
    pass
