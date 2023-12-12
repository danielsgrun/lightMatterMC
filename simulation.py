# New Monte-Carlo suite of codes #
# Includes definitions of ODT potential, RK4 time-evol.,
# AC-Stark shift, recoil from photon abs./re-emission etc. #
# as interfacing with the lightMatterMC Fotran 90 suite of codes
# Developed by Daniel S. Grun in 2023

# Daniel S. Grun, Innsbruck 2023

import numpy as np
from tqdm import tqdm
from joblib import Parallel, delayed 
import sys
from physical_definitions import *

import myLib # Fortran 90 suite of light-matter interaction codes

# arg_input = int(sys.argv[1])
# n_jobs=arg_input
n_jobs = 4
 
conversion = 0.1482e-24 * 1.113e-16 # conversion from a.u. to S.I.
alpha_GS = 430 * conversion # polarizability

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

def runSimulation_loss(P,T,w0,titf, s0=[0], lambd=[lambd583], Gamma=[Gamma583], absProj=[[1.0,0,0]], delta=[0], alpha=[alpha_GS], n_samples=int(1)):
    print(n_jobs)
    initialCond = generateInitialCond(P,T,w0,n_samples=n_samples,alpha=alpha_GS,lambd=lambd_trap)
    results = Parallel(n_jobs=n_jobs, backend="threading")(delayed(myLib.create_simulation_loss)(initialCond[i], P,T,w0,titf,s0,lambd,Gamma,absProj,delta, alpha)
                                               for i in tqdm(range(n_samples)))
    return results
        

def survivalProb(P,T,w0,titf,s0=[0],lambd=[lambd583],Gamma=[Gamma583],absProj=[[1.0,0,0]],delta=[0],alpha=[alpha_GS], n_samples=1):
    results = runSimulation_loss(P,T,w0,titf,s0,lambd,Gamma,absProj,delta,alpha,n_samples)
    results = np.array(results)
    surv = results[:,0]
    phScatt = results[:,2]
    survProb = len(surv[surv==1])/len(surv)
    return [survProb, phScatt]

if __name__ == "__main__":
    
    Gamma1 = Gamma583
    lambd1 = lambd583

    Gamma2 = Gamma841
    lambd2 = lambd841
    
    P = 2.4e-3 # trap power, in W
    T = 20e-6 # temperature, in K
    w0 = 1.0e-6 # trap waist, in um
    s0 = 0.8 # near-resonant field, saturation parameter

    alpha_E = 1*alpha_GS

    dParam1 = 1
    dParam2 = -80
    delta1 = -dParam1*180e3 # detuning from resonance, in Hz
    delta2 = -dParam2*8e3 # detuning from resonance, in Hz

    t0 = 0
    tf = 70e-3
    dt = 1e-8
    titf = [t0, tf, dt]
    
    n_samples = int(300)

    def checkSurvAbsZ(absZ):
        absProj = [[1.0, 0, absZ]]
        result = survivalProb(P,T,w0,titf,[s0],[lambd1],[Gamma1],absProj,[delta1], [alpha_GS], n_samples)
        return result
    
    def checkSurvPDelta(s0,delta1):
        absProj = [[1.0, 0, 0.18]]
        result = survivalProb(P,T,w0,titf,[s0],[lambd1],[Gamma1],absProj,[delta1], [alpha_GS], n_samples)
        return result

    absZscan = np.arange(0,0.4,0.05)
    s0Scan = np.arange(0.5,1.7,0.25)
    deltaScan = -np.arange(1,2.7,0.5) * 180e3

    print(n_jobs)
    
    # result = np.array([checkSurvAbsZ(absZscan[i]) 
    #             for i in tqdm(range(len(absZscan)))])

    result = np.array([[checkSurvPDelta(s0Scan[i], deltaScan[j]) 
                        for i in tqdm(range(len(s0Scan)))]
                       for j in tqdm(range(len(deltaScan)))])
    
    # direct = 'C:\\Users\\x2241135\\Desktop\\PhD\\codes\\new_monteCarlo\\results\\'
    # name = 'result_delta_m{0:1.0f}_sat{1:1.0f}_tf{2:1.0f}.txt'.format(dParam, s0, 1e3*tf)
    
    # np.save('resultEr_1horBeam_1vert841Beam_dm1-5_dmFinerScan_s0-5_s100_70ms.npy', result)
    
    # np.savetxt(direct+name, np.c_[abszAngle, result])
