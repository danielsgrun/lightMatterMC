# New Monte-Carlo suite of codes #
# Includes definitions of ODT potential, RK4 time-evol.,
# AC-Stark shift, recoil from photon abs./re-emission etc. #
# as interfacing with the lightMatterMC Fotran 90 suite of codes


import numpy as np
from tqdm import tqdm
from joblib import Parallel, delayed 
import sys
from physical_definitions import *
from nearResonantLight import *
import os
import json

from scanFunctions import run_scan, getSurvivalPerAtom

C3_au_conversion = (electronCharge*a0)**2/(4*np.pi*e0) # C_3 conversion from a.u. to S.I.
alpha_conversion = 0.1482e-24 * 1.113e-16 # polarizability conversion from a.u. to S.I.

if __name__ == "__main__":
    
    # arg_input = int(sys.argv[1])
    # n_jobs=arg_input
    n_jobs = 5 

    spontCase = 's'   

    alpha_GS = 430 * conversion # polarizability

    nAtoms = 3
    
    C3u = 0.1 * C3_au_conversion # converting from a.u. to S.I.
    C3g = 0.1 * C3_au_conversion # converting from a.u. to S.I.

    C3Vals = np.array([C3u, -C3g]) * 1
    
    P = 1.6e-3 # trap power, in W
    T = 10e-6 # temperature, in K
    w0 = 0.8e-6 # trap waist, in um
    
    absProj = np.array([[1.0, 0.0, 0.18], [0,0,1.0]]) # fixed absorption projection    
    
    alpha_E = 1*alpha_GS # valid since we have a magic condition for the transition in the lab
    
    t0 = 0    
    tf = 30e-3
    dt = 1e-8
    
    titf = [t0, tf, dt]


    s0v = 0.1
    s0h = 0.5
    
    s0 = [s0h, s0v]
 
    n_samples = int(20) # we should do everything with the same # of samples. Choose either 3e2 or 1e3.
    
    nBeams = 2
    
    beamWavelengths = ['583', '583']
    beamDetunings = [-2, -2]
    
    lambdas, Gammas, deltas, s0, absProj = nearResonantBeams(beamWavelengths, s0, beamDetunings, absProj)
    
    modFreq = 1
    modulateAlpha = 0
    
    #%% 
    ## The following functions are defined as different scans (i.e. "different measurements") ##
        
    scan_vars = ['P']
    scan_values = [np.array([1,2,3])*1e-3]
    
    param_template = {
        'P': None,
        'C3Vals': C3Vals,
        'w0': w0,
        's0': s0,
        'deltas': deltas,
        'alphas': [alpha_GS]*nBeams,
        'lambd': lambdas,
        'Gammas': Gammas,
        'absProj': absProj,
        'nBeams': nBeams,
        'n_samples': n_samples,
        'nAtoms': nAtoms,
        'modFreq': modFreq,
        'modulateAlpha': modulateAlpha
    }
    
    static_params = {
        'T': T,
        'titf': titf
    }
    
    filename_template = "result_Pscan_w{w0}_s{s0}_d{deltas}_tf{titf}"
    scan_labels = ['Tweezer power (mW)']
    
    runner_args = {
        'n_jobs': 4,
        'spontCase': 's'
    }
    
    results, fileName, paramsDict = run_scan(
        scan_vars,
        scan_values,
        param_template,
        static_params,
        filename_template,
        scan_labels,
        runner_args
    )
    
    survivals, phScatt = getSurvivalPerAtom(results)


    
    
    
    # direct = 'C:\\Users\\c7041356\\Desktop\\LeonardoBG\\Science\\Project_TReqsSim\\code\\lightMatterMC_fortran90_v0102224\\results\\'
    
    # os.chdir('results')
    
    # with open(fileName+'.json', "w") as savingFile:
    #     json.dump(resultDict, savingFile)
    
    # with open(fileName+'_params.json', "w") as savingParamFile:
    #     json.dump(paramsDict, savingParamFile)
    
    # os.chdir('..')
    
#os.remove(myLib.cpython-310-x86_64-linux-gnu.so) 
