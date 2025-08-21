import numpy as np

Gamma583 = 2*np.pi * 186e3 # Gamma of the 583 nm transition (in Hz)
Gamma841 = 2*np.pi * 8e3 # Gamma of the 841 nm transition (in Hz)
Gamma626 = 2*np.pi * 135e3 # Gamma of the 626 nm transition (Dy) (in Hz)
Gamma401 = 2*np.pi * 28e6 # Gamma of the 401 nm transition (in Hz)
Gamma631 = 2*np.pi * 28e3 # Gamma of the 631 nm transition (in Hz)

lambd401 = 401e-9
lambd583 = 583e-9
lambd841 = 841e-9
lambd626 = 626e-9
lambd631 = 631e-9

wavelengths = [401e-9, 583e-9, ]

lines = {'401':[lambd401, Gamma401],
         '583':[lambd583, Gamma583],
         '631':[lambd631, Gamma631],
         '626':[lambd626, Gamma626],
         '841':[lambd841, Gamma841]
         }


def normalizeProjections(absProj):
    absProj_norm = np.array([absProj[i]/np.sqrt(sum(abs(absProj[i])**2))
                             for i in range(len(absProj))])
    
    return absProj_norm



def nearResonantBeams(wavelengths, intensities, detunings, absProj):
    
    trueMask = (len(wavelengths) == len(intensities) == len(detunings) == len(absProj))
    
    if trueMask==False:
        
        print("Length of 'wavelengths', 'intensities', 'detunings' and 'absProj' should be the same.")
        exit 
        
    else:
        
        lambdas = []; gammas = []; s0 = []; deltas = []
        
        for i in range(len(wavelengths)):
            lam, gam = lines[wavelengths[i]]
            lambdas.append(lam)
            gammas.append(gam)
            deltas.append(detunings[i]*gam)
            s0.append(intensities[i])
            
        lambdas = np.array(lambdas)
        gammas = np.array(gammas)
        deltas = np.array(deltas)
        s0 = np.array(s0)
        
        absProj = normalizeProjections(absProj)
        
        
        
        return [lambdas, gammas, deltas, s0, absProj]