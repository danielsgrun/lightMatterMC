import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

save_files = False

nPoints = int(1e5)

xIso = np.zeros((nPoints))
yIso = np.zeros_like(xIso)
zIso = np.zeros_like(xIso)

xSigma = np.zeros((nPoints))
ySigma = np.zeros_like(xSigma)
zSigma = np.zeros_like(xSigma)

counter = 0

def Fiso(u):
    return 0.5*(1-np.cos(u))

def FSigma(u):
    return 1/8 * (4-3*np.cos(u)-np.cos(u)**3)

u = np.linspace(0,np.pi,1000)

interpIso = interp1d(Fiso(u), u)
interpSigma = interp1d(FSigma(u), u)

for i in range(len(xIso)):
    phi = np.random.rand() * 2*np.pi
    thetaIso = interpIso(np.random.rand())
    thetaSigma = interpSigma(np.random.rand())
    
    xIso[i] = np.cos(phi)*np.sin(thetaIso)
    yIso[i] = np.sin(phi)*np.sin(thetaIso)
    zIso[i] = np.cos(thetaIso)
    
    xSigma[i] = np.cos(phi)*np.sin(thetaSigma)
    ySigma[i] = np.sin(phi)*np.sin(thetaSigma)
    zSigma[i] = np.cos(thetaSigma)

# theta = np.random.rand(nPoints) * np.pi
# phi = 2*np.pi * np.random.rand(nPoints)
# x = np.cos(phi)*np.sin(theta)
# y = np.sin(phi)*np.sin(theta)
# z = np.cos(theta)

fs = 13
nBins = 50

plt.figure()
thetaIso = np.arccos(zIso)
freqsIso, binsIso = np.histogram(thetaIso, bins=nBins)
freqsIso = freqsIso / max(freqsIso)
binsIso = [0.5*(binsIso[i]+binsIso[i+1]) for i in range(nBins)]
plt.bar(binsIso, freqsIso, 1.0*(binsIso[1]-binsIso[0]), color='royalblue', ec='k')
plt.plot(binsIso, np.sin(binsIso), lw=2, alpha=0.6, color='r')
plt.xlabel("$\\theta$", fontsize=fs)
plt.ylabel("Normalized frequency (isotropic)", fontsize=fs)


plt.figure()
curveMax = 4/3 * np.sqrt(2/3)
thetaSigma = np.arccos(zSigma)
freqsSigma, binsSigma = np.histogram(thetaSigma, bins=nBins)
freqsSigma = freqsSigma / max(freqsSigma)
binsSigma = [0.5*(binsSigma[i]+binsSigma[i+1]) for i in range(nBins)]
plt.bar(binsSigma, freqsSigma, 1.0*(binsSigma[1]-binsSigma[0]), color='royalblue', ec='k')
plt.plot(binsSigma, 1/curveMax*np.sin(binsSigma)*(1+np.cos(binsSigma)**2), lw=2, alpha=0.6, color='r')
plt.xlabel("$\\theta$", fontsize=fs)
plt.ylabel("Normalized frequency ($\sigma$ emission)", fontsize=fs)

if save_files:
    np.savetxt('spontEmission_isotropic.txt', np.c_[xIso, yIso, zIso])
    np.savetxt('spontEmission_Sigma.txt', np.c_[xSigma, ySigma, zSigma])




