import myLib  # Fortran 90 suite of light-matter interaction codes
from physical_definitions import *
from nearResonantLight import *
import numpy as np
from joblib import Parallel, delayed
from tqdm import tqdm

class simClass:
    def __init__(self, P, T, w0, titf, C3Vals, s0=[0], nBeams=1, lambd=[lambd583], 
                 Gammas=[Gamma583], absProj=[[1.0,0,0]], deltas=[0], alphas=[alpha_GS],
                 n_samples=int(1e3), nAtoms=3, n_jobs=1, spontCase='s', modFreq=1, modulateAlpha=0):
        
        '''
        
        Parameters
        ----------
        P : float
            Power of a single optical tweezer (in W).
        T : float
            Temperature of the ensemble of samples (in K).
        w0 : float
            Minimum waist of a single optical tweezer (in m).
        titf : list/array of float.
            [Initial time, final time, time-step] (in s).
        C3Vals : list/array of float.
            [Repulsive C3 (pos.), attractive C3 (neg.)].
        s0 : list/array of float, optional
            list/array containing saturation parameters for each near-resonant beam.
            [s0 (beam 1), s1 (beam 2), ... sn (beam n)]. The default is [0].
        nBeams : integer, optional
            Number of near-resonant beams. The default is 1.
        lambd : list/array of float, optional
            list/array containing the wavelength of each near-resonant beam (in m).
            [lambda1 (beam 1), lambda2 (beam 2), ... lambdan (beam n)]. The default is [lambd583].
        Gammas : list/array of float, optional
            list/array containing the Gamma of the transition caused by each near-resonant beam (in 1/s).
            [Gamma1 (beam 1), Gamma2 (beam 2), ... Gamman (beam n)]. The default is [Gamma583].
        absProj : list/array of lists/arrays, optional
            list/array containing the unitary projection of each near-resonant beam. 
            [[1,0,0] (beam 1), [0,0,1] (beam 2), ... [0,1,0] (beam n)]. The default is [[1.0,0,0]].
        deltas : list/array of float, optional
            list/array containing the detuning from resonance of each near-resonant beam (in 1/s).
            [Delta1 (beam 1), Delta2 (beam 2), ... Deltan (beam n)]. The default is [0].
        alphas : list/array of float, optional
            list/array containing the polarizability of the excited state corresp. to the a transition 
            from each near-resonant beam (in S.I.).
            [alpha1 (beam 1), alpha2 (beam 2), ... alphan (beam n)]. The default is [alpha_GS].
        n_samples : int, optional
            Number of samples for averaging. The default is int(1e3).
        nAtoms : int, optional
            Number of atoms per simulation run. The default is 3.
        n_jobs : int, optional
            Number of cores to be used when parallellizing. The default is 1.
        spontCase : str, optional
            Dipole emission pattern for spontaneous emissions to consider. The default is 's' ('sigma').
        modFreq : float, optional
            Frequency of modulation for the polarizabilities. The default is 1.
        modulateAlpha : int, optional
            Modulate the polarizabilities in time (1 or 0). The default is 0.

        Returns
        -------
        None.

        '''
        
        self.P = P
        self.T = T
        self.w0 = w0
        self.titf = titf
        self.C3Vals = C3Vals
        self.s0 = s0
        self.nBeams = nBeams
        self.lambd = lambd
        self.Gammas = Gammas
        self.absProj = absProj
        self.deltas = deltas
        self.alphas = alphas
        self.n_samples = n_samples
        self.nAtoms = nAtoms
        self.n_jobs = n_jobs
        self.spontCase = spontCase
        self.modFreq = modFreq
        self.modulateAlpha = modulateAlpha

    def _generateInitialCond(self):
        alpha = self.alphas[0]
        lambd = self.lambd[0]

        u0 = self.P * alpha / (np.pi * c * e0 * self.w0 ** 2)
        zR = np.pi * self.w0 ** 2 / lambd
        T = 0.1 * u0 / kB

        omega_perp = np.sqrt(4 * u0 / (m * self.w0 ** 2))
        omega_par = np.sqrt(2 * u0 / (m * zR ** 2))

        dx_par = np.sqrt(kB * T / (m * omega_par ** 2))
        dx_perp = np.sqrt(kB * T / (m * omega_perp ** 2))
        dv = np.sqrt(kB * T / m)

        vz0 = np.random.normal(0, dv, size=(self.n_samples, self.nAtoms))
        vy0 = np.random.normal(0, dv, size=(self.n_samples, self.nAtoms))
        vx0 = np.random.normal(0, dv, size=(self.n_samples, self.nAtoms))

        x0 = np.random.normal(0, dx_perp, size=(self.n_samples, self.nAtoms))
        y0 = np.random.normal(0, dx_perp, size=(self.n_samples, self.nAtoms))
        z0 = np.random.normal(0, dx_par, size=(self.n_samples, self.nAtoms))

        initialCond = np.array([x0, y0, z0, vx0, vy0, vz0], order='F')
        initialCond = np.moveaxis(initialCond, 0, 2)

        return initialCond

    def _simulation_wrapper(self, initialCond):
        return myLib.create_simulation_loss(
            initialCond, self.P, self.T, self.w0, self.titf,
            self.C3Vals, self.s0, self.lambd, self.Gammas,
            self.absProj, self.deltas, self.alphas, self.spontCase, 
            self.nAtoms, self.nBeams, self.modulateAlpha, self.modFreq
        )

    def run_simulation_loss(self):
        initialConds = self._generateInitialCond()

        results = Parallel(n_jobs=self.n_jobs)(
            delayed(self._simulation_wrapper)(initialConds[i])
            for i in tqdm(range(self.n_samples))
        )
        return np.array(results)

    def survival_prob(self, spontCase='s'):
        results = self.run_simulation_loss()

        survProb = [np.sum(results[:, i]) / len(results[:, i]) for i in range(self.nAtoms)]
        phScatt = list(results[:, -1])  # Photon scattering

        return [survProb, phScatt]

    @staticmethod
    def convert_point_to_dash(value):
        return str(value).split('.')

    @staticmethod
    def cptd(value):
        return SimulationRunner.convert_point_to_dash(value)
