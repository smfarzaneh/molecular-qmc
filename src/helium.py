import numpy as np
from nbody import *

class Helium():

	def __init__(self):

		self.nbody = None
		self.defineNbody()
		self.initializeElectronPositions()

	def defineNbody(self):

		numNuclei = 1
		numElectrons = 2
		nucleiCharges = 2
		nucleiPositions = np.array([[0.0, 0.0, 0.0]]).transpose()
		electronPositions = np.zeros((3, numElectrons))
		nbody = Nbody(numNuclei, numElectrons, nucleiCharges, nucleiPositions, electronPositions)

		self.nbody = nbody

	def initializeElectronPositions(self):

		N = self.nbody.N
		electronPositions = np.random.rand(3, N)*5.0 - 2.5
		self.nbody.setElectronPositions(electronPositions)

	def trialWavefunction(self, alpha):

		r = self.nbody.distanceElectronNucleus()
		ree = self.nbody.distanceElectronElectron() + np.identity(self.nbody.N)
		slater = np.sum(np.exp(-1.0*(alpha[0])*r), axis = 1)
		jastrow = np.sqrt(np.prod(np.exp(ree/(2.0*(1.0 + alpha[1]*ree)))))
		singleElectronWavefunction = slater*jastrow
		totalWavefunction = np.prod(singleElectronWavefunction)

		return totalWavefunction

	def kinetic(self, alpha):

		h = 0.001
		rx = np.copy(self.nbody.rx).transpose()
		ry = np.copy(self.nbody.ry).transpose()
		rz = np.copy(self.nbody.rz).transpose()
		r = np.array([rx, ry, rz])
		self.nbody.setElectronPositions(r + h)
		wavefunction_r_forward = self.trialWavefunction(alpha)
		self.nbody.setElectronPositions(r - h)
		wavefunction_r_backward = self.trialWavefunction(alpha)
		self.nbody.setElectronPositions(r)
		wavefunction_r = self.trialWavefunction(alpha)
		secondDerivative = (wavefunction_r_forward + wavefunction_r_backward - 2.0*wavefunction_r)/h**2
		laplacian = np.sum(secondDerivative)
		kineticEnergy = -0.5*laplacian/wavefunction_r

		return kineticEnergy

	def localEnergy(self, alpha):

	 	return self.kinetic(alpha) + self.nbody.potential()
