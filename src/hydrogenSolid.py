import numpy as np
from nbody import *

class HydrogenSolid():

	def __init__(self, initialNucleiDistance, numNuclei):

		self.nucleiDistance = initialNucleiDistance
		self.numNuclei = numNuclei
		self.nbody = None
		self.defineNbody()
		self.initializeElectronPositions()

	def defineNbody(self):

		d = self.nucleiDistance
		K = self.numNuclei
		N = K
		nucleiCharges = np.ones(K)
		x_vals = np.linspace(-(np.sqrt(K) - 1)*d/2.0, (np.sqrt(K) - 1)*d/2.0, np.sqrt(K))
		y_vals = x_vals
		z_vals = x_vals*0
		nucleiPositions = np.zeros((3, K))
		nucleiPositions[0, :] = x_vals
		nucleiPositions[1, :] = y_vals
		nucleiPositions[2, :] = z_vals
		electronPositions = np.zeros((3, N))
		nbody = Nbody(numNuclei, numElectrons, nucleiCharges, nucleiPositions, electronPositions)

		self.nbody = nbody

	def initializeElectronPositions(self):

		d = self.nucleiDistance
		N = self.nbody.N
		electronPositions = np.random.rand(3, N)*20.0*d - 10*d
		self.nbody.setElectronPositions(electronPositions)

	def trialWavefunction(self, alpha):
		
		r = self.nbody.distanceElectronNucleus()
		singleElectronWavefunction = np.sum(np.exp(-alpha*r), axis = 1)
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
