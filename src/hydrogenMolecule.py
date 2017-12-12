import numpy as np
from nbody import *

class HydrogenMolecule():

	def __init__(self, initialNucleiDistance):

		self.nucleiDistance = initialNucleiDistance
		self.nbody = None
		self.defineNbody()
		self.initializeElectronPositions()
		# self.wavefunction = self.trialWavefunction(0.1)

	def defineNbody(self):

		d = self.nucleiDistance
		numNuclei = 2
		numElectrons = 2
		nucleiCharges = np.array([1, 1])
		nucleiPositions = np.array([[-d/2.0, 0.0, 0.0], [d/2.0, 0, 0]]).transpose()
		electronPositions = np.zeros((3, 2))
		nbody = Nbody(numNuclei, numElectrons, nucleiCharges, nucleiPositions, electronPositions)

		self.nbody = nbody

	def initializeElectronPositions(self):

		d = self.nucleiDistance
		N = self.nbody.N
		electronPositions = np.random.rand(3, N)*2.0*d - d
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
