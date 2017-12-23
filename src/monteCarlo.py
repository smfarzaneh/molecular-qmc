import numpy as np

class MonteCarlo(object):

	def __init__(self, molecule, stepSize):

		self.molecule = molecule
		self.stepSize = stepSize
		self.numAccepted = 0

	def sample(self, alpha):

		N = self.molecule.nbody.N
		whichElectron = np.random.randint(N)
		deltaPosition = self.stepSize*(np.random.rand(3) - 0.5)
		rx = np.copy(self.molecule.nbody.rx).transpose()
		ry = np.copy(self.molecule.nbody.ry).transpose()
		rz = np.copy(self.molecule.nbody.rz).transpose()
		r = np.array([rx, ry, rz])
		dr = np.zeros((3, N))
		dr[:, whichElectron] = deltaPosition
		
		# Calculate probabilities
		currentWavefunction = self.molecule.trialWavefunction(alpha)
		currentProbability = currentWavefunction**2
		self.molecule.nbody.setElectronPositions(r + dr)
		nextWavefunction = self.molecule.trialWavefunction(alpha)
		nextProbability = nextWavefunction**2
		
		# Metropolis Test
		s = np.random.rand()
		ratio = nextProbability/currentProbability
		if ratio >= 1:
			accepted = True
			self.numAccepted += 1
		elif ratio > s:
			self.numAccepted += 1
			accepted = True
		else:
			self.molecule.nbody.setElectronPositions(r)
			accepted = False 