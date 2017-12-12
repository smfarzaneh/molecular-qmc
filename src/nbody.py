from trialWavefunction import *
import numpy as np

class Nbody(object):

	def __init__(self, numNuclei, numElectrons, nucleiCharges, nucleiPositions, electronsPositions):
		
		self.K = numNuclei
		self.N = numElectrons
		self.Z = nucleiCharges
		self.setNucleiPositions(nucleiPositions)
		self.setElectronPositions(electronsPositions)

	def setNucleiPositions(self, nucleiPositions):
		
		self.Rx = nucleiPositions[0, :]
		self.Ry = nucleiPositions[1, :]
		self.Rz = nucleiPositions[2, :]

	def setElectronPositions(self, electronsPositions):
		
		self.rx = electronsPositions[0, :]
		self.ry = electronsPositions[1, :]
		self.rz = electronsPositions[2, :]

	def distanceElectronNucleus(self):

		# calculate all electron-nucleus distances 
		# return them as a matrix with N rows and K columns

		# tile positions into matrices to avoid nested loops
		Rx_col = np.tile(self.Rx, (self.N, 1))
		Ry_col = np.tile(self.Ry, (self.N, 1))
		Rz_col = np.tile(self.Rz, (self.N, 1))
		rx_row = np.tile(self.rx, (self.K, 1)).transpose()
		ry_row = np.tile(self.ry, (self.K, 1)).transpose()
		rz_row = np.tile(self.rz, (self.K, 1)).transpose()

		# calculate distances
		dx = np.subtract(rx_row, Rx_col)
		dy = np.subtract(ry_row, Ry_col)
		dz = np.subtract(rz_row, Rz_col)
		distance = np.sqrt(np.power(dx, 2) + np.power(dy, 2) + np.power(dz, 2))

		return distance

	def distanceElectronElectron(self):

		# calculate all electron-electron distances 
		# return them as a matrix with N rows and N columns

		# tile positions into matrices to avoid nested loops
		rx_col = np.tile(self.rx, (self.N, 1))
		ry_col = np.tile(self.ry, (self.N, 1))
		rz_col = np.tile(self.rz, (self.N, 1))
		rx_row = np.tile(self.rx, (self.N, 1)).transpose()
		ry_row = np.tile(self.ry, (self.N, 1)).transpose()
		rz_row = np.tile(self.rz, (self.N, 1)).transpose()

		# calculate distances
		dx = np.subtract(rx_row, rx_col)
		dy = np.subtract(ry_row, ry_col)
		dz = np.subtract(rz_row, rz_col)
		distance = np.sqrt(np.power(dx, 2) + np.power(dy, 2) + np.power(dz, 2))

		return distance

	def distanceNucleusNucleus(self):

		# calculate all nucleus-nucleus distances 
		# return them as a matrix with K rows and K columns

		# tile positions into matrices to avoid nested loops
		Rx_col = np.tile(self.Rx, (self.K, 1))
		Ry_col = np.tile(self.Ry, (self.K, 1))
		Rz_col = np.tile(self.Rz, (self.K, 1))
		Rx_row = np.tile(self.Rx, (self.K, 1)).transpose()
		Ry_row = np.tile(self.Ry, (self.K, 1)).transpose()
		Rz_row = np.tile(self.Rz, (self.K, 1)).transpose()

		# calculate distances
		dx = np.subtract(Rx_row, Rx_col)
		dy = np.subtract(Ry_row, Ry_col)
		dz = np.subtract(Rz_row, Rz_col)
		distance = np.sqrt(np.power(dx, 2) + np.power(dy, 2) + np.power(dz, 2))

		return distance

	def potentialElectronNucleus(self):

		D = self.distanceElectronNucleus()
		RD = np.reciprocal(D)
		Z = np.tile(self.Z, (self.N, 1))
		potential = -np.sum(np.multiply(Z, RD))

		return potential

	def potentialElectronElectron(self):

		D = self.distanceElectronElectron()
		I = np.identity(self.N)
		RD = np.reciprocal(D + I) - I
		potential = np.sum(RD)/2.0 		
		# divide by 2 is for counting electron-electron distances twice

		return potential

	def potentialNucleusNucleus(self):

		D = self.distanceNucleusNucleus()
		I = np.identity(self.K)
		RD = np.reciprocal(D + I) - I
		potential = np.sum(RD)/2.0 		
		# divide by 2 is for counting nucleus-nucleus distances twice

		return potential

	def potential(self):

		return self.potentialElectronNucleus() \
			+ self.potentialElectronElectron() \
			+ self.potentialNucleusNucleus()

	def kinetic(self):

		h = 0.001
		r = np.array([[self.rx], [self.ry], [self.rz]])
		r_forward = r + h
		r_backward = r - h
		wavefunction_r = TrialWavefunction.evaluate(r)
		wavefunction_r_forward = TrialWavefunction.evaluate(self.r_forward)
		wavefunction_r_backward = TrialWavefunction.evaluate(self.r_backward)
		secondDerivative = (wavefunction_r_forward + wavefunction_r_backward - 2.0*wavefunction_r)/h**2
		laplacian = np.sum(secondDerivative)
		kinetic = -0.5*laplacian/wavefunction_r

		return kinetic
















