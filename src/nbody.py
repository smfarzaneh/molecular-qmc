import numpy as np

class Nbody(object):

	def __init__(self, numNuclei, numElectrons, nucleiCharges, nucleiPositions, electronsPositions):
		
		self.K = numNuclei
		self.N = numElectrons
		self.Z = nucleiCharges
		self.Rx = nucleiPositions[0, :]
		self.Ry = nucleiPositions[1, :]
		self.Rz = nucleiPositions[2, :]
		self.rx = electronsPositions[0, :]
		self.ry = electronsPositions[1, :]
		self.rz = electronsPositions[2, :]

	def distanceToNuclei(self):

		# calculate all electrons-nuclei distances 
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

	def distanceToElectrons(self):

		# calculate all electrons-electrons distances 
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

	def nucleiPotential(self):

		D = self.distanceToNuclei()
		RD = np.reciprocal(D)
		Z = np.tile(self.Z, (self.N, 1))
		potential = -np.sum(np.multiply(Z, RD))

		return potential

	def electronsPotential(self):

		D = self.distanceToElectrons()
		I = np.identity(self.N)
		RD = np.reciprocal(D + I) - I
		potential = np.sum(RD)/2.0 		
		# divide by 2 is for counting electron-electron distances twice

		return potential

	def potential(self):

		return self.nucleiPotential() + self.electronsPotential()

	