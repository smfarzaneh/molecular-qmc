# Python 2.7 class 
# S. M. Farzaneh, farzaneh@nyu.edu

import numpy as np
import scipy.constants as sc

class Potential(object):

	hbar = sc.hbar
	me = sc.electron_mass
	eps = sc.epsilon_0
	e = sc.elementary_charge
	pi = sc.pi
	k = 1.0/(4*pi*eps)
	a0 = hbar**2/(me*k*e**2)	# Bohr radius 
	E0 = k*e/a0 				# in units of [eV]

	def __init__(self, nucleusCharge, numElectrons, distancesToNucleus):
		
		self.Z = nucleusCharge
		self.n = numElectrons
		self.r = distancesToNucleus 	# in units of Bohr radius [a0]
		self.V = -1.0*self.E0

	def nucleusPotential(self):
		
		potentialEnergy = -self.E0*self.Z*np.sum(1.0/self.r)

		return potentialEnergy

	def setPotential(self):

		self.V = self.nucleusPotential() # TODO: + electron-electron interaction
	