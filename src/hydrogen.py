import numpy as np
import nbody

class Hydrogen():

	def __init__(self):

		# define Hydrogen
		numNuclei = 1
		numElectrons = 1
		nucleiCharges = 1
		nucleiPositions = np.array([0.0, 0.0, 0.0]).transpose()
		electronPositions = np.array([10.0, 0.0, 0.0]).transpose()
		Nbody.__init__(numNuclei, numElectrons, nucleiCharges, nucleiPositions, electronPositions)
