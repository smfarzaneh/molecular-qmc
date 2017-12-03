import unittest
from ..src.potential import *
import numpy as np

class TestPotential(unittest.TestCase):

	@classmethod
	def setUpClass(cls):
		pass

	def test_bohr_energy(self):
		hydrogen = Potential(nucleusCharge = 1, numElectrons = 1, distancesToNucleus = 1.0)
		self.assertAlmostEqual(hydrogen.nucleusPotential(), -27.2113860)
		hydrogen.setPotential()
		self.assertAlmostEqual(hydrogen.V, -27.2113860)

	def test_two_electrons(self):
		helium = Potential(nucleusCharge = 2, numElectrons = 2, distancesToNucleus = np.array([1.0, 1.0]))
		self.assertAlmostEqual(helium.nucleusPotential(), -108.8455441)

if __name__ == '__main__':
    unittest.main()