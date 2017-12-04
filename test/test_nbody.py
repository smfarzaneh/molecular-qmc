import unittest
from ..src.nbody import *
import numpy as np

class TestNbody(unittest.TestCase):

	@classmethod
	def setUpClass(cls):
		numNuclei = 2
		numElectrons = 2
		nucleiCharges = np.array([1, 2])
		nucleiPositions = np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]]).transpose()
		electronPositions = np.array([[4.0, 0.0, 0.0], [8.0, 0.0, 0.0]]).transpose()
		cls._molecule = Nbody(numNuclei, numElectrons, nucleiCharges, nucleiPositions, electronPositions)

	def test_init(self):
		self.assertItemsEqual(self._molecule.Rx, np.array([0.0, 10.0]))
		self.assertItemsEqual(self._molecule.rx, np.array([4.0, 8.0]))

	def test_potential(self):
		self.assertAlmostEqual(self._molecule.nucleiPotential(), -41.0/24.0)
		self.assertAlmostEqual(self._molecule.electronsPotential(), 1.0/4.0)

if __name__ == '__main__':
    unittest.main()