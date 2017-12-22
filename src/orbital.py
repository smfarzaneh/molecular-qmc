import numpy as np
from scipy.special import sph_harm as sh

class Orbital(object):

	def __init__(self):

		pass

	@classmethod
	def sphericalHarmonic(self, l, ml, theta phi):

		return sh(ml, l, theta, phi)