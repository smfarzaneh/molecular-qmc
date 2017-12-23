import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from monteCarlo import *
from hydrogen import *
from helium import *
from hydrogenMoleculeIon import *
from hydrogenMolecule import *

def testA():

	# initialize Hydrogen atom and Monte Carlo
	initialStepSize = 2.5
	h = Hydrogen()
	mc = MonteCarlo(h, initialStepSize)
	
	# calculate energy and variance vs. variational params
	N = 2**18 					# number of samples
	m = 20 						
	alpha = np.linspace(0.1, 2.0, m)
	en = np.zeros(m)
	en2 = np.zeros(m)
	var = np.zeros(m)
	accepted = np.zeros(m)

	print 'Test A: Hydrogen Atom'
	print '# of Samples: ' + str(N)
	print '# of Variations: ' + str(m)

	for a in range(m):

		mc.stepSize = initialStepSize/alpha[a]
		for iteration in range(N):
			mc.sample(alpha[a])
			localEnergy = mc.molecule.localEnergy(alpha[a])
			en[a] += localEnergy
			en2[a] += localEnergy**2
		# calculate energy and variance
		en[a] /= float(N)
		var[a] = en2[a]/float(N) - en[a]**2
		print 'alpha('+ str(a + 1) +')=' + str(alpha[a]) + ', energy=' + str(en[a]) + ', variance=' + str(var[a]) + ', accepted: ' + str(mc.numAccepted*100/N) + '%' 
		mc.reset()

	# save data to file
	data = np.array([alpha.transpose(), en.transpose(), var.transpose()]).transpose()
	np.savetxt('test-a-alpha.out', data, delimiter='\t')
	print 'test-a-alpha.out was saved.'

def testBalpha():

	# initialize Helium atom and Monte Carlo
	initialStepSize = 2.5
	h = Helium()
	mc = MonteCarlo(h, initialStepSize)
	
	# calculate energy and variance vs. variational params
	N = 2**16 					# number of samples
	m = 20 						
	alpha = np.linspace(0.1, 2.0, m)
	beta = 1.0
	en = np.zeros(m)
	en2 = np.zeros(m)
	var = np.zeros(m)
	accepted = np.zeros(m)

	print 'Test B: Helium Atom, beta=' + str(beta)
	print '# of Samples: ' + str(N)
	print '# of Variations: ' + str(m)

	for a in range(m):

		varpar = np.array([alpha[a], beta])
		mc.stepSize = initialStepSize/alpha[a]
		for iteration in range(N):
			mc.sample(varpar)
			localEnergy = mc.molecule.localEnergy(varpar)
			en[a] += localEnergy
			en2[a] += localEnergy**2
		# calculate energy and variance
		en[a] /= float(N)
		var[a] = en2[a]/float(N) - en[a]**2
		print 'alpha('+ str(a + 1) +')=' + str(alpha[a]) + ', energy=' + str(en[a]) + ', variance=' + str(var[a]) + ', accepted: ' + str(mc.numAccepted*100/N) + '%' 
		mc.reset()

	# save data to file
	data = np.array([alpha.transpose(), en.transpose(), var.transpose()]).transpose()
	np.savetxt('test-b-alpha.out', data, delimiter='\t')
	print 'test-b-alpha.out was saved.'

def testBbeta():

	# initialize Helium atom and Monte Carlo
	initialStepSize = 2.5
	h = Helium()
	mc = MonteCarlo(h, initialStepSize)
	
	# calculate energy and variance vs. variational params
	N = 2**16 					# number of samples
	m = 20 						
	alpha = 1.8
	beta = np.linspace(0.1, 2.0, m)
	en = np.zeros(m)
	en2 = np.zeros(m)
	var = np.zeros(m)
	accepted = np.zeros(m)

	print 'Test B: Helium Atom, alpha=' + str(alpha)
	print '# of Samples: ' + str(N)
	print '# of Variations: ' + str(m)

	for a in range(m):

		varpar = np.array([alpha, beta[a]])
		mc.stepSize = initialStepSize/alpha
		for iteration in range(N):
			mc.sample(varpar)
			localEnergy = mc.molecule.localEnergy(varpar)
			en[a] += localEnergy
			en2[a] += localEnergy**2
		# calculate energy and variance
		en[a] /= float(N)
		var[a] = en2[a]/float(N) - en[a]**2
		print 'beta('+ str(a + 1) +')=' + str(beta[a]) + ', energy=' + str(en[a]) + ', variance=' + str(var[a]) + ', accepted: ' + str(mc.numAccepted*100/N) + '%' 
		mc.reset()

	# save data to file
	data = np.array([beta.transpose(), en.transpose(), var.transpose()]).transpose()
	np.savetxt('test-b-beta.out', data, delimiter='\t')
	print 'test-b-beta.out was saved.'

def testC():

	# initialize Hydrogen molecule ion and Monte Carlo
	initialStepSize = 2.5
	h = HydrogenMoleculeIon(2.0)
	mc = MonteCarlo(h, initialStepSize)
	
	# calculate energy and variance vs. variational params
	N = 2**16 					# number of samples
	m = 20 						
	alpha = np.linspace(0.1, 2.0, m)
	en = np.zeros(m)
	en2 = np.zeros(m)
	var = np.zeros(m)
	accepted = np.zeros(m)

	print 'Test C: Hydrogen Molecule Ion'
	print '# of Samples: ' + str(N)
	print '# of Variations: ' + str(m)

	for a in range(m):

		varpar = alpha[a]
		mc.stepSize = initialStepSize/alpha[a]
		for iteration in range(N):
			mc.sample(varpar)
			localEnergy = mc.molecule.localEnergy(varpar)
			en[a] += localEnergy
			en2[a] += localEnergy**2
		# calculate energy and variance
		en[a] /= float(N)
		var[a] = en2[a]/float(N) - en[a]**2
		print 'alpha('+ str(a + 1) +')=' + str(alpha[a]) + ', energy=' + str(en[a]) + ', variance=' + str(var[a]) + ', accepted: ' + str(mc.numAccepted*100/N) + '%' 
		mc.reset()

	# save data to file
	data = np.array([alpha.transpose(), en.transpose(), var.transpose()]).transpose()
	np.savetxt('test-c-alpha.out', data, delimiter='\t')
	print 'test-c-alpha.out was saved.'

 
def testDalpha():

	# initialize Hydrogen molecule ion and Monte Carlo
	initialStepSize = 2.5
	h = HydrogenMolecule(1.4)
	mc = MonteCarlo(h, initialStepSize)
	
	# calculate energy and variance vs. variational params
	N = 2**16 					# number of samples
	m = 20 						
	alpha = np.linspace(0.1, 2.0, m)
	beta = 1.0
	en = np.zeros(m)
	en2 = np.zeros(m)
	var = np.zeros(m)
	accepted = np.zeros(m)

	print 'Test D: Hydrogen Molecule, beta=' + str(beta)
	print '# of Samples: ' + str(N)
	print '# of Variations: ' + str(m)

	for a in range(m):

		varpar = np.array([alpha[a], beta])
		mc.stepSize = initialStepSize/alpha[a]
		for iteration in range(N):
			mc.sample(varpar)
			localEnergy = mc.molecule.localEnergy(varpar)
			en[a] += localEnergy
			en2[a] += localEnergy**2
		# calculate energy and variance
		en[a] /= float(N)
		var[a] = en2[a]/float(N) - en[a]**2
		print 'alpha('+ str(a + 1) +')=' + str(alpha[a]) + ', energy=' + str(en[a]) + ', variance=' + str(var[a]) + ', accepted: ' + str(mc.numAccepted*100/N) + '%' 
		mc.reset()

	# save data to file
	data = np.array([alpha.transpose(), en.transpose(), var.transpose()]).transpose()
	np.savetxt('test-d-alpha.out', data, delimiter='\t')
	print 'test-d-alpha.out was saved.'

def testDbeta():

	# initialize Hydrogen molecule ion and Monte Carlo
	initialStepSize = 2.5
	h = HydrogenMolecule(1.4)
	mc = MonteCarlo(h, initialStepSize)
	
	# calculate energy and variance vs. variational params
	N = 2**16 					# number of samples
	m = 20 						
	alpha = 1.3
	beta = np.linspace(0.1, 2.0, m)
	en = np.zeros(m)
	en2 = np.zeros(m)
	var = np.zeros(m)
	accepted = np.zeros(m)

	print 'Test D: Hydrogen Molecule, alpha=' + str(alpha)
	print '# of Samples: ' + str(N)
	print '# of Variations: ' + str(m)

	for a in range(m):

		varpar = np.array([alpha, beta[a]])
		mc.stepSize = initialStepSize/alpha
		for iteration in range(N):
			mc.sample(varpar)
			localEnergy = mc.molecule.localEnergy(varpar)
			en[a] += localEnergy
			en2[a] += localEnergy**2
		# calculate energy and variance
		en[a] /= float(N)
		var[a] = en2[a]/float(N) - en[a]**2
		print 'beta('+ str(a + 1) +')=' + str(beta[a]) + ', energy=' + str(en[a]) + ', variance=' + str(var[a]) + ', accepted: ' + str(mc.numAccepted*100/N) + '%' 
		mc.reset()

	# save data to file
	data = np.array([beta.transpose(), en.transpose(), var.transpose()]).transpose()
	np.savetxt('test-d-beta.out', data, delimiter='\t')
	print 'test-d-beta.out was saved.'

# testA()
# testBalpha()
# testBbeta()
# testC()
# testDalpha()
testDbeta()