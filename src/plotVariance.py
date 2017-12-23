import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from monteCarlo import *
from hydrogenMolecule import *

hm = HydrogenMolecule(1.4)
mc = MonteCarlo(hm, 2.0)

alpha = 1.1
powers = np.arange(8, 17)
numSample = np.power(2, powers)
energy = np.zeros(len(numSample))
var = np.zeros(len(numSample))


# prepare figure
fig, ax = plt.subplots(1,1)

for i in range(len(numSample)):
	
	N = numSample[i]
	en = np.zeros(N)
	for iteration in range(N):
		mc.sample(alpha)
		en[iteration] = mc.molecule.localEnergy(alpha)
	print 'N = ' + str(N) + ', accepted: ' + str(mc.numAccepted/(N + 0.0)*100) + '%'
	# calculate mean energy and variance 
	energy[i] = np.sum(en)/float(N)
	var[i] = np.sum((en - energy[i])**2)/float(N*(N - 1))

ax.plot(numSample, var)

ax.set_xlabel(r'N')
ax.set_ylabel(r'variance $\sigma^2$')
ax.set_xscale('log')

fig.savefig('variance.pdf')