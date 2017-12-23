import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from monteCarlo import *
from hydrogen import *

# initialization
hm = Hydrogen()
mc = MonteCarlo(hm, 2.5)

alpha = 1.0
powers = np.arange(8, 18)
numSample = np.power(2, powers)
energy = np.zeros(len(numSample))
energy2 = np.zeros(len(numSample))
var = np.zeros(len(numSample))
N = 0
lastEnergy = 0
lastEnergy2 = 0

# prepare figure
fig, ax = plt.subplots(1,1)

for i in range(len(numSample)):
	
	N = numSample[i] - N
	for iteration in range(N):
		mc.sample(alpha)
		localEnergy = mc.molecule.localEnergy(alpha)
		energy[i] += localEnergy
		energy2[i] += localEnergy**2
	# calculate mean energy and variance 
	energy[i] += lastEnergy
	energy2[i] += lastEnergy2
	energy[i] /= float(numSample[i])
	var[i] = energy2[i]/float(numSample[i]) - energy[i]**2
	print 'N = ' + str(numSample[i]) + ', energy=' + str(energy[i]) + ', variance=' + str(var[i]) + ', accepted: ' + str(mc.numAccepted*100/N) + '%'
	lastEnergy = energy[i]*numSample[i]
	lastEnergy2 = energy2[i]
	

ax.plot(numSample, var, 'k-*')

ax.set_xlabel(r'N')
ax.set_ylabel(r'variance $\sigma^2$')
ax.set_xscale('log')

fig.savefig('variance.pdf')