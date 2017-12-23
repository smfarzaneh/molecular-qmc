from mpl_toolkits.mplot3d import Axes3D
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import numpy as np
from monteCarlo import *
from hydrogenMolecule import *

hm = HydrogenMolecule(1.4)
mc = MonteCarlo(hm, 2.0)
alpha = 1.1
N = 256
en = 0

# prepare figure
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim(-10, 10)
ax.set_ylim(-10, 10)
ax.set_zlim(-10, 10)

def plotElectrons(ax, mc):

	N = mc.molecule.nbody.N 
	rx = mc.molecule.nbody.rx[0]
	ry = mc.molecule.nbody.ry[0]
	rz = mc.molecule.nbody.rz[0]
	p1 = ax.scatter(rx, ry, rz, c='red', s = 10.0, cmap=plt.cm.autumn)
	rx = mc.molecule.nbody.rx[1]
	ry = mc.molecule.nbody.ry[1]
	rz = mc.molecule.nbody.rz[1]
	p2 = ax.scatter(rx, ry, rz, c='green', s = 10.0, cmap=plt.cm.autumn)

def plotNuclei(ax, mc):
	 
	 K = mc.molecule.nbody.K
	 for i in range(K):
	 	Rx = mc.molecule.nbody.Rx[i]
	 	Ry = mc.molecule.nbody.Ry[i]
	 	Rz = mc.molecule.nbody.Rz[i]
	 	ax.scatter(Rx, Ry, Rz, c=0.1, s = 40.0)

plotNuclei(ax, mc)

for iteration in range(N):
	plotElectrons(ax, mc)
	mc.sample(alpha)
	en += mc.molecule.localEnergy(alpha)

en = en/N
print 'alpha = ' + str(alpha) + ', energy = ' + str(en) + ', accepted= ' + str(mc.numAccepted/(N + 0.0)*100) + '%'


fig.savefig('scatter.pdf')