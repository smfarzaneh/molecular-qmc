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
N = 1024
en = 0

# prepare figure
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim(-10, 10)
ax.set_ylim(-10, 10)
ax.set_zlim(-10, 10)

def samplingRate(numMonteCarloSamples):

	fps = 24
	duration = 20
	numFrames = fps*duration
	k = numMonteCarloSamples/duration/fps

	return k, numFrames

def prepareCamera():

	FFMpegWriter = ani.writers['ffmpeg']
	metadata = dict( title  = 'QMC', artist='Matplotlib', 
		comment= 'Quantum Monte Carlo for Hydrogen Molecule')
	writer = FFMpegWriter(fps=24, metadata=metadata)

	return writer

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

k, numFrames = samplingRate(N)
writer = prepareCamera()

with writer.saving(fig, "test.mp4", numFrames):

    for frame in range(0, numFrames):
        print "frame " + str(frame)

        ax.cla()
        ax.set_xlim(-10, 10)
        ax.set_ylim(-10, 10)
        ax.set_zlim(-10, 10)
        plotNuclei(ax, mc)
       	plotElectrons(ax, mc)
        titleString = 'N =' + str(frame*k)
        plt.title(titleString)
        writer.grab_frame()

        for i in np.arange(k):
			mc.sample(alpha)
			en += mc.molecule.localEnergy(alpha)
        
        
# fig.savefig('d.pdf')
print en/N
print 'accepted: ' + str(mc.numAccepted/(N + 0.0)*100) + '%'
