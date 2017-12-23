import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# load data 
dataAalpha = np.genfromtxt('test-a-alpha.out', delimiter='\t')
dataBalpha = np.genfromtxt('test-b-alpha.out', delimiter='\t')
dataBbeta = np.genfromtxt('test-b-beta.out', delimiter='\t')
dataCalpha = np.genfromtxt('test-c-alpha.out', delimiter='\t')
dataDalpha = np.genfromtxt('test-d-alpha.out', delimiter='\t')
dataDbeta = np.genfromtxt('test-d-beta.out', delimiter='\t')

# test A
fig, ax = plt.subplots(1,1)
ax.plot(dataAalpha[:, 0], dataAalpha[:, 1], 'k-o')
ax2 = ax.twinx()
ax2.plot(dataAalpha[:, 0], dataAalpha[:, 2], 'r:.')
ax.set_title(r'Test A: Hydrogen Atom, $N = 2^{18}$ = 262,144')
ax.set_xlabel(r'$\alpha$')
ax.set_ylabel(r'Local Energy [a.u.]')
ax.tick_params('y', colors='k')
ax2.tick_params('y', colors='r')
ax2.set_ylabel(r'Variance', color='r')
fig.savefig('test-a.pdf')
print 'test-a.pdf was saved'


# test B
fig, ax = plt.subplots(1,1)
ax.plot(dataBalpha[:, 0], dataBalpha[:, 1], 'k-o')
ax2 = ax.twinx()
ax2.plot(dataBalpha[:, 0], dataBalpha[:, 2], 'r:.')
ax.set_title(r'Test B: Helium Atom, $N = 2^{16}$ = 65,536, $\beta=1.0$')
ax.set_xlabel(r'$\alpha$')
ax.set_ylabel(r'Local Energy [a.u.]')
ax.tick_params('y', colors='k')
ax2.tick_params('y', colors='r')
ax2.set_ylabel(r'Variance', color='r')
fig.savefig('test-b-alpha.pdf')
print 'test-b-alpha.pdf was saved'

fig, ax = plt.subplots(1,1)
ax.plot(dataBbeta[:, 0], dataBbeta[:, 1], 'k-o')
ax2 = ax.twinx()
ax2.plot(dataBbeta[:, 0], dataBbeta[:, 2], 'r:.')
ax.set_title(r'Test B: Helium Atom, $N = 2^{16}$ = 65,536, $\alpha=1.8$')
ax.set_xlabel(r'$\beta$')
ax.set_ylabel(r'Local Energy [a.u.]')
ax.tick_params('y', colors='k')
ax2.tick_params('y', colors='r')
ax2.set_ylabel(r'Variance', color='r')
fig.savefig('test-b-beta.pdf')
print 'test-b-beta.pdf was saved'


# test C
fig, ax = plt.subplots(1,1)
ax.plot(dataCalpha[:, 0], dataCalpha[:, 1], 'k-o')
ax2 = ax.twinx()
ax2.plot(dataCalpha[:, 0], dataCalpha[:, 2], 'r:.')
ax.set_title(r'Test C: Hydrogen Molecule Ion, $N = 2^{16}$ = 65,536')
ax.set_xlabel(r'$\alpha$')
ax.set_ylabel(r'Local Energy [a.u.]')
ax.tick_params('y', colors='k')
ax2.tick_params('y', colors='r')
ax2.set_ylabel(r'Variance', color='r')
fig.savefig('test-c.pdf')
print 'test-c.pdf was saved'

# test D
fig, ax = plt.subplots(1,1)
ax.plot(dataDalpha[:, 0], dataDalpha[:, 1], 'k-o')
ax2 = ax.twinx()
ax2.plot(dataDalpha[:, 0], dataDalpha[:, 2], 'r:.')
ax.set_title(r'Test B: Hydrogen Molecule, $N = 2^{16}$ = 65,536, $\beta=1.0$')
ax.set_xlabel(r'$\alpha$')
ax.set_ylabel(r'Local Energy [a.u.]')
ax.tick_params('y', colors='k')
ax2.tick_params('y', colors='r')
ax2.set_ylabel(r'Variance', color='r')
fig.savefig('test-d-alpha.pdf')
print 'test-d-alpha.pdf was saved'

fig, ax = plt.subplots(1,1)
ax.plot(dataDbeta[:, 0], dataDbeta[:, 1], 'k-o')
ax2 = ax.twinx()
ax2.plot(dataDbeta[:, 0], dataDbeta[:, 2], 'r:.')
ax.set_title(r'Test B: Hydrogen Molecule, $N = 2^{16}$ = 65,536, $\alpha=1.3$')
ax.set_xlabel(r'$\beta$')
ax.set_ylabel(r'Local Energy [a.u.]')
ax.tick_params('y', colors='k')
ax2.tick_params('y', colors='r')
ax2.set_ylabel(r'Variance', color='r')
fig.savefig('test-d-beta.pdf')
print 'test-d-beta.pdf was saved'
