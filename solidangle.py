import matplotlib.pyplot as plt
import numpy as np
import math

plt.style.use('ggplot')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


temp = []

with open("dt=half_angles.dat", "r") as infile:	#reading file of all masses with respective fallback times
	for line in infile:
		if line.startswith("#"):
			continue
		else:
			temp.append(map(float, line.split()))





phi = [line[0] for line in temp]
costheta = [line[1] for line in temp]



#in C, save all particle data
#atan2(y,x)

#plt.plot(phi, costheta, 'ro')

#plot log(mass/unit surface) = log(mass) *4*pi/(number of bins)


#plt.ylabel(r'$\lamba$')


plt.xlabel(r'$\phi$ (angle from $x$-axis in radians)')
plt.ylabel(r'$cos(\theta)$ ($\theta$ is angle from $z$ axis)')
plt.title(r'$\Delta t = \frac12, 2$: Directions of velocities of the unbound particles')


from matplotlib.colors import LogNorm
from pylab import *

#normal distribution center at x=0 and y=5
#x = randn(100000)
#y = randn(100000)+5

hist2d(phi, costheta, bins=100, norm=LogNorm())
colorbar()

show()