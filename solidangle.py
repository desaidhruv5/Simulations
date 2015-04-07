import matplotlib.pyplot as plt
import numpy as np
import math
plt.style.use('ggplot')


temp = []

with open("dt=2_angles.dat", "r") as infile:	#reading file of all masses with respective fallback times
	for line in infile:
		if line.startswith("#"):
			continue
		else:
			temp.append(map(float, line.split()))





tanphi = [line[0] for line in temp]
costheta = [line[1] for line in temp]

phi = np.arctan(tanphi)

#plt.plot(phi, costheta, 'ro')

plt.xlabel('phi (angle from x-axis in radians)')
plt.ylabel('cos(theta) (theta is angle from z axis)')
plt.title('Directions of velocities of the unbound particles')


from matplotlib.colors import LogNorm
from pylab import *

#normal distribution center at x=0 and y=5
#x = randn(100000)
#y = randn(100000)+5

hist2d(phi, costheta, bins=40, norm=LogNorm())
colorbar()

show()