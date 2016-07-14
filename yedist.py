import matplotlib.pyplot as plt
import numpy as np
import math

plt.style.use('ggplot')
#plt.rc('text', usetex = True)
plt.rc('font', family = 'serif')


temp = []

with open("solidangle.dat", "r") as infile: #reading file of all masses with respective fallback times
  for line in infile:
    if line.startswith("#"):
      continue
    else:
      temp.append(map(float, line.split()))


phi = [line[0] for line in temp]
costheta = [line[1] for line in temp]
ye = [line[3] for line in temp]
mass = [line[2] for line in temp]
ye_ = []

for i in range(0, len(ye)):
     ye_.append(ye[i]*mass[i])

#ye_ = [float(i) for i in ye_]

print len(ye_)
print len(ye)


mass_ = sum(mass)


#actual number of bins is the square of the below value
numofbins = 100

#scaling, so that we get mass/steradian,   for each square
ye = [elem*(numofbins**2)/(4.*math.pi*mass_) for elem in ye_]


plt.xlabel(r'$\phi$ (angle from $x$-axis in radians)')
plt.ylabel(r'$cos(\theta)$ ($\theta$ is angle from $z$ axis)')
plt.title(r'$\Delta t = 10$: Distribution of electron fractions for unbound particles')


from matplotlib.colors import LogNorm
from pylab import *

#normal distribution center at x=0 and y=5
#x = randn(100000)
#y = randn(100000)+5

hist2d(phi, costheta, bins=numofbins, norm=LogNorm(), weights = ye)

#show the log scale color legend
colorbar()

show()