
import matplotlib.pyplot as plt
import numpy as np
'''dt = 	[0, 0.0001, 		0.0005,			0.001, 		0.002,		0.005, 		0.01,		0.02]
error = [0, 0.00117285, 0.000679399, 0.000780327, 0.000813678, 0.00182792, 0.00681031, 0.0168093]
plt.plot(dt, error)
plt.xlabel('dt')
plt.ylabel('error')
plt.title('Error vs. time step')
plt.show()'''



time = []

with open("halfdt corrected velocities.dat", "r") as infile:
	for line in infile:
		if line.startswith("#"):
			continue
		else:
			time.append(map(float, line.split()))

temp = [line[8] for line in time]



#scaling time
# temp = [elem * 4.96*1e-6 for elem in temp]



log_temp = np.log(temp)


mass = [line[0] for line in time]

#defining cummulative mass
cu_mass = []
cu_mass.append(mass[0])

for i in range(1, len(mass)):
	cu_mass.append(cu_mass[i-1] + mass[i] )


log_cu_mass = np.log(cu_mass)

#defining mass fall back RATE
from scipy.interpolate import interp1d

f2 = interp1d(temp, cu_mass, kind='cubic')

'''massrate = []
for i in range(1, len(mass)):

	print i, temp[i], temp[i-1]
	massrate.append( (cu_mass[i] - cu_mass[i-1]) / (temp[i] - temp[i-1])  )
'''


#1507518

particle = (range(1, 1507519) )

plt.plot(temp, cu_mass)
#plt.plot(temp, massrate)
plt.xlabel('time')
plt.ylabel('mass')
plt.title('Mass as a function of time')
plt.show()











