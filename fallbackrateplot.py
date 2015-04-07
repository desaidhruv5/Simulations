import matplotlib.pyplot as plt
import numpy as np
plt.style.use('ggplot')

temp = []

with open("dt=2_corrected.dat", "r") as infile:	#reading file of all masses with respective fallback times
	for line in infile:
		if line.startswith("#"):
			continue
		else:
			temp.append(map(float, line.split()))

time = [line[11] for line in temp]

#scaling time
time = [elem +791 for elem in time]
time = [elem * 4.96*1e-6 for elem in time]

log_time = np.log10(time)



mass = [line[0] for line in temp]

#for each time interval, add up total mass in mass vector



#defining cummulative mass
cu_mass = []
cu_mass.append(mass[0])

for i in range(1, len(mass)):
	cu_mass.append(cu_mass[i-1] + mass[i] )



log_cu_mass = np.log10(cu_mass)

log_mass = np.log10(mass)

#defining mass fall back RATE
from scipy.interpolate import interp1d

#f2 = interp1d(time, cu_mass, kind='cubic')

'''massrate = []
for i in range(1, len(mass)):

	print i, time[i], time[i-1]
	massrate.append( (cu_mass[i] - cu_mass[i-1]) / (time[i] - time[i-1])  )
'''


N=80
log_tbin = (log_time[-1] - log_time[0])  /N

print "total time: " , (log_time[-1] - log_time[0])
print "t_bin: " , log_tbin

i=0
binmass = []
for j in range(1, N):
	#print j
	if  i == len(log_time):
		break
	binmass.append(0)
	while log_time[i] <= log_time[0] + j* log_tbin:			#while time is less than bin size + intial time,
		binmass[j-1] = binmass[j-1] + mass[i]
		i= i+1
		if i == len(log_time):
			break

tbin = []
for j in range(1,N):
	tbin.append(      time[0]*     (  10**((j)*log_tbin) - 10**((j-1)*log_tbin)  )           )


mrate = []
for j in range(N-1):
	mrate.append(binmass[j]/tbin[j])

log_mrate = np.log10(mrate)




log_bintime = []
for j in range(1,N):
	print j
	log_bintime.append(log_time[0]+j*log_tbin)


#print "final time is actually ", log_time[-1], "vs. ", log_time[0]+N*t_bin


#1507518

#particle = (range(1, 1507519) )

#plt.plot(time, cu_mass)



#plotting time vs massrate



#''' This is mass rate plot:
print "bin_time length: ", len(log_bintime)
plt.plot(log_bintime, log_mrate)
plt.ylabel('massrate (log Solar Masses)')
plt.title('Log-Log: Massrate as a function of time')
#'''


'''This is cummulative mass plot:

plt.xlabel('time (seconds)')
plt.ylabel('cummulative mass (Solar Masses)')
plt.title('Mass as a function of time')
plt.plot(log_time, log_cu_mass)
'''

plt.xlabel('time (log seconds)')
plt.show()




#plot mass on y axis vs ye 






