import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib.cm as cm
plt.style.use('ggplot')
plt.rc('font', family='serif')

initial = []
#USER INPUT############################
#tos is added to make evolution begin at correct time after merger (time offset)
tos = 346
#######################################


#1 solar mass = 1.99e33 grams
#1 km^3 = 1e9 cm^3

#1Msun/km^3 = 1.99e33g/1e9cm^3 = 1.99e24 g/cm^3

filename = []
for i in os.listdir(os.getcwd()):
  if i.startswith("densevo") and i.endswith(".dat"):
    filename.append(i)

colors = iter(cm.rainbow(np.linspace(0, 1, len(filename))))

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
for i in filename:

  with open(i, "r") as file:
    temp = []
    for line in file:
      if line.startswith("#"):
        continue
      else:
        temp.append(map(float, line.split()))

  time = [line[12] for line in temp]

  #scaling time
  time = [elem +tos for elem in time]
  time = [elem * 4.96*1e-6 for elem in time]

  log_time = np.log10(time)


  density = [line[8]/((1.48)**3) for line in temp]


  #convert to g/cm^3
  density = [1.99e24*x for x in density]


  x = [1.48*line[1] for line in temp]
  y = [1.48*line[2] for line in temp]
  z = [1.48*line[3] for line in temp]
  r = [np.sqrt(x[k]**2 + y[k]**2 + z[k]**2) for k in range(len(x))]


  log_density = np.log10(density)
  log_r = np.log10(r)
  #plt.plot(log_r, log_density)
  plt.plot(r, density)
  plt.yscale('log')
  plt.xscale('log')
#--------------------------------------------------------------
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


#plotting actual theoretical curve


#

#'''

p=[1e2,1e5]
q= [1e20, 1e11]
plt.plot(p,q, linewidth = 3, label = r' $r^{-3}$')


plt.xlabel(r'Distance from Black Hole [$km$]')
#'''


'''
plt.plot(log_time, log_density)
plt.xlabel('Time (\log(seconds)')
#'''

plt.ylabel(r'Density [$g$ $cm^{-3}$]')

plt.legend()

##USER EDIT############
plt.title(r'M5_S9: Density of Ejecta (r)')
#'''


#




plt.show()




#plot mass on y axis vs ye 






