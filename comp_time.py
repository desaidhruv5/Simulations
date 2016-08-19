import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pylab import *
from matplotlib.colors import LogNorm

plt.style.use('ggplot')
plt.rc('font', family='serif')
pi=np.pi




#EXTRACT DATA
def extract(filename):
  matrix = []
  with open(filename, "r") as infile: #reading file of all masses with respective fallback times
    for line in infile:
      if line.startswith("#"):
        continue
      else:
        matrix.append(map(float, line.split()))
  return matrix

#----------------------------------------------------------1
def set_variables(file):
  
  time = [line[1] for line in FILE]
  percent = [line[0] for line in FILE]

  
  return time, percent

#-----------------------------------------------









#''' This is mass rate plot:

##USER EDIT############
##############################################################################################################
###############################     PLOTTING FUNCTIONS            #############################################
##############################################################################################################




#bins masses according respective radial velocities
def plot_velocity_spread(velocities, masses, num_of_bins, label):
  plt.hist(velocities, weights=masses, bins = num_of_bins, alpha=.5, label=label)
  plt.title(r'%s: Spread of ejecta velocities' %NAME)
  plt.xlabel(r'Velocities $[km/s]$')
  plt.ylabel(r'Mass $[g]$')
  plt.savefig('%s_velocity_spread.png' %NAME)

def plot_2dhist(phi, vr, m, num_of_bins):
  hist2d(phi, vr, bins= num_of_bins, weights=m, norm=LogNorm())
  plt.colorbar().set_label(label=r'Mass [$g$]')
  plt.xlabel(r'$\phi$ [angle from $x$-axis]')
  plt.ylabel(r'Radial velocity $v_{r}$ $[km/s]$')
  plt.title(r'%s: Distribution of mass over velocity & polar angle' %NAME)
  plt.savefig('%s_mass_distribution_vr_phi.png' %NAME)


def plot_avg_velocity_vs_polar_angle(velocity, mass, polar_angle, num_of_bins, cut = None):


#average v (mass weighted) for a given phi
#binmass[j] = total mass in jth bin

#Find total binmass

  weighted_velocity  = []   #weighted radial velocities
  for j in range(len(mass)):
    weighted_velocity.append(mass[j]*velocity[j])

  

  B = np.linspace(min(phi), max(phi), num_of_bins+1)
  binmass=np.zeros(len(B))
  for i in range(len(polar_angle)):
    j=0
    while not (B[j] <= polar_angle[i] <= B[j+1]): #sort ith particle in proper phi bin
      j=j+1
    binmass[j]=binmass[j]+mass[i]

  #Divide each particle by total bin mass
  avg_v=[]
  for i in range(len(mass)):
    j=0
    while not (B[j] <= phi[i] <= B[j+1]):
      j=j+1
    avg_v.append(weighted_velocity[i]/binmass[j])

  plt.hist(phi, weights=avg_v, bins = B, alpha=.5)
  plt.xlabel(r'$\phi$ [angle from x-axis]')
  plt.ylabel(r'Mass-weighted average $v_{r}$ [km/s]')
  plt.title(r'%s: Average ejecta velocities over polar angle' %NAME)
  plt.savefig('%s_avg_ejecta_velocity_over_phi.png' %NAME)

##############################################################################################################
##############################################################################################################
####################################                                    #####################################
####################################          USER INPUTS BELOW         ######################################
####################################                                    ######################################
##############################################################################################################
##############################################################################################################



#USER INPUT#############
NAME= "M5_S9"
########






FILE = extract("comp_time.dat")

time, computation = set_variables(FILE)


plt.plot(time, computation, label = "actual computation time")


#look back 10 steps to approximate slope
index = -30
slope= (computation[-1]-computation[index])/(time[-1]-time[index])
intercept = computation[-1]-slope*time[-1]

zero = -intercept/slope
final_time = (100-intercept)/slope
X= [zero,final_time]
Y= [0,100]
hours = int(floor(time[-1]/3600))
#plt.plot(X,Y, label = "Predicted time to completion")

print "Time taken so far: ",hours, ":", int((time[-1]%3600)/60), ":", int(time[-1]%60)
print "Predicted total time for simulation: ", final_time/3600, " hours."


plt.xlabel(' time [seconds]')
plt.ylabel('Percentage of Evolution Complete [%]')
plt.title(r'Computation time')
  
#plt.xscale('log')
#plt.yscale('log')

plt.legend(loc = 'upper left')
plt.show()

#figure out amount of time it takes to loop through N particles using method 1 and method 2.


#plot mass on y axis vs ye 






