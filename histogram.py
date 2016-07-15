import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pylab import *
from matplotlib.colors import LogNorm
import os
import matplotlib.cm as cm
from numpy import linalg as la
from numpy import matrix as mat



plt.style.use('ggplot')
plt.rc('font', family='serif')
pi=np.pi

MASS = 6.0552018


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
  
  #position components
  x = [line[1] for line in FILE]
  y = [line[2] for line in FILE]
  z = [line[3] for line in FILE]

  #radial vector
  radius = [np.sqrt(x[k]**2 + y[k]**2 + z[k]**2) for k in range(len(x))]

  #velocity components
  vx = [line[4] for line in FILE]
  vy = [line[5] for line in FILE]
  vz = [line[6] for line in FILE]

  v = [[vx[i], vy[i], vz[i]] for i in range(len(vx))  ]
  #radial velocity
  speed = la.norm(v, axis =1)



  #energy
  energy   = [line[7] for line in FILE]      #actually specific energy. units of c^2
  
  #angle from x-axis
  polar_angle = [line[0] for line in FILE]    #radians

  mass = [line[0]*1.99e33 for line in FILE]   #grams



  return radius, speed, energy, polar_angle , mass

#-----------------------------------------------

def remove_bound_particles(radius_, radial_velocity_, energy_, polar_angle_ , mass_):
  radial_velocity  = []
  polar_angle = []
  energy   = []
  mass=[]
  uncounted_mass=[]
  false_energy=[]
  radius = []
  for j in range(len(vr_)):
    if  energy_[j]>1:
      mass.append(mass_[j])
      radial_velocity.append(radial_velocity_[j])
      polar_angle.append(polar_angle_[j])
      energy.append(energy_[j])
      radius.append(radius_[j])

      if radial_velocity_[j] < 0:
        false_energy.append(mass_[j])
    else:
      uncounted_mass.append(mass_[j])
    

    """
    print "fraction of mass left out: ", sum(uncounted_mass)/sum(mass)
    print "fraction of mass with u>1, vr<0 (falsely marked): ", sum(false_energy)/sum(mass)
    """
  return radius, radial_velocity, energy, polar_angle, mass

#---------------------------------------------










#plt.xscale('log')
#plt.yscale('log')
#''' This is mass rate plot:

##USER EDIT############
#plt.title(r'M5_S9: Ejecta Velocity Distribution')
#'''


"""
B= 10**(B)

B= 100


SLICED=np.linspace(min(phi), max(phi),num=13)#, endpoint=False )
count=0

lower = 3
upper = 5

allm = []
allv=[]
for j in range(len(SLICED)-1):
  editmasses=[]
  editvr=[]
  for i in range(len(vr)):
    if  SLICED[j] < phi[i] <=SLICED[j+1]:
      editmasses.append(masses[i])
      editvr.append(vr[i])
  allm.append(editmasses)
  allv.append(editvr)
  if SLICED[lower+6]<=SLICED[j]<=SLICED[upper+6]:
    plt.hist(editvr, weights=editmasses, bins = B, alpha=.5, label= r"$ %f \leq$ $\phi$ $< %f$" %(SLICED[j], SLICED[j+1]) )
    count =count+1
print count, " slices plotted."
"""
##############################################################################################################
###############################     PLOTTING FUNCTIONS            #############################################
##############################################################################################################

def predict_velocities_from_energies(u, r):
  predicted_velocities = []
  for i in range(len(u)):
    predicted_velocities.append(   np.sqrt( (u[i]-1+MASS/r[i])*2  )   )
    #predicted_velocities[i] = 3e5*predicted_velocities[i]
  return predicted_velocities




#bins masses according respective radial velocities
def plot_velocity_spread(velocities, masses, num_of_bins, label):
  hist, bin_edges, ignored = plt.hist(velocities, weights=masses, bins = num_of_bins, alpha=.5, label=label)
  plt.title(r'%s: Spread of ejecta velocities' %NAME)
  plt.xlabel(r'Velocities $[km/s]$')
  plt.ylabel(r'Mass $[g]$')
  plt.savefig('%s_velocity_spread.png' %NAME)
  return hist, bin_edges

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
B=50      #number of bins
NAME= "M5_S9"
########EDIT LINE IMMEDIATELY BELOW THIS

FILE = extract("solidangle.dat")

r_, vr_, u_, phi_ , m_ = set_variables(FILE)


r, vr, u, phi, m = remove_bound_particles(r_, vr_, u_, phi_ , m_)




#"""

"""plt.xscale('log')
plt.yscale('log')"""

pv = predict_velocities_from_energies(u ,r)

B= np.linspace(min(vr+pv), max(vr+pv), 101)
counts1, bins = plot_velocity_spread(vr, m, B, "actual")
counts2, bins = plot_velocity_spread(pv, m, B, "predicted by initial energies")
#plt.close()



bins = [((bins[i+1]-bins[i])/2).tolist() for i in range(len(bins.tolist())-1)]

#plt.plot(bins, (np.array(counts1)-np.array(counts2)).tolist() )#, bins = len(bins)+1)


#"""

#plot_avg_velocity_vs_polar_angle(vr, m, phi, B)


"""
#plot the distribution of energies for those with velocities below the cutoff

plt.hist(velocities, weights=masses, bins = num_of_bins, alpha=.5, label=label)
"""
#want to know how much mass is at each energy, for those above v_cut, and below 

#then we can also plot mass-avged u 




#plot_2dhist(phi, vr, m, num_of_bins = 150)


plt.legend()
plt.show()


# remove saving feature from code, see if that speeds up evolution
#plot Energy(t) for unbound particles
#if bad, decrease time step. do adaptive time step
#for adaptive, make step proportional to distance from BH. Have all particles update with large dt,
#but closer ones have smaller updates

#plot mass on y axis vs ye 






