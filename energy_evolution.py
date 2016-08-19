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

#MASS of BH = 


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
  #in km
  x = [line[1] for line in FILE]
  y = [line[2] for line in FILE]
  z = [line[3] for line in FILE]
  r = [ [x[k], y[k], z[k]] for k in range(len(x))] 

  #radial vector
  #radius = [np.sqrt(x[k]**2 + y[k]**2 + z[k]**2) for k in range(len(x))]

  #velocity components
  vx = [line[4] for line in FILE]
  vy = [line[5] for line in FILE]
  vz = [line[6] for line in FILE]
  v = [[vx[i], vy[i], vz[i]] for i in range(len(vx))  ]

  t = [line[12]/2.27e5 for line in FILE]
  
  #energy
  energy   = [line[7] for line in FILE]      #actually specific energy. units of c^2
  
  #angle from x-axis

  mass = [line[0]*1.99e33 for line in FILE]   #grams



  return r, v, energy, mass, t

#-----------------------------------------------

def remove_bound_particles(radius_, velocity_, energy_, mass_, t_):
 
  velocity  = []
  energy   = []
  mass=[]
  uncounted_mass=[]
  false_energy=[]
  radius = []
  t = []

  for j in range(len(mass_)):
    if  energy_[j]>1:
      mass.append(mass_[j])
      velocity.append(velocity_[j])
      energy.append(energy_[j])
      radius.append(radius_[j])
      t.append(t_[j])

      if velocity_[j] < 0:
        false_energy.append(mass_[j])
    else:
      uncounted_mass.append(mass_[j])
    

    """
    print "fraction of mass left out: ", sum(uncounted_mass)/sum(mass)
    print "fraction of mass with u>1, vr<0 (falsely marked): ", sum(false_energy)/sum(mass)
    """
  return radius, velocity, energy, mass, t

#---------------------------------------------









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



##############################################################################################################
##############################################################################################################
####################################                                    #####################################
####################################          USER INPUTS BELOW         ######################################
####################################                                    ######################################
##############################################################################################################
##############################################################################################################



#USER INPUT#############
B=50      #number of bins
NAME= "M7_S9_I20"
########EDIT LINE IMMEDIATELY BELOW THIS

filename = []
for i in os.listdir(os.getcwd()):
  if i.startswith("densevo") and i.endswith(".dat"):
    filename.append(i)

colors = iter(cm.rainbow(np.linspace(0, 1, len(filename))))
MASS= 8.14


DIFFS = []
for i in range(len(filename)):
  if (1):
    FILE = extract(filename[i])

    #r_, v_, u_, m_, t_ = set_variables(FILE)
    r, v, u, m, t = set_variables(FILE)


    #r, v, u, m, t = remove_bound_particles(r_, v_, u_, m_, t_)

    speed = la.norm(v, axis=1).tolist()
    dist  = la.norm(r, axis=1).tolist()



    #plt.xscale('log')
    #plt.yscale('log')
    #"""

    predicted_speed = [ np.sqrt(2*(u[k] - 1 + MASS/dist[k])) for k in range(len(dist))]

    predicted_dist  = [MASS/(0.5*speed[k]**2 +1- u[k]) for k in range(len(speed))]
    #to do next: compare corrected energies to initial energies AND u, the given energy

    computed_u = [1 + 0.5*speed[k]**2 - MASS/dist[k] for k in range(len(dist))]

    u_diff = [(u[k]-u[0])/u[0] for k in range(len(dist)) ]
    speed_diff = [(speed[k] - predicted_speed[k]) for k in range(len(speed)) ]
    dist_diff = [(dist[k] - predicted_dist[k]) for k in range(len(dist)) ]
    DIFFS.append(dist_diff)
    #B = 10**(B)

    points = 2
    pS  = [predicted_speed[x-points] for x in range(points)]
    S  = [speed[x-points] for x in range(points)]
    T  = [t[x-points] for x in range(points)]


    #plt.plot(t, speed_diff)
    #plt.plot(T,pS,label= "predicted Speed")
    #plt.plot(T,S, label= "Speed")

    i = -1
    v_t = la.norm(v[i])
    derr_n = 24*MASS*v_t
    derr_d = (v_t**2 -  2*(u[i]-1)   )**2
    derr   = derr_n/derr_d*(5e-7)
    fudge = derr/6
    MAX = max(dist_diff)

    delta_u = np.array(u)-u[0]

    #plt.xlim(.8, 1)
    #plt.xlim(t[-points], t[4-points])
    #plt.plot(t,speed_diff, alpha = .5)#, label= "dist differences")
    plt.plot(t,delta_u, alpha = .5)#, label= "computed distance")
    #plt.plot(t,dist, alpha = .5, label= "computed distance")
    plt.title(" Evolution of unbound particles")
    plt.xlabel("Time (s)")
    #plt.ylabel("Simulated distance minus predicted distance")


plt.legend()
plt.show()
#we find that velocities through the simulation are greater than velocities that are predicted.

#difference between simulated total energy and predicted energy dampens

#difference between simulated distances and predicted distance


#plt.xlim(0,2000)
#plt.ylim(.9999,1.005)
#"""

#plot_avg_velocity_vs_polar_angle(vr, m, phi, B)


"""
#plot the distribution of energies for those with velocities below the cutoff

plt.hist(velocities, weights=masses, bins = num_of_bins, alpha=.5, label=label)
"""
#want to know how much mass is at each energy, for those above v_cut, and below 

#then we can also plot mass-avged u 




#plot_2dhist(phi, vr, m, num_of_bins = 150)




# remove saving feature from code, see if that speeds up evolution
#plot Energy(t) for unbound particles
#if bad, decrease time step. do adaptive time step
#for adaptive, make step proportional to distance from BH. Have all particles update with large dt,
#but closer ones have smaller updates

#plot mass on y axis vs ye 






