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
pi = np.pi

def check_file_existence(FILENAME):
  not_found = FILENAME+ " not found."
  found = FILENAME+ " found."
  for i in os.listdir(os.getcwd()):
    if i == FILENAME:
      not_found = found
      print found
      return FILENAME
      break

  if not_found != found:
    print not_found
    exit()

def extract(filename):
  DATA = []
  with open(filename, "r") as infile: #reading file of all masses with respective fallback times
    for line in infile:
      if line.startswith("#"):
        continue
      else:
        DATA.append(map(float, line.split()))
  return DATA


def set_variables(DATA):


  #position components
  #in km
  x = [line[1] for line in DATA]
  y = [line[2] for line in DATA]
  z = [line[3] for line in DATA]
  r = [ [x[k], y[k], z[k]] for k in range(len(x))] 
  
  #radial vector
  #radius = [np.sqrt(x[k]**2 + y[k]**2 + z[k]**2) for k in range(len(x))]
  
  #velocity components
  vx = [line[4] for line in DATA]
  vy = [line[5] for line in DATA]
  vz = [line[6] for line in DATA]
  v = [[vx[i], vy[i], vz[i]] for i in range(len(vx))  ]

  energy = [line[7] for line in DATA]   #actually specific energy. units of c^2
  
  
  #mass = [line[0]*1.99e33 for line in DATA]  #grams
  mass = [line[0] for line in DATA] #grams

  return mass, r, v, energy






def keep_bound_particles(mass_, r_, v_, energy_):

  print "Keeping only bound particles..."
  mass             = []

  r                = []
  v                = []
  #energy = (1+.5*la.norm(v_,axis=1)**2-MASS/la.norm(r_,axis=1)).tolist()
  energy =[]
  #print vr.shape
  #print mat(radial_velocity[:,0]).T.shape

  for j in range(len(mass_)):
    if  energy_[j]<1:
      mass.append(mass_[j])
      r.append(r_[j])
      v.append(v_[j])
      energy.append(energy_[j])


  print "Bound particles found."
  return mass, r, v, energy

#############################################################







def follow_single_particle(r,v):
  
  #convert to SI
  r = 1.477e3*np.array(r)
  v = 3e8*np.array(v)

  dist  = la.norm(r)
  speed = la.norm(v)

  #defining eccentricity vector
  e1 = (speed**2/mu - 1/dist)*r
  e2 = (dot(r,v)/mu)*v
  e = e1 - e2
  ecc   = la.norm(e)


  #defining true anomaly
  if dot(v,r) >0:
    true_anomaly = np.arccos(dot(e,r)/(dist*ecc))
  else:
    true_anomaly= 2*pi -np.arccos(dot(e,r)/(dist*ecc))
    #print "true anomaly: ", true_anomaly




  factor = np.sqrt((1+ecc)/(1-ecc))

  E = 2*np.arctan2(np.sin(true_anomaly/2), factor*np.cos(true_anomaly/2)  )

  atan_.append(E)

  if E <0:
    E = E + 2*2*pi
    print "An E value is negative. Corrected."


  #E is the eccentric anomaly

  #Kepler's equation
  M  = E - ecc*np.sin(E)
  
  delta_M = 2*pi - M
  #delta_M = - M


  a  = 1/(2/dist - (speed**2)/mu)
  n  = np.sqrt(mu/a**3)   #radians per second
  if 2*pi/n<1e-5:
    print "n: ", n
  fallback_time = delta_M/n

  return fallback_time


def generate_fallback_times(rlist, vlist):    #in seconds
  t_fall = range(len(rlist))
  
  for j in range(len(rlist)):
    t_fall[j] = follow_single_particle(rlist[j],vlist[j])

  return t_fall


def calc_fallback_rate(mass, time):
  


  log_t_start = np.log10(min(time))
  log_t_end = np.log10(max(time))

  B = np.logspace(log_t_start,log_t_end, num=200)
  dm, bins, ignored = plt.hist(time, bins = B, weights = mass)
  plt.close()
  T  = []
  dt = []
  B_log = np.log10(B)
  for i in range(len(B)-1):
    dt.append(B[i+1]-B[i])
    T.append( 10**(  (B_log[i+1]+B_log[i])/2  )       )

  mdot = np.array(dm)/np.array(dt)

  return T, mdot



def plot_fallback(T, mdot):
  plt.xscale('log')
  plt.yscale('log')
  plt.plot(T, mdot, label = "Simulated fallback rate")
  plt.xlabel('time [seconds]')
  plt.ylabel(r'Fallback rate [$M_{\odot}/s$]')


  ##USER EDIT############
  plt.title(r'M5_S9: Predicted Mass fallback rate')
  #'''
  plt.legend()
  plt.savefig('%s_orbit_fallback.png' %NAME)
  plt.show()  




#############################################################
#USE SI UNITS
MASS = 6.0552018     #solar mass in units of solar mass
G = 6.67408e-11       #Gravitational constant in SI
mu = MASS*G*1.989e30  #used for orbit calculations
NAME = "M5_S9"
TIMEAFTERMERGER= 346 + 30*2e5


filename = check_file_existence("solidangle.dat")
data = extract(filename)
m_, r_, v_, u_ = set_variables(data)

#m, r, v, u  = keep_bound_particles(m_, r_, v_, u_)


m, r, v  = m_, r_, v_

atan_=[]
t_ = generate_fallback_times(r,v)   #in seconds


t=[]
matter = []
r__ =[]
v__ =[]
atan =[]
for j in range(len(t_)):
  if t_[j]==t_[j]:
    t.append(t_[j])
    matter.append(m[j])
    r__.append(r[j])
    v__.append(v[j])
    atan.append(atan_[j])

"""
r__ = 1.48e3*np.array(r__)
v__ = 3e8*np.array(v__)
dist  = la.norm(r__, axis = 1)
speed = la.norm(v__, axis = 1)
a  = 1/(2/dist - (speed**2)/mu)
n  = np.sqrt(mu/a**3)
print 'min orbital period: ' , min(2*pi/n)

print "max ", max(t)
print "min ", min(t)
print "median ", np.median(t)

print "max atan ", max(atan)
print "min atan ", min(atan)
print "median atan", np.median(atan)
"""


t = [x +TIMEAFTERMERGER/2.27e5 for x in t]


T, Mrate = calc_fallback_rate(matter, t)


plot_fallback(T,Mrate)
#mass of BH = 6.022






"""
x = 152.6e9
y = 0
z = 0
vx = 0
vy = 29300
vz = 0

"""



"""
g = mu/dist**2
t_fall = sqrt(2*dist/g)

print "time to fallback under orbit: ", dt
print "time to fallback under constant initial acceleration:", t_fall

"""





