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

  """  
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
  """
  energy = [line[7] for line in DATA]		#actually specific energy. units of c^2


  t = [line[12] for line in DATA]
  t = [x +TIMEAFTERMERGER for x in t]
  t = [x/2.27e5 for x in t]

  #mass = [line[0]*1.99e33 for line in DATA]	#grams
  mass = [line[0] for line in DATA]	#grams

  return mass, t, energy







#####################################################

#BEGIN FALLBACK CALCULATIONS

###########################################
"""
logt = np.log10(t)

B = np.linspace(logt[0], logt[-1], 81)
dm, bins, ignored = plt.hist(logt, bins = B, weights = m)
plt.show()
T  = []
dt = []
for i in range(len(B)-1):
  dt.append(10**B[i+1]-10**B[i])
  T.append( (B[i+1]+B[i])/2)

T= 10**np.array(T)
"""

def calc_fallback(mass, time, res = 200):
  
  print "calculating fallback rate..."

  log_t_start = np.log10(time[0])
  log_t_end = np.log10(time[-1])

  B = np.logspace(log_t_start,log_t_end, num=res)
  dm, bins, ignored = plt.hist(time, bins = B, weights = mass)
  plt.close()
  T  = []
  dt = []
  B_log = np.log10(B)
  for i in range(len(B)-1):
    dt.append(B[i+1]-B[i])
    T.append( 10**(  (B_log[i+1]+B_log[i])/2  )       )

  Mrate = np.array(dm)/np.array(dt)

  return Mrate, T

####################################################################


def plot_fallback(mdot, T):
  plt.xscale('log')
  plt.yscale('log')
  x=[1*10**(-1.7), 1*10**(-.5)]
  y=[ 1 , 1e-2]
  plt.plot(x, y, label = r'expected $t^{-5/3}$ slope')
  plt.plot(T, mdot, label = "Simulated fallback rate")
  plt.xlabel('time [seconds]')
  plt.ylabel(r'Fallback rate [$M_{\odot}/s$]')


def output_data(FILENAME, x, mass):
  #np.savetxt(FILENAME, np.transpose([mass, x,y]), fmt='%1.6g')

  file = open(FILENAME, "w")
  for i in range(len(x)):
    file.write("%1.6e " %mass[i]) #bin total mass
    file.write("%1.6e " %x[i]) #logr coordinates
    file.write("\n")


  ##USER EDIT############
  plt.title(r'%s: Mass fallback rate' %NAME)
  #'''
  plt.legend()
  plt.savefig('fallback.png')
  plt.show()  
#############################################################################









#USER INPUT##################################################################
NAME = "M5_S7"
TIMEAFTERMERGER = 346


data = extract("fallback.dat")
m, t, u = set_variables(data)


Mrate, T = calc_fallback(m, t)


plot_fallback(Mrate, T)

output_data('%s_binned_fallback.dat' %NAME, T, Mrate)
























