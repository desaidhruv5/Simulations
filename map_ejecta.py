import matplotlib.pyplot as plt
import numpy as np
import math
import os
from matplotlib.pylab import *
from matplotlib.colors import LogNorm
from numpy import linalg as la
from numpy import matrix as mat
from numpy.linalg import svd as svd
from mpl_toolkits.mplot3d import Axes3D

plt.style.use('ggplot')
plt.rc('font', family = 'serif')


# REMEMBER TO DO THIS
#cut -d " " -f 2- solidangle.dat > 1second_later.dat



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
  data = []
  with open(filename, "r") as infile: #reading file of all masses with respective fallback times
    for line in infile:
      if line.startswith("#"):
        continue
      else:
        data.append(map(float, line.split()))
  return data




def set_variables(DATA):
  print "Extracting data..."
  #position components
  x = [1.48*line[1] for line in DATA]
  y = [1.48*line[2] for line in DATA]
  z = [1.48*line[3] for line in DATA]


  #radial vector
  r = [ [x[k], y[k], z[k]] for k in range(len(x))]  
  


  #velocity components
  vx = [line[4] for line in DATA]
  vy = [line[5] for line in DATA]
  vz = [line[6] for line in DATA]
  v = [[vx[i], vy[i], vz[i]] for i in range(len(vx))  ]


  #energy
  energy   = [line[7] for line in DATA]      #actually specific energy. units of c^2
  

  mass = [line[0]*1.99e33 for line in DATA]   #grams


  print "Variables set."
  return r, v, energy, mass


def remove_bound_particles(r_, v_, energy_, mass_):

  print "Removing bound particles..."
  energy           = []
  mass             = []
  uncounted_mass   = []
  false_energy     = []
  r                = []
  v                = []
  #energy_ = (1+.5*la.norm(v_,axis=1)**2-MASS/la.norm(r_,axis=1)).tolist()

  for j in range(len(mass_)):
    if  energy_[j]>1:
      mass.append(mass_[j])
      v.append(v_[j])
      energy.append(energy_[j])
      r.append(r_[j])
  #vr = np.delete(vr, 0, 1)


  #print "should have been deleted: ", count
  #print "New radial_velocity.shape, ", vr.shape
  #print "should agree with mass length: ", len(mass)
    

  """
  print "fraction of mass left out: ", sum(uncounted_mass)/sum(mass)
  print "fraction of mass with u>1, vr<0 (falsely marked): ", sum(false_energy)/sum(mass)
  """


  print "Bound particles removed."
  return r, v, energy, mass


def pick_most_massive_marticles(velocity, mass):
  pass
  MED = np.percentile(mass, 90)
  m   = []
  v   = []
  #print vr.shape
  #print mat(radial_velocity[:,0]).T.shape
  print "90th percentile mass: ", MED

  for j in range(len(mass)):
    if mass[j] > MED:
      v.append(velocity[j])
      m.append(mass[j])
  return v, m


def calculate_U(velocity, mass):


  v, m = pick_most_massive_marticles(velocity, mass)
  #w_v = np.delete(w_v, 0, 1)    #deleting first column


  v= np.array(v)          #Nx3
  w_v = np.array(v).T*m   #3xN



  c = w_v.sum(axis=1)/sum(mass)
  #print "approx offset: ", c
  M = mat(w_v)

  print "commence svd..."
  U, s, V = svd(M)
  print "svd complete."
  

  return U
  #return c, k



def transform_coordinates(velocity, mass):


  v = np.array(velocity)  #Nx3 array
  weighted_v = v.T*mass   #3xN array
  offset = mat(weighted_v.sum(axis=1)).T/sum(mass)
  U = calculate_U(velocity, mass)

  M = mat(v.T)
  new_v = (U.T*M).T
  print "U"
  print U

  return new_v




def find_angles(vr):

  vr = mat(vr).tolist()
  print "len of vr ", len(vr)

  def polar(v):
    return math.atan2(v[1],v[0])


  polar_angle = [polar(v) for v in vr]   #radians

  cos_theta = [v[2]/la.norm(v) for v in vr]

  return polar_angle, cos_theta



def rotate(a, v):
  #a = angle to rotate by
  R = [[np.cos(a), np.sin(a), 0] ,[-np.sin(a), np.cos(a), 0], [0,0,1]]
  R = mat(R)

  rot_v = (R*mat(v).T )
  rot_v = (rot_v.T).tolist()

  return rot_v

def plot_v_map(phi, v, m, num_of_bins):
  speed = la.norm(v, axis=1)
  counts, xedges, yedges, Image = hist2d(phi, speed, bins= num_of_bins, weights=m, norm=LogNorm())
  plt.colorbar().set_label(label=r'Mass [$g$]')
  plt.xlabel(r'$\phi$ [angle from $x$-axis]')
  plt.ylabel(r'Speed $[km/s]$')
  plt.title(r'%s: Distribution of mass over velocity & polar angle' %NAME)
  plt.savefig('%s_mass_distribution_vr_phi.png' %NAME)
  #hist, bin_edges, ignored = plt.hist(
  plt.show()  

  x = [(xedges[i+1]+xedges[i])/2 for i in range(len(xedges)-1)]
  y = [(yedges[i+1]+yedges[i])/2 for i in range(len(yedges)-1)]
  return x, y, counts


def plot_solidangle(phi, cos, m, num_of_bins):

  counts, xedges, yedges, Image = hist2d(phi, cos, bins= num_of_bins, weights=m, norm=LogNorm())
  plt.colorbar().set_label(label=r'Mass [$g$]')
  plt.xlabel(r'$\phi$ (angle from $x$-axis in radians)')
  plt.ylabel(r'$\cos(\theta)$ ($\theta$ is angle from $z$ axis)')
  plt.title('%s: Directions of velocities of the unbound particles'%NAME)
  plt.savefig('%s_solidangle.png' %NAME)
  plt.show()

  x = [(xedges[i+1]+xedges[i])/2 for i in range(len(xedges)-1)]
  y = [(yedges[i+1]+yedges[i])/2 for i in range(len(yedges)-1)]
  return x, y, counts


def output_data(FILENAME, x, y, mass):
  #np.savetxt(FILENAME, np.transpose([mass, x,y]), fmt='%1.6g')

  file = open(FILENAME, "w")
  for i in range(len(x)):
    for j in range(len(y)):
      file.write("%1.6e " %mass[i][j]) #bin total mass
      file.write("%1.6e " %x[i]) #logr coordinates
      file.write("%1.6e " %y[j]) #r/z coordinates
      file.write("\n")


#---------------------------------------------

#USER INPUT#############
B=150      #number of bins
NAME= "M5_S7"
MASS = 6.11
########EDIT LINE IMMEDIATELY BELOW THIS


filename = check_file_existence("solidangle.dat")

data = extract(filename)

r_, v_, u_, m_ = set_variables(data)


r, v, u, m = remove_bound_particles(r_, v_, u_, m_)
m = [elem*100./(4.*math.pi) for elem in m]

v1= transform_coordinates(v, m)

#v2= transform_coordinates(v1, m)
#v3= transform_coordinates(v2, m)
#print v_new.shape
#print vr.type



#v1 = rotate(np.pi, v1)
phi, cos = find_angles(v1)

print "Behold!"
binphi, bincos, binmass = plot_solidangle(phi, cos, m, B)
output_data('%s_binned_solidangle.dat' %NAME, binphi, bincos, binmass)

binphi,  binv1, binmass = plot_v_map(phi, v1, m, B)
output_data('%s_binned_velocity_map.dat' %NAME, binphi, binv1, binmass)



#scaling, so that we get mass/steradian,   for each square

"""
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(np.array(v1).T[0].tolist(), np.array(v1).T[1].tolist(), np.array(v1).T[2].tolist() ,norm=LogNorm() )#, color=norm(m))#, weight= m)
"""


