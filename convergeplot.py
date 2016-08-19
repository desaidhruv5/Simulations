import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pylab import *
from matplotlib.colors import LogNorm

plt.style.use('ggplot')
plt.rc('font', family='serif')
pi=np.pi

"""#dt   = 0.628319
#error= 0.107401

dt =   [0.628319, 0.314159,    0.125664,   0.0628319,   0.0314159,    0.0125664,   0.00628319]
error = [0.107401, 0.00390535,  6.15661e-05, 3.0481e-06,  1.65412e-07,   3.84832e-09,   2.32465e-10]

plt.xscale('log')
plt.yscale('log')

dt = [elem * 4.96*1e-6 for elem in dt]

log_dt = np.log10(dt)
log_error = np.log10(error)




plt.plot(x,y, label='slope = 4 (theoretical RK4 Global order of convergence)')

plt.plot(dt, error, linewidth=1)
plt.xlabel(r'$\Delta t$ [seconds]')
plt.ylabel('Error in radius')
plt.title('Global Error vs. time step')
plt.legend()
plt.show()


print (-5 - -9)/(-1 - -2)

"""





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
def set_variables(FILE):
  
  dt = [line[0] for line in FILE]
  error = [line[1] for line in FILE]

  
  return dt, error

#-----------------------------------------------





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




file = extract("convergeplot.dat")

time_step, relative_error = set_variables(file)

plt.xscale('log')
plt.yscale('log')
plt.plot(time_step, relative_error, label = "actual error")
x=[1e-4, 1e-3]
y=[1e-9,1e-5]
plt.plot(x,y, label='slope = 4 (theoretical RK4 Global order of convergence)')

plt.xlabel(' time step [code units]')
plt.ylabel('Relative error')
plt.title(r'Global computational error')
  


plt.legend(loc = 'lower right')
plt.show()

#figure out amount of time it takes to loop through N particles using method 1 and method 2.


#plot mass on y axis vs ye 






