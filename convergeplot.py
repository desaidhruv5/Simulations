import matplotlib.pyplot as plt
import numpy as np
plt.style.use('ggplot')
# plt.rc('text', usetex=True)
plt.rc('font', family='serif')

#dt   = 0.628319
#error= 0.107401

dt = 	[0.628319, 0.314159,		0.125664, 	0.0628319,	 0.0314159,    0.0125664,	 0.00628319]
error = [0.107401, 0.00390535,  6.15661e-05, 3.0481e-06,  1.65412e-07,	 3.84832e-09, 	2.32465e-10]

plt.xscale('log')
plt.yscale('log')

dt = [elem * 4.96*1e-6 for elem in dt]

log_dt = np.log10(dt)
log_error = np.log10(error)


x=[1e-7, 1e-6]
y=[1e-6,1e-2]

plt.plot(x,y, label='slope = 4 (theoretical RK4 order of convergence)')

plt.plot(dt, error, linewidth=1)
plt.xlabel(r'$\Delta t$ [seconds]')
plt.ylabel('Error in radius')
plt.title('Global Error vs. time step')
plt.legend()
plt.show()


print (-5 - -9)/(-1 - -2)