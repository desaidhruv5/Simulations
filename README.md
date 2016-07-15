#POST-PROCESSING CODE FOR BLACK HOLE-NEUTRON STAR MERGER SIMULATIONS


#_Background_

Data is taken from General Relativistic (GR) - Hydrodynamic Simulations
run by the Spectral Einstein Code (SpEC). ( Visit http://www.black-holes.org/SpEC.html).

These simulations are of binary systems which consist of a neutron star (NS) orbiting around a central black hole (BH). Under appropriate conditions, such as the ones in these simulations, the NS is tidally disrupted as it falls into the BH. Low-energy matter falls in, energetic matter effectively slingshots around the BH and is ejected from the system, and matter in between forms an accretion disk.

To learn more about the cool physics of these systems, visit http://www.black-holes.org/!

Data is taken at approximately ~10 ms (simulation time) into the simulation. This is near the beginning of the merger and during the formation of the accretion disk. Only particles that lie ~40 km away from the BH are used in this post-processing code. In this regime, we have deemed GR effects to be negligible and the Newtonian approximation of the gravitational potential to be sufficient.


The code further evolves particles typically for ~1 second, and is primarily concerned with (1) the rate at which bound matter falls back to the black hole during the evolution, and (2) the distribution of ejecta, or unbound matter, at the end of the evolution.


#_Main Evolution Code_

'density.C' is the main post-processing evolution code, which evolves particles in a Newtonian potential. Densities are approximated, and a heating scheme is also incorporated.

*IMPORTANT*:

+ When compiling, use optimization option '-02' to drastically reduce runtimes.

+ This code requires C++11 to run. Use option -std=c++11 when compiling.

+ The 'Eigen Library' is also required for the code to run. Once installed, place a symbolic link named, 'Eigen_lib', pointing to /usr/.../eigen/3.2.8/include/eigen3/, in the same directory.

INPUT FILE FOR EVOLUTION CODE:

The main evolution code requires the input file, 'density.input', to run. The following parameters can be changed in the input file:

* Step size
-- This is the step size of the simulation, which is constant throughout. This is in code units (1 second = 2e5 units of time in the code). Under 10 is recommended. Typically 2.

* Mass
-- This is the mass of the central black hole, in solar masses.

* Filename
-- This is the name of the file of initial data, taken from the SpEC simulation, which 'density.C' reads from. Format of the data should correspond to that described at the top of 'density.C'.

* Domain restrictions
-- To filter out any particles in the initial data based on their properties, this must be set to 1, and 'density.C' must be editted appropriately, inside the 'read' function. For no filter (to read in all particles), set to 0.

* Number of nearest neighbors, _n_
-- This is related to density calculations made throughout the simulation. Density is approximated at a handful of particles based on their proximities to _n_ nearest neighbors. Typically 10-20. (This number includes the central particle itself, at which density is approximated, so _n_ - 1 other particles are looked at.)

* Total evolution time
-- This is the total time that the simulation runs for, in code units (conversion in "Step size", above).

* Heating
-- This number turns heating on (= 1) or off (= 0). This effectively translates to a boost in velocity for particles, and applies mostly to those traveling outwards, away from the BH. Heating is due to radioactive decay of heavy elements formed by r-process nucleosynthesis.

OUTPUT FILES:

The code produces the following files, which are subsequently used for analysis with their respective Python scripts (detailed below):

* 'fallback.dat' -- used to measure fallback rate

* 'solidangle.dat' -- used to map ejecta distribution at the end of the simulation

* 'densevo*.dat' -- this set of files consists of the evolution of a handful of particles at which densities are approximated. 


* 'comp_time.dat' -- this outputs rate at which the simulation progresses


#_Python Scripts_

MEASURING THE RATE OF FALLBACK OF BOUND MATERIAL

Data is read from 'fallback.dat'.

'fallback.py' measures and plots the rate at which matter falls back on to the black hole. For all particles in bound orbits, A particle is considered to have fallen back if it reaches the point of closest approach in its orbit (pericenter).

'orbit.py' predicts the fallback time for particles that have not yet fallen back, via Kepler's equations and orbital mechanics.

===========

PLOTTING THE DISTRIBUTION OF EJECTA

Data is read from 'solidangle.dat'.

'mapejecta.py' plots two 2d histograms, where color represents mass: (1) a mercator projection of the directions of velocities of unbound matter; and (2) a profile of velocities as a function of polar angle. Since particles may not all be lying in the xy-plane, a best-fit plane is found (via Singular Value Decomposition [SVD]), and coordinates are rotated such that most of the matter thereafter lies in a new x'y'-plane. SVD calculations increase runtime to just under a minute.

'solidangle.py' has the same function as (1) of 'mapejecta.py', except that it does not transform coordinates and thus gives more immediate results.

===========

PLOTTING APPROXIMATE DENSITIES THROUGHOUT EVOLUTION

Data is read from all 'densevo*.dat' files.

'densities.py' plots the evolution of densities at a handful of tracer particles in the simulation, approximated by the proximity of their _n_ nearest neighbors. Density is calculated by taking the mass of _n_ nearest neighbors, and dividing it by the volume of the smallest sphere enclosing the _n_ particles, centered at the tracer particle.









