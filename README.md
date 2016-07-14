Post-processing code for black hole-neutron star merger simulations
===========
BACKGROUND
Data is taken from General Relativistic (GR) - Hydrodynamic Simulations
run by the Spectral Einstein Code (SpEC). ( Visit http://www.black-holes.org/SpEC.html).

These simulations are of binary systems which consist of a neutron star (NS) orbiting around a central black hole (BH). Under appropriate conditions, such as the ones in these simulations, the NS is tidally disrupted as it falls into the BH. Low-energy matter falls in, energetic matter effectively slingshots around the BH and is ejected from the system, and matter in between forms an accretion disk.

To learn more cool physics regarding these systems, visit http://www.black-holes.org/!

Data is taken at approximately ~10 ms (simulation time) into the simulation. This is near the beginning of the merger and during the formation of the accretion disk. Only particles that lie ~40 km away from the BH are used in this post-processing code. In this regime, we have deemed GR effects to be negligible and the Newtonian approximation of the gravitational potential to be sufficient.

This code is primarily concerned with the rate of fall back of bound matter, and the distribution of ejecta, or unbound matter.

===========
MAIN EVOLUTION CODE

'density.C' is the main post-processing evolution code, which evolves particles in a Newtonian potential.

The code also relies on C++11. Use option -std=c++11 when compiling.

The code requires the input file, 'density.input', to run. Parameters are detailed below.



The following parameters can be changed via the input file:

Step size
-- This is the step size of the simulation, which is constant throughout. This is in code units (1 second = 2e5 units of time in the code). Under 10 is recommended. Typically 2.

Mass
-- This is the mass of the central black hole, in solar masses.

Filename
-- This is the name of the file of initial data, taken from the SpEC simulation, which 'density.C' reads from. Format of the data should correspond to that described at the top of 'density.C'.

Domain restrictions
-- To filter out any particles in the initial data based on their properties, this must be set to 1, and 'density.C' must be editted appropriately, inside the 'read' function. For no filter (to read in all particles), set to 0.

Number of nearest neighbors, _n_
-- This is related to density calculations made throughout the simulation. Density is approximated at a handful of particles based on their proximities to _n_ nearest neighbors. Typically 10-20. (This number includes the central particle itself, at which density is approximated, so _n_ - 1 other particles are looked at.)

Total evolution time
-- This is the total time that the simulation runs for, in code units (conversion in "Step size", above).

Heating
-- This number turns heating on (= 1) or off (= 0). This effectively translates to a boost in velocity for particles, and applies mostly to those traveling outwards, away from the BH. Heating is due to radioactive decay of heavy elements formed by r-process nucleosynthesis.


===========


















