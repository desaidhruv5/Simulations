#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include "Eigen_lib/Eigen/Dense"
#include <cstdio>
#include <time.h>

//////////////////////////////////////////////////////////////////*

/*
1) The density is greater than 1.e-10 
2) The radius is larger than 40
A note on units: we set G=c=1 and solar mass = 1, hence
Density: 1 = 6.1e17 g/cc
Distance: 1 = 1.477 km
Velocity: 1 = 3e10 cm/s = c
Time: 1 = 4.92703806e-6 s
The first few lines give you the list of variables. 
"Minus_ut" is the specific energy (in case we want to use that instead of velocities). A priori, all particles with Minus_ut>1 should be unbound.
Ye is the electron fraction - ratio of electrons to protons+neutrons (or the fraction of baryons that are protons). 

File copied from:
scp -r briaree:/RQexec/ffoucart/   etc.

read_file       >> m >> x >> y >> z >> vx >> vy >> vz >> u >> rho >> temp >> ye >> e;

fallback_file   << m << x << y << z << vx << vy << vz << u << rho << temp << ye << e << t;


solid_angle     << id << m << x << y << z << vx << vy << vz << u << rho << temp << ye << e;
solid_angle     << m << x << y << z << vx << vy << vz << u << rho << temp << ye << e;

*/


//MASS OF BLACK HOLE
double MASS=1.e200;
const double G = 1.;
const double pi = 4.*atan(1);
int nnn=-1;// = 20; //Number of nearest neighbors


double EVOLTIME=-10000.;

int i_D_count = 20;
int Dcount = i_D_count; //Dcount is number of unbound particles being tracked


//If set to 1, domain will be restricted according to formula in read() function.
//If set to 0, bhsim will incorporate all particles.
int domain_restricted=-1;



double Heat=1.e300;    //if 1, add heating
double T_Heat=-1.e300;    //if 1, add heating




//checks direction of particle velocity
//returns 1 if radial velocity is positive, else returns 0
bool check (const Eigen::Vector3d & r, const Eigen::Vector3d & v)
{

  bool b;
  if (r.dot(v) > 0) b = 1;
  else b = 0;

  return b;
}













///////////////

class Particle
{      //this is a vector of particles
public:
  double m;           //mass of particle
  Eigen::Vector3d r;  //position vector of particle
  Eigen::Vector3d v;  //velocity vector of particle
  double u;           //energy
  double ye;          //electron fraction
  double e;           //entropy
  double rho;         //density
  double temp;        //temperature
  bool vcheck;        //1 if particle is moving outwards, otherwise 0
  int id;             //initial index of particle in original list
  bool dflag   = 0;   //(default is 0) set to 1 for particles that we track densities for
  double tapo  = 0;
  double dt    = 2;
  double t     = 0;
  

  Particle (double m_, Eigen::Vector3d r_, Eigen::Vector3d v_) :    //constructor for object of Particle class
               m(m_),               r(r_),              v(v_){}



  void
  print_position(std::ofstream & outf)
  {                              //report particle info
    outf << r[0] << " " << r[1]<< " "<< r[2] << " ";
  }
  
  void
  print_velocity(std::ofstream & outf)
  {                              //report particle info
    outf << v[0] << " " << v[1]<< " "<< v[2] << " ";
  }

  void
  print_quantities(std::ofstream & outf)
  {
    outf << u << " " << rho << " " << temp << " " << ye << " " << e << " ";
  }

};
///////////////////////////////////////////





//This prints out 1 particle info line to output file 
void save(Particle & p, std::ofstream & outf, const double & t)
{
  
  outf << p.m << " ";
  p.print_position(outf);
  p.print_velocity(outf);
  p.print_quantities(outf);
  outf <<t << " "<< std::endl;
}

//particledat >> m >> x >> y >> z >> vx >> vy >> vz >> u >> rho >> temp >> ye >> e;











//This is function which compares distances of two particles from that of densityevo
bool compdist (const Particle & p1, const Particle & p2, const Particle & center)
{
  
  return ( (p1.r-center.r).norm() < (p2.r-center.r).norm() );  //This returns 1 if r1-r < r2-r

}




//Make separate velocity correction function here
void correct_velocity(Particle & p)//, std::ofstream & outf)
{
  //now correct the velocities based on energy ( E = KE + PE )
  //outf << p.m << " "<< p.r[0] << " " << p.r[1] << " " << p.r[2] <<" "<< p.v[0]<< " " << p.v[1] << " " << p.v[2]<< " ";
  //outf << p.u << " ";
  double correctv = (   p.u - 1 + MASS/p.r.norm()   )*2;
  correctv = sqrt(correctv);
  p.v = p.v * correctv/p.v.norm();

  double new_energy = 1+ .5*p.v.norm()*p.v.norm()-MASS/p.r.norm();
  //outf << p.v[0]<< " " << p.v[1] << " " << p.v[2]<< " "<< new_energy << std::endl;
  //std::cout << "correct: " << correctv << std::endl;
}
































//FUNCTION TO READ IN PARTICLE DATA--------------------------------------------------------------
//-------------------------------------------------------------------------------------------------


//'read' takes in file, outputs vector of particles, a.k.a. 'tail'
//also corrects velocities based on energy
std::vector <Particle> read(std::ifstream & particledat)
{


  std::cout << "Reading in particles from tail..." << std::endl;
//STEPS:
//reads each line in file containing original list of particles
//1: creates object of type particle, 2:labels it, 3: assigns position, 4: assigns velocity, 5: ...

  //particle id - initial index in the particle list, which will stick with particle throughout evolution
  int id;
  //here we define what each element on each line represents
  double m, x, y, z, vx, vy, vz, u, rho, temp, ye, e;

  Eigen::Vector3d r;
  Eigen::Vector3d v;
  std::vector<Particle> tail;
  Particle p(m, r, v);

  int i = 0;

  //std::ofstream corrected("corrected.dat");
  while (!particledat.eof())
  {
        //peek at first character.   
    if (particledat.peek() == '#') particledat.ignore(256, '\n');  //does the line start with "#"? If so, ignore until new line
    
    else
    {
      particledat >> m >> x >> y >> z >> vx >> vy >> vz >> u >> rho >> temp >> ye >> e;
      
      //RESTRICTION ON DOMAIN
      //If domain_restricted=0, enter block.
      //If domain is restricted, but particle lies in subdomain, enter block. Else do nothing.

      double r = sqrt(x*x + y*y + z*z);
      double vnorm = sqrt(vx*vx + vy*vy + vz*vz);
      double phi = atan2(y,x);

      //(pi/-2 < phi && phi < pi/3 && r < 40 ) || ( phi > pi/3 || phi < pi/-2 )  

      if
      (
        ( domain_restricted != 1)    ||
        ( (u>1) &&   (sqrt(2*(u -1))< (8000/3e5))     )           //insert domain restriction here
      )
      
      {
        tail.push_back(p);              //creating entry for ith particle in the tail
        tail[i].m    = m;                  //mass of ith particle in the tail
        tail[i].r[0] = x;          //ith particle of the tail, position vector, create x component
        tail[i].r[1] = y;          //...
        tail[i].r[2] = z; 
        tail[i].v[0] = vx;         //ith particle of the tail, position vector, create x component
        tail[i].v[1] = vy;         //...
        tail[i].v[2] = vz;
        tail[i].u    = u;
        tail[i].ye   = ye;               //tail[i] is ith particle, ye is composition of particle
        tail[i].e    = e;
        tail[i].rho  = rho;
        tail[i].temp = temp;
        tail[i].id   = i;
        //std::cout << "id: " << tail[i].id << std::endl;
        tail[i].vcheck = check (tail[i].r, tail[i].v); //checking if intial radial velocity is positive/negative

        //now correct the velocities based on energy ( E = KE + PE )
        correct_velocity(tail[i]);//, corrected);

        //-pi/4 to pi/4, r > 70


        //std::cout << "actual: " << tail[i].v.norm() << std::endl;


        //std::cout << "i=" << i << std::endl;
        //std::cout << "m=" << tail[i].m << std::endl;
        int mod = 10000;
        if (i % mod==0) std::cout << "Approximate number of particles read: "<< i << "\r"<<std::flush;
        ++i;
      }
      
    }
    
  }
  std::cout << "There are " << tail.size()<<" particles." << std::endl;
  return tail;
}

//__________________________________________________________________________
//__________________________________________________________________________




















//COMPUTING FIELD at position r away from black hole, which is at origin
Eigen::Vector3d fieldon(const Particle & p)
{
  double R = p.r.norm();

  //know mass, position
  return p.r*(-1*G*MASS/(R*R*R));
}





//This updates 1 particle via RK 4th order method
void RKupdate(Particle & p)
{ 
                                              //after a time step of 'dt'
  double dt = p.dt;
  Eigen::Vector3d dr, dv;   //declare changes in position/velocity. To be added to present state variables
  Particle p_ = p;

  //declaring RK variables k1, k2, etc.
  Eigen::Vector3d kv1, kv2, kv3, kv4;   //velocities
  Eigen::Vector3d ka1, ka2, ka3, ka4;   //accelerations

  kv1  = p.v;
  ka1  = fieldon(p_);

  p_.r = p.r + .5*dt*kv1;                       //approximates field at the 
  kv2  = p.v + .5*dt*ka1;    //velocity at midpt using field initial
  ka2  = fieldon(p_);     //field at midpt using velocity initial

  p_.r = p.r + .5*dt*kv2;
  kv3  = p.v + .5*dt*ka2;    //velocity at midpt. using field 
  ka3  = fieldon(p_);

  p_.r = p.r + dt*kv3;
  kv4  = p.v + dt*ka3;
  ka4  = fieldon(p_);


  dv = dt/6.* ( ka1 + 2*ka2  + 2*ka3  + ka4  );
  dr = dt/6.* ( kv1 + 2*kv2  + 2*kv3  + kv4  );


  //update this info back into particle info:
  p.r = p.r + dr; 
  p.v = p.v + dv;
  p.t = p.t + dt;

}









//choosing particles to be flagged, to measure density at
//outputs list of indices
void set_dflags(std::vector<Particle> & tail)
{
  std::cout << "Finding candidate particles where density can be approximated." << std::endl;
  std::string densityfile;
  std::ofstream densitydat(densityfile);

  //FIND ALL PARTICLES SATISFYING REQUIREMENT
  std::vector<Particle> unbound;
  for (int j = 0; j < int(tail.size()); ++j)
  {
    if (tail[j].u > 1 && tail[j].r.norm() > 100/1.48)
    {
      unbound.push_back(tail[j]);
    }
  }
  if (Dcount> int(unbound.size() ) ) Dcount= unbound.size();
  std::cout << unbound.size() << " particles meet tracer requirement. Choosing "<< Dcount<<" particles." << std::endl;

  int j = 0;  //number of tracers found
  //TAKE PARTICLES EVENLY SPACED FROM ABOVE LIST
  for (int num_of_tracers_saved = 0; num_of_tracers_saved < Dcount; ++num_of_tracers_saved)
  {

    int ID = unbound[j].id;
    //std::cout<< "energy of particle stored "<< tail[dindex].u << std::endl;
    tail[ID].dflag= 1;

    densityfile = "densevo"+ std::to_string(tail[ID].id) + ".dat";
    densitydat.open(densityfile);
    densitydat.close();

    std::cout << "Distance of particle " << tail[ID].id<< ": " << tail[ID].r.norm()*1.48 << " km." << std::endl;
    //std::cout << "Should be " << 100/1.48<<std::endl;
    //std::cout<< "Creating new tracker particle: "<< std::endl;
    j = j + int(unbound.size())/Dcount;

  }
}





//we want it to only calculate density if it is also a tracer particle


void density_calculation(std::vector<Particle> & tail)
{
  std::string densityfile;
  std::ofstream densitydat(densityfile);
  //std::cout << "begin denscalc at time t = " << t << std::endl;
  int i = 0;
  int j = 0;

  while (i < Dcount)    //loop through all particles, only increment i when we find tracer particle
  {
    if (tail[j].dflag == 1)     //make sure it's a tracer particle && make sure it hasn't fallen back
    {
      Particle center = tail[j];

      //std::cout << "dflag: "<< tail[j].dflag <<", for tail[j].id: "<< tail[j].id << std::endl;
      
      //std::cout << "i: "<< i <<", j: "<< j << ", and tail id of flagged particle: " << center.id << std::endl;
      densityfile = "densevo"+ std::to_string(tail[j].id) + ".dat"; //setting correct file for info to be stored
      densitydat.open(densityfile, std::ofstream::app);
      //once Neighbors is filled with N-1 particles (in order, closest to farthest),
      //check if the next new tail particle is closer than Nth farthest particle (last one in list)
      //If not, move on to next tail particle

      //std::cout << "number of nearest neighbors= " << nnn << std::endl;
      //std::cout << "last tail id"<< tail[tail.size()-1].id << std::endl;

      //define comparison function which is modified each time depending on centered particle

      auto compare = [center] (Particle part1, Particle part2)
      {
        return compdist(part1, part2, center);
      };

      double NNN = fmin(nnn, tail.size());
      //making copy of tail and sorting that, so we can maintain order in original tail
      std::vector<Particle> ordered_tail = tail;
      std::partial_sort(ordered_tail.begin(), ordered_tail.begin()+NNN, ordered_tail.end(), compare);
      //std::cout << "after sort, last tail id" << tail[tail.size()-1].id << std::endl;

      //std::cout << "after sort id" << tail[0].id << " "<<(densityevo[denseinc].r- tail[0].r).norm() << std::endl;
      //std::cout << "after sort id" << tail[nnn-1].id << " "<<(densityevo[denseinc].r- tail[nnn-1].r).norm() << std::endl;

      //Calculating radius of sphere
      double rad = (tail[j].r- ordered_tail[NNN-1].r ).norm(); //distance of particle farthest away from center


      //Volume of sphere
      double vol = 4*pi*rad*rad*rad/3.;
      

      //Calculating total mass in sphere
      double mtot = 0;
      for (int k = 0; k < NNN-1; k++)
      {
        mtot = mtot + ordered_tail[k].m;
        //std::cout << "id: " << ordered_tail[k].id<< std::endl;
      }


      //Filling density into tracker particle list (densityevo)
      tail[j].rho = mtot/vol;

      /*
      if (t!=0)
      {
        tail[j].u = 1 + .5*tail[j].v.norm()*tail[j].v.norm() - MASS/tail[j].r.norm();
      }
      */

      save(tail[j], densitydat, tail[j].t);
      densitydat.close();
      ++i;
    }
    ++j;
  }

  //std::cout << "End denscalc at time t = " << t << std::endl;

}   //end looping over particles






//OUTPUTS RATE OF COMPUTATION PLOT
void computation_time(std::string str, clock_t initial, std::ofstream & outf, double t)
{
  clock_t current = clock();
  float seconds = (float(current-initial))/(CLOCKS_PER_SEC);
  std::cout << "Percent complete: "<< str <<", Minutes elapsed: "<< seconds/60 <<" Simulation time (s): "<< t/2.27e5<< "____ \r"<<std::flush;
  outf << str <<" "<<(float(current-initial))/(CLOCKS_PER_SEC) <<std::endl;
}











//Heating scheme:
const double X_s = 0.1;     //seed nuclei mass fraction
const double X_n = 1-X_s;   //neutron mass fraction

const double B_s = 8.7;     //in MeV, binding energy for seed nuclei
const double B_r = 8;       //in MeV, binding energy for r-process nuclei
const double diff = 1.293;  //in MeV, this is neutron proton rest mass difference
const double neutrino_f = .5;        //fraction of energy lost to neutrinos

const double proton_mass = 1.673e-24; //atomic mass unit to grams

//possibly variable in future. 2.03e5 code time in a second.
const double erg_in_joules = 1e-7;   //1erg = 1e-7 joule


void add_heat(Particle & p)
{
  //const double E_r = (1-neutrino_f)*(B_r - X_s*B_s - X_n*diff);

  const double E_r = Heat;
  const double total_E   = E_r * 1.602e-13;    //joules per nucleon
  const double one_sec = 2.27e5;  //1s = 2.27e5 codetime units
  double t_tot = T_Heat*one_sec;  //total heating time in code time
  //const double time_step = 2;
  double heat_rate = total_E/t_tot;

  //convert back to cgs (c=1)
  double delta_u = heat_rate/proton_mass/1e-3/(3e8*3e8);
  /*
  double t_calc;
  if (p.tapo !=0 && p.t< one_sec) t_calc = (one_sec+p.t-p.tapo+p.dt/2)/one_sec;
  else t_calc = (p.t+p.dt/2)/one_sec;
  //heating rate at a point halfway in the interval [t, t+dt], in MeV
  //double specific_e_dot = (Heat/2.6)*2e18*erg_in_joules*1e3*pow(  (.5 - 1/pi*atan((t_calc-1.3)/.11)  ),1.3)/(one_sec*3e8*3e8);
  //above is code energy per kg per codetime.

  p.v = sqrt(p.v.norm()*p.v.norm() + 2*specific_e_dot*p.dt)/p.v.norm()*p.v;
  */


  if (p.t < t_tot)
  {
    //std::cout << "Delta_u: "<<delta_u <<std::endl;
    double utest = p.u;
    p.v = (sqrt(p.v.norm()*p.v.norm() + 2*delta_u*p.dt) )/p.v.norm()*p.v;


  }


}





//MAIN EVOLUTION FUNCTION---------------------------------------------------------------
//--------------------------------------------------------------------------------------
//Evolves all particles, automatically stops each particle after reaching periapse, or at end of simulation
//also evolves density







bool generate_smaller_tail(std::vector<Particle> & tail, double tol, double dcalc, std::ofstream & tofile)
{
  std::vector<Particle> next_tail;

  //number of deleted particles
  int deleted  = 0;

  //set to 1 if all particles are caught up to same time and ready to calculate densities. else, 0.
  bool advance_densities = 1;
  for (size_t n = 0; n < tail.size(); n++)
  {

    if (tail[n].t < dcalc)          //if particle not reached density checkpoint,
    {


      //CALCULATE TIME STEP
      tail[n].dt = tol*2*M_PI*pow(tail[n].r.norm(), 1.5)/sqrt(MASS);  //calculate minimum time step for given tolerance
      if (tail[n].dt + tail[n].t > dcalc)                                 //if time step + current time exceeds final time,
      {
        tail[n].dt = dcalc- tail[n].t;                                    //force time step to be such that particle ends at final time
        //std::cout << "tail[n].dt: "<< tail[n].dt << std::endl;
      }
      else
      {
        advance_densities = 0;                               //else, still more evolution to go.

      }
      ///////////////////////



      if (check(tail[n].r, tail[n].v) == 1)       //currently moving outwards
      {

        if (tail[n].vcheck == 0)                  //used to move inwards
        {
          
          if (tail[n].dflag == 1)                 //if density tracker is about to be deleted,
          {
            Dcount = Dcount-1;
            std::cout << "Initially tracked "<< i_D_count <<" particles. Now down to "<< Dcount << std::endl; 
            tail[n].dflag=0;     
          }
          save(tail[n], tofile, tail[n].t);       //save particle data to file, along with time halted
          
          ++deleted;
        }

        else                                      // if moving outwards, keep updating
        {
          //heat
          add_heat(tail[n]);
          RKupdate(tail[n]);
          next_tail.push_back(tail[n]);
        }
      }




      else //if particle is MOVING INWARDS,
      {
        //save fallback time for each particle. 

        /*
        //enter if JUST REACHED APOAPSE and beginning to move inwards
        if (tail[n].vcheck == 1 && t>2.27e5)
        {

          if (t<1)
          {
            tail[n].tapo = t;
            add_heat(tail[n]);
          }
        }
        */

        tail[n].vcheck = 0;       //affirms that the particle _has_ in the past had negative radial velocity
        RKupdate(tail[n]);    //update for the case of falling particle
        next_tail.push_back(tail[n]);
      }
      //end the stopping condition loop
    }


    else
    {
      next_tail.push_back(tail[n]);
    }

  }   //end looping over particles




  /*
  double MIN = tail[0].dt;
  double MAX = tail[0].dt;  

  for (size_t n = 0; n < tail.size()-1; n++)
  {

    if (tail[n+1].dt <MIN)
    {
      MIN = tail[n+1].dt;
    }
    if (tail[n+1].dt >MAX)
    {
      MAX = tail[n+1].dt;
    }
     
  }
  */


  tail = next_tail;


  //std::cout << "MIN: " <<MIN<< " MAX: "<< MAX << std::endl;
  return advance_densities;
}















bool delete_particles(std::vector<Particle> & tail, double tol, double dcalc, std::ofstream & tofile)
{

  //number of deleted particles
  int deleted  = 0;

  //set to 1 if all particles are caught up to same time and ready to calculate densities. else, 0.
  bool advance_densities = 1;

  for (size_t n = 0; n < tail.size(); n++)
  {
    if (tail[n].t < dcalc)
    {

      tail[n].dt = tol*2*M_PI*pow(tail[n].r.norm(), 1.5)/sqrt(MASS);  //calculate minimum time step for given tolerance
      if (tail[n].dt + tail[n].t > dcalc)                                 //if time step + current time exceeds final time,
      {
        tail[n].dt = dcalc- tail[n].t;                                    //force time step to be such that particle ends at final time
        //std::cout << "tail[n].dt: "<< tail[n].dt << std::endl;
      }
      else
      {
        advance_densities = 0;                               //else, still more evolution to go.

      }





  //enter if CURRENTLY MOVING OUTWARDS (radial velocity is positive)
      if (check(tail[n].r, tail[n].v) == 1)
      {

        //USED TO MOVE INWARDS (enter this case only if velocity ever used to be negative)
        if (tail[n].vcheck == 0)
        {
          
          if (tail[n].dflag == 1)
          {
            Dcount = Dcount-1;
            std::cout << "Initially tracked "<< i_D_count <<" particles. Now down to "<< Dcount << std::endl; 
            tail[n].dflag=0;     
          }
          save(tail[n], tofile, tail[n].t);                   //save particle data to file, along with time halted
          
          tail.erase(tail.begin() + n);   //when a particle is deleted
          n = n - 1; 
          ++deleted;
          //std::cout << "p" << n+1 << " was deleted" << std::endl;
        }

        //enter if within range but moving outwards, keep updating
        else
        {
          //heat
          add_heat(tail[n]);
          RKupdate(tail[n]);
          //next_tail.push_back(tail[n]);
        }
      }



      //enter if particle is MOVING INWARDS
      else
      {
        //save fallback time for each particle. 

        /*
        //enter if JUST REACHED APOAPSE and beginning to move inwards
        if (tail[n].vcheck == 1 && t>2.27e5)
        {

          if (t<1)
          {
            tail[n].tapo = t;
            add_heat(tail[n], dt, t);
          }
        }
        */

        tail[n].vcheck = 0;       //affirms that the particle _has_ in the past had negative radial velocity
        RKupdate(tail[n]);    //update for the case of falling particle
        //next_tail.push_back(tail[n]);
      }
      //end the stopping condition loop
    }

  }   //end looping over particles




  //std::cout << deleted<< " particles were deleted." << std::endl;
  return advance_densities;
}















std::vector<Particle> autevolve(std::vector<Particle> & tail, const double & tol, std::ofstream & tofile, const double T = 2e5)
{ 

  //CANNOT CALCULATE DENSITY WHILE
  //LOOPING THROUGH ALL PARTICLES UPDATING POSITIONS, SINCE SOME PARTICLES WOULD BE UPDATED WHILE OTHERS WOULDN'T BE
  
  //SO MUST CALCULATE DENSITIES, THEN UPDATE POSITIONS
  //LOOP THROUGH tail. IF DCHECK =1, THEN CALCULATE DENSITY.


  double NUM = 100.;      //so we calculate densities/energies etc. 100 times
  std::vector<double> dtime;
  int exp = 0;
  while (pow(T, double(exp)/NUM) <= T)
  {

    dtime.push_back(pow(T, double(exp)/NUM));
    //std::cout << "dtime["<< exp <<"] = "<< dtime[exp]<< std::endl;


    ++exp;
  }


  int h_index = 0;
  while (dtime[h_index]<T_Heat*2.27e5)
  {
    ++h_index;
  }
  dtime.insert(dtime.begin()+h_index, T_Heat*2.27e5);




  std::ofstream comp_time("comp_time.dat");

  std::cout << "Finding particles to track..." << std::endl;
  set_dflags(tail);

  std::cout << "Calculating initial densities..." << std::endl;
  density_calculation(tail);


  std::cout << "Number of total density updates: " << dtime.size() << std::endl;  
  std::cout << "Commencing evolution... " << std::endl;
  clock_t initial = clock();





  double t;

  exp = 0;   //index for density calculations
  while (exp < int(dtime.size()) )
  {
    t = dtime[exp];   //catching all particles up to this time

  



    //UPDATE ENERGIES
    for (int j = 0; j < int(tail.size()); ++j)
    {
      tail[j].u = 1 + .5*tail[j].v.norm()*tail[j].v.norm() - MASS/tail[j].r.norm();
    }
    


    bool advance = generate_smaller_tail(tail, tol, t, tofile);
    
    if (advance==1)
    {
      density_calculation(tail);
      ++exp;


          //PROGRESS BAR  ///////////////////////////
      std::string percent_complete = std::to_string(t/T*100);
      computation_time(percent_complete, initial, comp_time, t);
      ///////////////////////////////////////////


    }
    /*
    std::cout << "exp: " << exp <<std::endl;
    std::cout << "tail.size(): " << tail.size() <<std::endl;
    */

    //faster if lots of particles being deleted
    /*
    if (t<1.6e4)
    {
      generate_smaller_tail(tail, t, dt, tofile);
    }
    //faster if not many particles being deleted
    

    else
    {
      delete_particles_same_tail(tail, t, dt, tofile);
    }
    */



    

    //CONSIDER SWITCHING TO LIST, AND STOP USING STD::VECTOR

    //std::cout << deleted<< " particles were deleted." << std::endl;



  }     //finish time whileloop




  clock_t final = clock();
  std::cout << std::endl;
  std::cout << "Complete. Minutes elapsed: "<<(final - initial)/(CLOCKS_PER_SEC)/60. << std::endl;

  return tail;

}
/////////////////////////////-----------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------











//This function TAKES OUTPUT LIST OF ONLY UNBOUND PARTICLES and calculates phi, cos(theta). Also prints masses and ye;
    //Creates output file to be used in solidangle.py
void solidangle(std::vector<Particle> & tail, std::ofstream & outf)
{

  for (size_t n = 0; n < tail.size(); n++)
  {
    outf << tail[n].m << " ";
    tail[n].print_position(outf);
    tail[n].print_velocity(outf);
    tail[n].print_quantities(outf);
    outf << std::endl;
  }
}

/////////////////////////////////////////////







//This will evolve all particles for a predetermined time, N*dt
void manevolve(std::vector<Particle> tail, double dt, double T = 2e5){

  double t = 0;
  int   i = 0;
  while (t < T)
  {    //begin time loop

    //std::cout <<"t"<< i << "\t";
    //tail[0].preport();              //report just the first particle at each time


//begin looping over all particles for given time
    for (size_t n = 0; n < tail.size(); n++)
    {
      RKupdate(tail[n]);    //update particle
    }

    t = t + dt;
    i = i + 1;
  }
}

//write code to track 

//Evolve should take in dt, tail vectorlist
//RKupdate should just take in one particle, dt



//////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////



int main() {

  double TOL;

  std::string TOL_;
  std::string MASS_;
  std::string particledat_;
  std::string domain_restricted_;
  std::string nnn_;
  std::string EVOLTIME_;
  std::string Heat_;
  std::string T_Heat_;

  std::string inputfile_name = "density.input";
  std::ifstream inputfile(inputfile_name);
  
  while (!inputfile.eof()) {         //while the end of file is not reached,
                                    //peek at first character.   
    if (inputfile.peek() == '#') inputfile.ignore(256, '\n');  //does the line start with "#"? If so, ignore until new line
      
    else { //reading bhsim.input and assigning each of its entries to appropriate objects
      std::getline(inputfile, TOL_);
      TOL = std::stod(TOL_);

      std::getline(inputfile, MASS_);
      MASS = std::stod(MASS_);

      std::getline(inputfile, particledat_);

      std::getline(inputfile, domain_restricted_);
      domain_restricted = std::stoi(domain_restricted_);

      std::getline(inputfile, nnn_);
      nnn = std::stoi(nnn_);

      std::getline(inputfile, EVOLTIME_);
      EVOLTIME = std::stod(EVOLTIME_);

      std::getline(inputfile, Heat_);
      Heat = std::stod(Heat_);            

      std::getline(inputfile, T_Heat_);
      T_Heat = std::stod(T_Heat_);

    }
  }

  std::cout << "Confirm the following:" << std::endl;
  std::cout << "Input file: "<< inputfile_name << std::endl;  
  std::cout << "Reading from: "<< particledat_<< std::endl;
  std::cout << "Tolerance: "<< TOL<< std::endl;
  std::cout << "Mass of black hole: "<< MASS<< std::endl;
  std::cout << "Is domain restricted? 1 = yes, 0 = no: "<< domain_restricted<< std::endl;
  std::cout << "Number of nearest neighbors per tracker particle: "<< nnn<< std::endl;
  std::cout << "Total evolution time: "<< EVOLTIME << std::endl;


  std::cout << "Heating: " <<Heat << " MeV/nucleon."<< std::endl;
  std::cout<< "T_Heat: "<<T_Heat<<std::endl;
  /*
  std::cout<< "one_sec: "<<one_sec<<std::endl;
  std::cout<< "T_Heat*one_sec: "<<T_Heat*one_sec<<std::endl;
  std::cout<< "t_tot: "<<t_tot<<std::endl;
  */


  std::ifstream particledat(particledat_);
  std::ofstream fallbackdat("fallback.dat");
  std::ofstream solidangledat("solidangle.dat");

  //USER INPUT: If 0, then domain unchanged. If 1, it's restricted according to read(). 

  if (domain_restricted==1) std::cout << "Confirm that domain is restricted as desired." << std::endl;
  std::cout << "Confirm that file is in the format:" << std::endl;
  std::cout << "read_file >> m >> x >> y >> z >> vx >> vy >> vz >> u >> rho >> temp >> ye >> e" << std::endl;






  std::vector<Particle> tail;
  tail  = read(particledat);

///////////////////

  fallbackdat << "# m << x << y << z << vx << vy << vz << u << rho << temp << ye << e << t" << std::endl;
  solidangledat << "# m << x << y << z << vx << vy << vz << u << rho << temp << ye << e" << std::endl;


  std::vector<Particle> ejecta;
  ejecta = autevolve(tail, TOL, fallbackdat, EVOLTIME);
  solidangle(ejecta, solidangledat);


  //manevolve(tail, dt);
  //std::cout << newtail[0].m << ", " << newtail[0].r << ", " << newtail[0].v << std::endl;


  return 0;
}





