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
double MASS;
const double G = 1;
const double pi = 4*atan(1);
int nnn;// = 20; //Number of nearest neighbors


double EVOLTIME;

int i_D_count = 20;
int Dcount = i_D_count; //Dcount is number of unbound particles being tracked


//If set to 1, domain will be restricted according to formula in read() function.
//If set to 0, bhsim will incorporate all particles.
int domain_restricted;



int Heat;    //if 1, add heating




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
  bool dflag    = 0;  //(default is 0) set to 1 for particles that we track densities for
  

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
void correct_velocity(Particle & p, std::ofstream & outf)
{
  //now correct the velocities based on energy ( E = KE + PE )
  outf << p.m << " "<< p.r[0] << " " << p.r[1] << " " << p.r[2] <<" "<< p.v[0]<< " " << p.v[1] << " " << p.v[2]<< " ";
  outf << p.u << " ";
  double correctv = (   p.u - 1 + MASS/p.r.norm()   )*2;
  correctv = sqrt(correctv);
  p.v = p.v * correctv/p.v.norm();

  double new_energy = 1+ .5*p.v.norm()*p.v.norm()-MASS/p.r.norm();
  outf << p.v[0]<< " " << p.v[1] << " " << p.v[2]<< " "<< new_energy << std::endl;
  //std::cout << "correct: " << correctv << std::endl;
}
































//FUNCTION TO READ IN PARTICLE DATA--------------------------------------------------------------
//-------------------------------------------------------------------------------------------------


//'read' takes in file, outputs vector of particles, a.k.a. 'disk'
//also corrects velocities based on energy
std::vector <Particle> read(std::ifstream & particledat)
{


  std::cout << "Reading in particles from disk..." << std::endl;
//STEPS:
//reads each line in file containing original list of particles
//1: creates object of type particle, 2:labels it, 3: assigns position, 4: assigns velocity, 5: ...

  //particle id - initial index in the particle list, which will stick with particle throughout evolution
  int id;
  //here we define what each element on each line represents
  double m, x, y, z, vx, vy, vz, u, rho, temp, ye, e;

  Eigen::Vector3d r;
  Eigen::Vector3d v;
  std::vector<Particle> disk;
  Particle p(m, r, v);

  int i = 0;

  std::ofstream corrected("corrected.dat");
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
        disk.push_back(p);              //creating entry for ith particle in the disk
        disk[i].m    = m;                  //mass of ith particle in the disk
        disk[i].r[0] = x;          //ith particle of the disk, position vector, create x component
        disk[i].r[1] = y;          //...
        disk[i].r[2] = z; 
        disk[i].v[0] = vx;         //ith particle of the disk, position vector, create x component
        disk[i].v[1] = vy;         //...
        disk[i].v[2] = vz;
        disk[i].u    = u;
        disk[i].ye   = ye;               //disk[i] is ith particle, ye is composition of particle
        disk[i].e    = e;
        disk[i].rho  = rho;
        disk[i].temp = temp;
        disk[i].id   = i;
        //std::cout << "id: " << disk[i].id << std::endl;
        disk[i].vcheck = check (disk[i].r, disk[i].v); //checking if intial radial velocity is positive/negative

        //now correct the velocities based on energy ( E = KE + PE )
        correct_velocity(disk[i], corrected);

        //-pi/4 to pi/4, r > 70


        //std::cout << "actual: " << disk[i].v.norm() << std::endl;


        //std::cout << "i=" << i << std::endl;
        //std::cout << "m=" << disk[i].m << std::endl;
        int mod = 10000;
        if (i % mod==0) std::cout << "Approximate number of particles read: "<< i << "\r"<<std::flush;
        ++i;
      }
      
    }
    
  }
  std::cout << "There are " << disk.size()<<" particles." << std::endl;
  return disk;
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
void RKupdate(Particle & p, const double & dt) { 
                                              //after a time step of 'dt'

  Eigen::Vector3d dv;
  Eigen::Vector3d dr;
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


  dv = dt/6.* ( ka1 + 2*ka2  + 2*ka3  + ka4  );       //find dv
  dr = dt/6.* ( kv1 + 2*kv2  + 2*kv3  + kv4  );


  //update this info back into particle info:
  p.r = p.r + dr; 
  p.v = p.v + dv;

}









//choosing particles to be flagged, to measure density at
//outputs list of indices
void set_dflags(std::vector<Particle> & disk)
{
  std::cout << "Finding candidate particles where density can be approximated." << std::endl;
  std::string densityfile;
  std::ofstream densitydat(densityfile);

  //FIND ALL PARTICLES SATISFYING REQUIREMENT
  std::vector<Particle> unbound;
  for (int j = 0; j < int(disk.size()); ++j)
  {
    if (disk[j].u > 1 && disk[j].r.norm() > 100/1.48)
    {
      unbound.push_back(disk[j]);
    }
  }
  std::cout << unbound.size() << " particles meet tracer requirement. Choosing "<< Dcount<<" particles." << std::endl;

  int j = 0;  //number of tracers found
  //TAKE PARTICLES EVENLY SPACED FROM ABOVE LIST
  for (int num_of_tracers_saved = 0; num_of_tracers_saved < Dcount; ++num_of_tracers_saved)
  {

    int ID = unbound[j].id;
    //std::cout<< "energy of particle stored "<< disk[dindex].u << std::endl;
    disk[ID].dflag= 1;

    densityfile = "densevo"+ std::to_string(disk[ID].id) + ".dat";
    densitydat.open(densityfile);
    densitydat.close();

    std::cout << "Distance of particle " << disk[ID].id<< ": " << disk[ID].r.norm()*1.48 << " km." << std::endl;
    //std::cout << "Should be " << 100/1.48<<std::endl;
    //std::cout<< "Creating new tracker particle: "<< std::endl;
    j = j + int(unbound.size())/Dcount;

  }
}





//we want it to only calculate density if it is also a tracer particle


void density_calculation(std::vector<Particle> & disk, const double & t)
{
  std::string densityfile;
  std::ofstream densitydat(densityfile);
  //std::cout << "begin denscalc at time t = " << t << std::endl;
  int i = 0;
  int j = 0;

  while (i < Dcount)    //loop through all particles, only increment i when we find tracer particle
  {
    if (disk[j].dflag == 1)     //make sure it's a tracer particle && make sure it hasn't fallen back
    {
      Particle center = disk[j];

      //std::cout << "dflag: "<< disk[j].dflag <<", for disk[j].id: "<< disk[j].id << std::endl;
      
      //std::cout << "i: "<< i <<", j: "<< j << ", and disk id of flagged particle: " << center.id << std::endl;
      densityfile = "densevo"+ std::to_string(disk[j].id) + ".dat"; //setting correct file for info to be stored
      densitydat.open(densityfile, std::ofstream::app);
      //once Neighbors is filled with N-1 particles (in order, closest to farthest),
      //check if the next new disk particle is closer than Nth farthest particle (last one in list)
      //If not, move on to next disk particle

      //std::cout << "number of nearest neighbors= " << nnn << std::endl;
      //std::cout << "last disk id"<< disk[disk.size()-1].id << std::endl;

      //define comparison function which is modified each time depending on centered particle

      auto compare = [center] (Particle part1, Particle part2)
      {
        return compdist(part1, part2, center);
      };


      //making copy of disk and sorting that, so we can maintain order in original disk
      std::vector<Particle> ordered_disk = disk;
      std::partial_sort(ordered_disk.begin(), ordered_disk.begin()+nnn, ordered_disk.end(), compare);
      //std::cout << "after sort, last disk id" << disk[disk.size()-1].id << std::endl;

      //std::cout << "after sort id" << disk[0].id << " "<<(densityevo[denseinc].r- disk[0].r).norm() << std::endl;
      //std::cout << "after sort id" << disk[nnn-1].id << " "<<(densityevo[denseinc].r- disk[nnn-1].r).norm() << std::endl;

      //Calculating radius of sphere
      double rad = (disk[j].r- ordered_disk[nnn-1].r ).norm(); //distance of particle farthest away from center


      //Volume of sphere
      double vol = 4*pi*rad*rad*rad/3.;
      

      //Calculating total mass in sphere
      double mtot = 0;
      for (int k = 0; k < nnn-1; k++)
      {
        mtot = mtot + ordered_disk[k].m;
        //std::cout << "id: " << ordered_disk[k].id<< std::endl;
      }


      //Filling density into tracker particle list (densityevo)
      disk[j].rho = mtot/vol;

      /*
      if (t!=0)
      {
        disk[j].u = 1 + .5*disk[j].v.norm()*disk[j].v.norm() - MASS/disk[j].r.norm();
      }
      */

      save(disk[j], densitydat, t);
      densitydat.close();
      ++i;
    }
    ++j;
  }

  //std::cout << "End denscalc at time t = " << t << std::endl;

}   //end looping over particles






//OUTPUTS RATE OF COMPUTATION PLOT
void computation_time(std::string str, clock_t initial, std::ofstream & outf)
{
  clock_t current = clock();
  float seconds = (float(current-initial))/(CLOCKS_PER_SEC);
  std::cout << "Percent complete: "<< str <<", Minutes elapsed: "<< seconds/60 << "\r"<<std::flush;
  outf << str <<" "<<(float(current-initial))/(CLOCKS_PER_SEC)<<std::endl;
}











//Heating scheme:
const double X_s = 0.1;     //seed nuclei mass fraction
const double X_n = 1-X_s;   //neutron mass fraction

const double B_s = 8.7;     //in MeV, binding energy for seed nuclei
const double B_r = 8;       //in MeV, binding energy for r-process nuclei
const double diff = 1.293;  //in MeV, this is neutron proton rest mass difference
const double neutrino_f = .5;        //fraction of energy lost to neutrinos
const double E_r = (1-neutrino_f)*(B_r - X_s*B_s - X_n*diff);
const double total_E   = E_r * 1.602e-13;    //joules
const double proton_mass = 1.673e-27; //atomic mass unit to kg

//possibly variable in future. 2.03e5 code time in a second.
const double t_heat = 0.7/4.927e-6;  //total heating time in code time
//const double time_step = 2;
const double heat_rate = total_E/t_heat;    //in Joules, heat per nucleon over this time step

//convert back to cgs (c=1)
const double delta_u = heat_rate/proton_mass/(3e8*3e8);


void add_heat(Particle & p, double dt, double t)
{
  //QUESTION: how is 3MeV determined to be per nucleon? how does the equation work?

  if (t < t_heat)
  {
    p.v = (sqrt(p.v.norm()*p.v.norm() + 2*delta_u*dt) )/p.v.norm()*p.v;
  }

}





//MAIN EVOLUTION FUNCTION---------------------------------------------------------------
//--------------------------------------------------------------------------------------
//Evolves all particles, automatically stops each particle after reaching periapse, or at end of simulation
//also evolves density






void generate_smaller_disk(std::vector<Particle> & disk, double t, double dt, std::ofstream & tofile)
{
  std::vector<Particle> next_disk;
  int deleted=0;

  for (size_t n = 0; n < disk.size(); n++)
  {
    //enter if CURRENT radial velocity is positive
    if (check(disk[n].r, disk[n].v) == 1)
    {

      //enter this case only if velocity ever used to be negative
      if (disk[n].vcheck == 0)
      {
        
        if (disk[n].dflag==1)
        {
          Dcount = Dcount-1;
          std::cout << "Initially tracked "<< i_D_count <<" particles. Now down to "<< Dcount << std::endl; 
          disk[n].dflag=0;     
        }
        save(disk[n], tofile, t);                   //save particle data to file, along with time halted


        ++deleted;
        //std::cout << "p" << n+1 << " was deleted" << std::endl;
      }

      //enter if within range but moving outwards, keep updating
      else
      {
        //heat
        
        if (Heat ==1)
        {
          add_heat(disk[n], dt, t);
        }
        RKupdate(disk[n], dt);
        next_disk.push_back(disk[n]);
      }
    } 
    //ENTER IF particle is MOVING INWARDS
    else
    {
      disk[n].vcheck = 0;       //affirms that the particle _has_ in the past had negative radial velocity
      RKupdate(disk[n], dt);    //update for the case of falling particle
      next_disk.push_back(disk[n]);
    }
   //} //end the stopping condition loop


  }   //end looping over particles


  //CONSIDER SWITCHING TO LIST, AND STOP USING STD::VECTOR
  disk = next_disk;
  //std::cout << deleted<< " particles were deleted." << std::endl;
}














void delete_particles_same_disk(std::vector<Particle> & disk, double t, double dt, std::ofstream & tofile)
{
//begin looping over all particles
  for (size_t n = 0; n < disk.size(); n++)
  {
    int deleted=0;
    //enter if CURRENT radial velocity is positive
    if (check(disk[n].r, disk[n].v) == 1)
    {

      //enter this case only if velocity ever used to be negative
      if (disk[n].vcheck == 0)
      {
        
        if (disk[n].dflag==1)
        {
          Dcount = Dcount-1;
          disk[n].dflag=0;              
        }
        save(disk[n], tofile, t);                   //save particle data to file, along with time halted

          
        disk.erase(disk.begin() + n);   //when a particle is deleted
        n = n - 1;               //subtract 1 from the index, so the next particle isn't skipped

        ++deleted;
        //std::cout << "p" << n+1 << " was deleted" << std::endl;
      }

      //enter if within range but moving outwards, keep updating
      else
      {
        //heat
        if (Heat ==1)
        {
          add_heat(disk[n], dt, t);
        }
        RKupdate(disk[n], dt);

      }
    } 
    //ENTER IF particle is MOVING INWARDS
    else
    {
      disk[n].vcheck = 0;       //affirms that the particle _has_ in the past had negative radial velocity
      RKupdate(disk[n], dt);    //update for the case of falling particle

    }
   //} //end the stopping condition loop


  }//end loop over particles

}

















std::vector<Particle> autevolve(std::vector<Particle> & disk, const double & dt, std::ofstream & tofile, const double T = 2e5)
{ 

  //CANNOT CALCULATE DENSITY WHILE
  //LOOPING THROUGH ALL PARTICLES UPDATING POSITIONS, SINCE SOME PARTICLES WOULD BE UPDATED WHILE OTHERS WOULDN'T BE
  
  //SO MUST CALCULATE DENSITIES, THEN UPDATE POSITIONS
  //LOOP THROUGH DISK. IF DCHECK =1, THEN CALCULATE DENSITY.
  double t = 0;
  std::ofstream comp_time("comp_time.dat");

  std::cout << "Finding particles to track..." << std::endl;
  set_dflags(disk);

  std::cout << "Calculating initial densities..." << std::endl;
  density_calculation(disk, t);


  t = dt;
  int i = 1;
  std::cout << "Number of total steps ~ " << T/dt << std::endl;  
  std::cout << "Commencing evolution... " << std::endl;
  clock_t initial = clock();

  while (t < T)
  {


    //PROGRESS BAR  ///////////////////////////
    std::string percent_complete = std::to_string(t/T*100);
    computation_time(percent_complete, initial, comp_time);
    ///////////////////////////////////////////



    //UPDATE ENERGIES
    for (int i = 0; i < disk.size(); ++i)
    {
      disk[j].u = 1 + .5*disk[j].v.norm()*disk[j].v.norm() - MASS/disk[j].r.norm();
    }
    

    //faster if lots of particles being deleted
    if (t<1.6e4)
    {
      generate_smaller_disk(disk, t, dt, tofile);
    }
    //faster if not many particles being deleted
    else
    {
      delete_particles_same_disk(disk, t, dt, tofile);
    }





    // DENSITY CALCULATION before evolving particles ///////////
    int mod = dt*10;
    if (    (i % mod == 0) && (nnn<= int(disk.size()))    )
    {

      density_calculation(disk, t);
    }
    ////////////////////////////////////////////////////////////


    //CONSIDER SWITCHING TO LIST, AND STOP USING STD::VECTOR

    t = t + dt;
    i = i + 1;
    //std::cout << deleted<< " particles were deleted." << std::endl;



  }     //finish time whileloop

  clock_t final = clock();
  std::cout << std::endl;
  std::cout << "Complete. Minutes elapsed: "<<(final - initial)/(CLOCKS_PER_SEC)/60. << std::endl;

  return disk;

}
/////////////////////////////-----------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------











//This function TAKES OUTPUT LIST OF ONLY UNBOUND PARTICLES and calculates phi, cos(theta). Also prints masses and ye;
    //Creates output file to be used in solidangle.py
void solidangle(std::vector<Particle> & disk, std::ofstream & outf)
{

  for (size_t n = 0; n < disk.size(); n++)
  {
    outf << disk[n].m << " ";
    disk[n].print_position(outf);
    disk[n].print_velocity(outf);
    disk[n].print_quantities(outf);
    outf << std::endl;
  }
}

/////////////////////////////////////////////







//This will evolve all particles for a predetermined time, N*dt
void manevolve(std::vector<Particle> disk, double dt, double T = 2e5){

  double t = 0;
  int   i = 0;
  while (t < T)
  {    //begin time loop

    //std::cout <<"t"<< i << "\t";
    //disk[0].preport();              //report just the first particle at each time


//begin looping over all particles for given time
    for (size_t n = 0; n < disk.size(); n++)
    {
      RKupdate(disk[n], dt);    //update particle
    }

  t = t + dt;
  i = i + 1;
  }
}

//write code to track 

//Evolve should take in dt, disk vectorlist
//RKupdate should just take in one particle, dt



//////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////



int main() {


  double dt;


  std::string dt_;
  std::string MASS_;
  std::string particledat_;
  std::string domain_restricted_;
  std::string nnn_;
  std::string EVOLTIME_;
  std::string Heat_;

  std::string inputfile_name = "density.input";
  std::ifstream inputfile(inputfile_name);
  
  while (!inputfile.eof()) {         //while the end of file is not reached,
                                    //peek at first character.   
    if (inputfile.peek() == '#') inputfile.ignore(256, '\n');  //does the line start with "#"? If so, ignore until new line
      
    else { //reading bhsim.input and assigning each of its entries to appropriate objects
      std::getline(inputfile, dt_);
      dt = std::stod(dt_);

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
      Heat = std::stoi(Heat_);            


    }
  }

  std::cout << "Confirm the following:" << std::endl;
  std::cout << "Input file: "<< inputfile_name << std::endl;  
  std::cout << "Reading from: "<< particledat_<< std::endl;
  std::cout << "Time step: "<< dt<< std::endl;
  std::cout << "Mass of black hole: "<< MASS<< std::endl;
  std::cout << "Is domain restricted? 1 = yes, 0 = no: "<< domain_restricted<< std::endl;
  std::cout << "Number of nearest neighbors per tracker particle: "<< nnn<< std::endl;
  std::cout << "Total evolution time: "<< EVOLTIME << std::endl;


  if (Heat==1) std::cout << "Heating is on." << std::endl;
  else std::cout << "Heating is off." << std::endl;



  std::ifstream particledat(particledat_);
  std::ofstream fallbackdat("fallback.dat");
  std::ofstream solidangledat("solidangle.dat");

  //USER INPUT: If 0, then domain unchanged. If 1, it's restricted according to read(). 

  if (domain_restricted==1) std::cout << "Confirm that domain is restricted as desired." << std::endl;
  std::cout << "Confirm that file is in the format:" << std::endl;
  std::cout << "read_file >> m >> x >> y >> z >> vx >> vy >> vz >> u >> rho >> temp >> ye >> e" << std::endl;






  std::vector<Particle> disk;
  disk  = read(particledat);

///////////////////

  fallbackdat << "# m << x << y << z << vx << vy << vz << u << rho << temp << ye << e << t" << std::endl;
  solidangledat << "# m << x << y << z << vx << vy << vz << u << rho << temp << ye << e" << std::endl;

  

  std::vector<Particle> ejecta;
  ejecta = autevolve(disk, dt, fallbackdat, EVOLTIME);
  solidangle(ejecta, solidangledat);


  std::vector<Particle> densities;
  densities = ejecta;
//  densities = densevolve(particledat, densities, dt, EVOLTIME);
  //manevolve(disk, dt);
  //std::cout << newdisk[0].m << ", " << newdisk[0].r << ", " << newdisk[0].v << std::endl;


  return 0;
}





