#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include "Eigen_lib/Eigen/Dense"
#include <cstdio>
#include <time.h>
#include <sstream>
#include <dirent.h>

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

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//Function Declarations
std::vector<std::string> split(std::string str, char delimiter);
std::vector<std::string> generate_files(std::string path, std::string ending);
class Heat_profile;
class Particle;
void read_heating(std::ifstream & heatingdat, std::vector<double> * time, std::vector<double> * dt, std::vector<double> * temp, 
                  std::vector<double> * dens, std::vector<double> * entr, std::vector<double> * heat_rate, std::vector<double> * ye);
void set_heat_profile(std::string dir, std::vector<Heat_profile> * all_heat_profiles);
void assign_heat_profile_to_particle(std::vector<Particle> * particle_list, 
                                     std::vector<Heat_profile> * all_heat_profiles);
void add_skynet_heat(Particle *p, Heat_profile *h);











//~~~~~~~~~~~~~~~~~~~~~~~~~Function definitions~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


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
  int heat_prof = 0;
  int heat_prof_i = 0;
  

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




class Heat_profile
{      //this is a vector of particles
public:
  std::vector<double>  time;
  std::vector<double>  dt;
  std::vector<double>  temp;
  std::vector<double>  dens;
  std::vector<double>  entr;
  std::vector<double>  heat_rate;
  std::vector<double>  ye;

  

  Heat_profile (std::vector<double>  time_, std::vector<double>  dt_, std::vector<double>  temp_, std::vector<double>  dens_,
                std::vector<double>  entr_, std::vector<double>  heat_rate_, std::vector<double>  ye_) :    //constructor for object of Particle class
               time(time_),  dt(dt_),  temp(temp_),  dens(dens_),  entr(entr_),   heat_rate(heat_rate_),   ye(ye_) {}


/*
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
  }*/

};







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





void read_heating(std::ifstream & heatingdat, std::vector<double> * time, std::vector<double> * dt, std::vector<double> * temp, 
                  std::vector<double> * dens, std::vector<double> * entr, std::vector<double> * heat_rate, std::vector<double> * ye){
  int j = 0;
  std::string data;
  std::string line;
  j=0;
  while (std::getline(heatingdat, data)) {
    std::vector<std::string> line = split(data, ' ');
    if (line.size()==0) break;
    char tab2[1024];
    strncpy(tab2, line[0].c_str(), sizeof(tab2));
    tab2[sizeof(tab2) - 1] = 0;
    //printf("size: %d\n", strncmp(tab2,"#",1) );
    if (strncmp(tab2,"#",1)==0 || line.size() == 0); 
    else
    {
      //printf("test: %s, %s, %s, %s\n", line[0].c_str(),line[1].c_str(),line[2].c_str(),line[3].c_str());
      time->push_back( std::stod(line[0])  );
      (*dt).push_back  ( std::stod(line[1])  );
      (*temp).push_back( std::stod(line[2])  );
      (*dens).push_back( std::stod(line[3])  );
      (*entr).push_back( std::stod(line[4])  );
      (*heat_rate).push_back( std::stod(line[5])  );
      (*ye).push_back(   std::stod(line[6])  );
      //printf("time[%i]: %e\n", j, (*time)[j]);
      ++j;

    }

  }

}


















std::string convert_next_line_to_string(std::ifstream & inputfile)
{
  std::string line;
  std::getline(inputfile, line);
  return line;
}

std::vector<std::string> split(std::string str, char delimiter)
{
  std::vector<std::string> internal;
  std::stringstream ss(str); // Turn the string into a stream.
  std::string tok;

  
  while(getline(ss, tok, delimiter)) {

    if (tok != "" ) {
      internal.push_back(tok);
    }
  }
  
  return internal;
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
  double m, x, y, z, vx, vy, vz, u, rho;

  Eigen::Vector3d r;
  Eigen::Vector3d v;
  std::vector<Particle> tail;
  Particle p(m, r, v);

  int i = 0;

  //std::ofstream corrected("corrected.dat");
  std::string data;
  std::string line;
  while (std::getline(particledat, data))
  {
        //peek at first character.   
    //if (particledat.peek() == '#') particledat.ignore(256, '\n');  //does the line start with "#"? If so, ignore until new line
    std::vector<std::string> line = split(data, ' ');
    if (line[0] == "#");  //does the line start with "#"? If so, do nothing in this iteration of loop.

    else
    {

      m   = std::stod(line[0]);
      x   = std::stod(line[1]);
      y   = std::stod(line[2]);
      z   = std::stod(line[3]);
      vx  = std::stod(line[4]);
      vy  = std::stod(line[5]);
      vz  = std::stod(line[6]);
      u   = std::stod(line[7]);
      rho = std::stod(line[8]);
      const double temp = std::stod(line[9]);
      const double ye   = std::stod(line[10]);
      const double e    = std::stod(line[11]);






      //particledat >> m >> x >> y >> z >> vx >> vy >> vz >> u >> rho >> temp >> ye >> e;
      
      //RESTRICTION ON DOMAIN
      //If domain_restricted=0, enter block.
      //If domain is restricted, but particle lies in subdomain, enter block. Else do nothing.
      double radius = sqrt(x*x + y*y + z*z);
      double vnorm = sqrt(vx*vx + vy*vy + vz*vz);
      double phi = atan2(y,x);

      //example of restriction
      //(pi/-2 < phi && phi < pi/3 && radius < 40 ) || ( phi > pi/3 || phi < pi/-2 )  





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



//std::vector<Heat_profile> all_heat_profiles;
//all_heat_profiles[1].time[0:-1]


std::vector<std::string> generate_files(std::string path, std::string ending)       //opening any folder and saving all file-names in a vector<string>
{
    DIR    *dir;
    dirent *pdir;
    std::string entry;
    std::vector<std::string> files;
    printf("%s\n", path.c_str());

    dir = opendir(path.c_str());
    std::cout << typeid(path.c_str()).name() << std::endl;
    while ((pdir = readdir(dir)) !=NULL   )
    {
      entry = pdir->d_name;
      //printf("first %s, %lu\n", entry.c_str(), entry.length());
      int min_char = ending.length();
      if (int(entry.length())>min_char && entry.substr( int(entry.length()) - min_char ) ==ending)
      {
        files.push_back(entry);
        //printf("pdir->d_name: %s\n", files[files.size()-1].c_str());

      }
    }

    return files;
}





void set_heat_profile(std::string dir, std::vector<Heat_profile> * all_heat_profiles)
{//create vector of heating profiles based on all relevant files in dir

  //move to dir;
  std::vector<std::string> files = generate_files(dir, "h5.dat");       //opening any folder and saving all file-names in a vector<string>
  std::string file_name;
  printf("number of files: %i\n", int(files.size()));
  //begin loop over files
    for (int i = 0; i < int(files.size()); ++i){ 
      file_name = dir+"/"+files[i];
      std::ifstream heatingdat(file_name);
      std::vector<double>  time; //seconds
      std::vector<double>  dt;    //seconds
      std::vector<double>  temp;
      std::vector<double>  dens;
      std::vector<double>  entr;
      std::vector<double>  heat_rate; //ergs/s/g
      std::vector<double>  ye;
      read_heating(heatingdat, &time, &dt, &temp, &dens, &entr, &heat_rate, &ye);
      
      //TO DO: remove initial vector elements for temp > 3e9 K
      bool keepgoing=1;
      int j =0;
      while(keepgoing)
      {
        if (temp[j] < 3.0)
        {
          keepgoing=0;
        } else {
          time.erase(time.begin());
          dt.erase(dt.begin());
          temp.erase(temp.begin());
          dens.erase(dens.begin());
          entr.erase(entr.begin());
          heat_rate.erase(heat_rate.begin());
          ye.erase(ye.begin());
        }
        ++j;
      }


      Heat_profile h(time, dt, temp, dens, entr, heat_rate, ye);
  
      double dh=0;
      double total_heat=0;
      for (int i = 0; i < int(dt.size()); ++i)
      {
        dh=dt[i]*heat_rate[i];
        total_heat += dh;
      }

      const double proton_mass = 1.673e-24; //atomic mass unit to grams
      double erg_per_g_2_MeV_per_nuc= 6.2e5*proton_mass;
      total_heat *= erg_per_g_2_MeV_per_nuc;
      printf("total_heat = %g\n", total_heat);
      all_heat_profiles->push_back(h);
    }
  //end loop over files


    //exit(1);



}


void assign_heat_profile_to_particle(std::vector<Particle> * particle_list, 
                                     std::vector<Heat_profile> * all_heat_profiles)
{//in beginning of simulation, assign each particle a heating profile based solely on Y_e

  printf("size of all_heat_profiles: %i\n",int(all_heat_profiles->size()) );
  for (int i = 0; i < int(particle_list->size()); ++i)
  {
    Particle p = (*particle_list)[i];
    double dist = 1;
    double current_dist = dist;
    for (int j = 0; j < int(all_heat_profiles->size()); ++j)
    {
      current_dist = fabs(p.ye-(*all_heat_profiles)[j].ye[0] );
      if (current_dist<dist)
      {
        p.heat_prof=j;
        dist = current_dist;
      }

    }


    printf("p.id = %i, p.heat_prof= %i, p.ye =%g, all_profs.ye =%g\n",p.id, p.heat_prof, p.ye, (*all_heat_profiles)[p.heat_prof].ye[0]);


    double min_ye=1;
    double max_ye =0;
    int min_index = 0;
    int max_index = 0;
    for (int i = 0; i < int(all_heat_profiles->size()); ++i)
    {
      if (min_ye>(*all_heat_profiles)[i].ye[0])
      {
        min_ye=(*all_heat_profiles)[i].ye[0];
        min_index = i;
      } if (max_ye<(*all_heat_profiles)[i].ye[0])
      {
        max_ye=(*all_heat_profiles)[i].ye[0];
        max_index = i;

      }
    }
    //printf("min:(all_profs.ye[%i])= %g, max: (all_profs.ye[%i])= %g\n", min_index, min_ye, max_index, max_ye);
    (*particle_list)[i] = p;
  }
  //loop over particles



}



void add_skynet_heat(Particle *p, Heat_profile *h){
  //at each time step, add heat to a particle p from its heat profile
  const double proton_mass = 1.673e-24; //atomic mass unit to grams
  const double one_sec = 2.27e5;  //1s = 2.27e5 codetime units

  //skynet time is in seconds
  double t_tot = T_Heat*one_sec;  //total heating time in code time
  //const double time_step = 2;
  double erg_per_g_2_MeV_per_nuc= 6.2e5*proton_mass;
  double MeV_2_J = 1.602e-13;


  double shift_avg_heating_down = 0.5;

  //do time unit conversion here:
  double t_sec = p->t/one_sec; //in seconds now


  bool keepgoing=1;
  int i = p->heat_prof_i;
  if (i>=int((h->time).size())-1)
  {
    keepgoing=0;
  }
  while(keepgoing){


    if (h->time[i] <= t_sec && h->time[i+1] >= t_sec) {

      double del_t = h->time[i+1]-h->time[i];
      double del_h = h->heat_rate[i+1] - h->heat_rate[i];



      double H_weighted = (t_sec - h->time[i])*(del_h/del_t) + h->heat_rate[i]; 
      double heat_rate = shift_avg_heating_down*H_weighted/one_sec* erg_per_g_2_MeV_per_nuc*MeV_2_J;//joules per nucleon per code unit of time
      //double heat_rate = total_E/t_tot;

      //convert back to cgs (c=1)
      //const double total_E   = E_r * 1.602e-13;
      double f_neutrino = 0.75;    
      double delta_u = f_neutrino*heat_rate/proton_mass/1e-3/(3e8*3e8);


      p->v = (sqrt(p->v.norm()*p->v.norm() + 2*delta_u*p->dt) )
              /p->v.norm()*p->v;
      //printf("heat added for particle %i, profile %i\n",p->id, p->heat_prof );

      p->heat_prof_i = i;
      keepgoing=0;
    } else if (i > int((h->time).size())-1) {


      keepgoing = 0;
    }
    ++i;

/*      double dh=0;
      double total_heat=0;
      for (int j = 0; j < int(dt.size()); ++j)
      {
        dh=dt[j]*heat_rate[j];
        total_heat += dh;
      }

      const double proton_mass = 1.673e-24; //atomic mass unit to grams
      double erg_per_g_2_MeV_per_nuc= 6.2e5*proton_mass;
      total_heat *= erg_per_g_2_MeV_per_nuc;
      printf("total_heat = %g\n", total_heat);*/

  }

}



/////////////////////////////////////////////////////////////////////////////////////////////////





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
double X_s = 0.1;     //seed nuclei mass fraction
double X_n = 1-X_s;   //neutron mass fraction

const double YE_avg = 0.0545;

const double Z = 36;
const double A = 118;
const double B_s = 8.7;     //in MeV, binding energy for seed nuclei
const double B_r = 8;       //in MeV, binding energy for r-process nuclei
const double dn = 1.293;  //in MeV, this is neutron proton rest mass difference
const double neutrino_f = .5;        //fraction of energy lost to neutrinos

const double proton_mass = 1.673e-24; //atomic mass unit to grams

//possibly variable in future. 2.03e5 code time in a second.
const double erg_in_joules = 1e-7;   //1erg = 1e-7 joule



void add_heat(Particle & p)//, Heat_profile * h)
{

/*
  X_s = A*p.ye/Z;
  X_n = 1-X_s;

  double E_r = Heat;
*/

  double A_Z = (2*Heat-B_r+dn )/(YE_avg*(dn-B_s));
  double X_s = A_Z*p.ye;
  double X_n = 1-X_s;

  double E_r = (1-neutrino_f)*(B_r - X_s*B_s - X_n*dn);

  const double total_E   = E_r * 1.602e-13;    //joules per nucleon
  const double one_sec = 2.27e5;  //1s = 2.27e5 codetime units
  double t_tot = T_Heat*one_sec;  //total heating time in code time
  //const double time_step = 2;
  double heat_rate = total_E/t_tot;

  //convert back to cgs (c=1)
//const double total_E   = E_r * 1.602e-13;    //joules per nucleon
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







bool generate_tail(std::vector<Particle> & tail, double tol, double dcalc, std::ofstream & tofile, std::vector<Heat_profile> * all_heat_profiles)
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
          add_skynet_heat(&tail[n], (&(*all_heat_profiles)[(tail[n]).heat_prof]));



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















bool reduce_tail(std::vector<Particle> & tail, double tol, double dcalc, std::ofstream & tofile, std::vector<Heat_profile> * all_heat_profiles)
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
          add_skynet_heat(&tail[n], (&(*all_heat_profiles)[(tail[n]).heat_prof]));
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


  //choosing times at which we compute density
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
  while (dtime[h_index]<T_Heat*2.27e5) {
    ++h_index;
  }
  dtime.insert(dtime.begin()+h_index, T_Heat*2.27e5);



  std::ofstream comp_time("comp_time.dat");

  //flagging particles which will be tracked for density
  std::cout << "Finding particles to track..." << std::endl;
  set_dflags(tail);

  std::cout << "Calculating initial densities..." << std::endl;
  density_calculation(tail);


  std::cout << "Number of total density updates: " << dtime.size() << std::endl;  
  std::cout << "Commencing evolution... " << std::endl;
  clock_t initial = clock();


  std::vector<Heat_profile> all_heat_profiles;
  std::string dir = "/Users/dhruv/Documents/Research/Visuals/Skynet_Run1";
  set_heat_profile(dir, &all_heat_profiles);
  assign_heat_profile_to_particle(&tail, &all_heat_profiles);



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
    

    bool advance;
    //faster if lots of particles being deleted
    if (t<1.6e4)
    {
      advance = generate_tail(tail, tol, t, tofile, &all_heat_profiles);
    }
    //faster if not many particles being deleted
    else
    {
      advance = reduce_tail(tail, tol, t, tofile, &all_heat_profiles);
    }

    //bool advance = generate_tail(tail, tol, t, tofile);    
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
  std::cout << "Is domain restricted? 1 = yes, 0 = no: "<< domain_restricted << std::endl;
  std::cout << "Number of nearest neighbors per tracker particle: "<< nnn<< std::endl;
  std::cout << "Total evolution time: "<< EVOLTIME << std::endl;


  std::cout << "Heating: " << Heat << " MeV/nucleon."<< std::endl;
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

  double avgye=0;
  double mtot = 0;
  for (int i = 0; i < int(tail.size()); ++i)
  {
    avgye += tail[i].ye*tail[i].m;
    mtot += tail[i].m;
  }
  avgye = avgye/mtot;
  printf("mass avg of ye: %e\n", avgye);

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





