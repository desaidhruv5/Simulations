#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include "Eigen_lib/Eigen/Dense"

//////////////////////////////////////////////////////////////////*

/*
That's about 1.5 million points, because it includes every single grid point for which
1) The density is greater than 1.e-10 
2) The radius is larger than 40
A note on units: we set G=c=1 and solar mass = 1, hence
Density: 1 = 6.1e17 g/cc
Distance: 1 = 1.48km
Velocity: 1 = 3e10 cm/s
The first few lines give you the list of variables. 
"Minus_ut" is the specific energy (in case we want to use that instead of velocities). A priori, all particles with Minus_ut>1 should be unbound.
Ye is the electron fraction. 
*/



//MASS OF BLACK HOLE
double MASS;

//If set to 1, domain will be restricted according to formula in read() function.
//If set to 0, bhsim will incorporate all particles.
int domain_restricted;


//checks direction of particle velocity
//returns 1 if radial velocity is positive, else returns 0
bool check (Eigen::Vector3d r, Eigen::Vector3d v) {


  bool b;
  if (r.dot(v) > 0) b = 1;
  else b = 0;

  return b;
}

///////////////



class Particle {    //this is a vector of particles
public:
  double m;
  Eigen::Vector3d r;
  Eigen::Vector3d v;
  double u;
  double ye;
  double e;
  double rho;
  bool vcheck;


  Particle (double m_, Eigen::Vector3d r_, Eigen::Vector3d v_, double u_, double ye_) :    //constructor for object of Particle class
               m(m_),               r(r_),              v(v_),     u(u_),    ye(ye_) {}



  void
  preport() {                              //report particle info
    std::cout << "Mass:\t"<< m
    << " x = " << r[0]
    << " y = " << r[1]
    << "\tz = " << r[2] << std::endl;
  }
  
  void
  vreport() {                              //report particle info
    std::cout << "Mass:\t"<< m
    << " vx = " << v[0]
    << " vy = " << v[1]
    << "\tvz = "<< v[2] << std::endl;
  }

};
///////////////////////////////////////////



//later this will be array
//DEFINING function 'Read' takes in file, outputs vector of particles, a.k.a. 'disk', also corrects velocities based on energy
std::vector <Particle> read(std::ifstream & particledat) {

//STEPS:
//reads each line, 1: creates object of type particle, 2:labels it, 3: assigns position, 4: assigns velocity, 5: ...

  //here we define what each element on each line represents
  double m, x, y, z, vx, vy, vz, u, ye, e, rho;

  Eigen::Vector3d r;
  Eigen::Vector3d v;
  std::vector<Particle> disk;
  Particle p(m, r, v, u, ye);

  int i = 0;


  int illegalvcount = 0;
  while (!particledat.eof()) {

    


        //peek at first character.
    if (particledat.peek() == '#') particledat.ignore(256, '\n');  //does the line start with "#"? If so, ignore until new line
      
    else {
      particledat >> m >> x >> y >> z >> vx >> vy >> vz >> u >> ye >> e >> rho;
      //RESTRICTION ON DOMAIN
      //If domain_restricted=0, enter block.
      //If domain is restricted, but particle lies in subdomain, enter block. Else do nothing.

      double r = sqrt(x*x + y*y + z*z);
      double phi = atan2(y,x);
      double pi = 4*atan(1);


      if(
      ( domain_restricted != 1)    ||/*Inside the following goes dom. restriction */
      (  (pi/-2 < phi && phi < pi/3 && r < 40 ) || ( phi > pi/3 || phi < pi/-2 )      )
      )
      
      {
        disk.push_back(p);              //creating entry for ith particle in the disk
        disk[i].m    = m;               //mass of ith particle in the disk
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
        disk[i].vcheck = check (disk[i].r, disk[i].v); //checking if intial radial velocity is positive/negative



        if (std::abs(disk[i].v.norm()) > 1)   ++illegalvcount;



        //now correct the velocities based on energy ( E = KE + PE )
        
        double correctv = (   disk[i].u - 1 + MASS/disk[i].r.norm()   )*2;

        correctv = sqrt(std::abs(correctv));
        disk[i].v = disk[i].v * correctv/disk[i].v.norm();


        //-pi/4 to pi/4, r > 70

        //std::cout << "correct: " << correctv << std::endl;
        //std::cout << "actual: " << disk[i].v.norm() << std::endl;


        //std::cout << "i=" << i << std::endl;
        //std::cout << "m=" << disk[i].m << std::endl;
        ++i;
      }
      
    }
    
  }
  std::cout << "# of particles with |v| > 1 BEFORE v correction: "<< illegalvcount<< std::endl;

  return disk;
}


//COMPUTING FIELD at position r away from black hole, which is at origin
Eigen::Vector3d fieldon(Particle p) {
  const double G = 1;
  //const double M = 7.942;    //mass of black hole, only source of gravitational field

  double R = p.r.norm();

  //know mass, position
  return p.r*(-1*G*MASS/(R*R*R));
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
Particle RKupdate(Particle p1, double dt)
{ //This updates 1 particle via RK 4th order method
                                              //after a time step of 'dt'

  Eigen::Vector3d dv;
  Eigen::Vector3d dr;

  Particle p = p1;
  //std::vector<double> j1 = p.v;                 //velocity at initial
  //Field k1    = fieldon(p);          //field at initial

  p.r = p1.r + .5*dt*p.v;                       //approximates field at the 
  Eigen::Vector3d j2 = p.v + dt*.5*fieldon(p1);    //velocity at midpt using field initial
  Eigen::Vector3d k2 = fieldon(p);     //field at midpt using velocity initial

  p.r = p1.r + .5*dt*j2;
  Eigen::Vector3d j3 = p.v + dt*.5*k2;    //velocity at midpt. using field 
  Eigen::Vector3d k3 = fieldon(p);

  p.r = p1.r + dt*j3;
  Eigen::Vector3d j4 = p.v + dt*k3;
  //std::vector<double> k4    = fieldon(p);

  dv = dt/6.* ( fieldon(p1) + 2*k2  + 2*k3  +  fieldon(p)  );       //find dv
  dr = dt/6.* ( p.v         + 2*j2 +  2*j3  +  j4   );       //find dr

  //update this info back into particle info:
  p.r = p1.r + dr; 
  p.v = p1.v + dv;

  return p;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

//adds particle data to output file (in cartesian)
void save(Particle p, std::ofstream & outf, double t) {

  outf << p.m << " " << p.r[0] << " " <<
  p.r[1] << " " << p.r[2] << " " <<
  p.v[0] << " " << p.v[1] << " " <<
  p.v[2] << " " << p.u << " " << p.ye << " " << p.e << " " << p.rho << " " << t << std::endl;


}

//This saves particle data in cylindrical coordinates, projecting everything on z vs. r graph
//(getting rid of phi dependence)

/*
void savecyl(Particle p, std::ofstream & outf, double t) {

  double r = sqrt(p.r[1]*p.r[1] + p.r[0]*p.r[0]);

  outf << p.m << " " << r << " " << p.r[2] << " " <<
  p.v[0] << " " << p.v[1] << " " <<
  p.v[2] << " " << p.u << " " << p.ye << " " << p.e << " " << p.rho << " " << t << std::endl;

}
*/


//This is main evolution function which automatically stops evolution of each particle appropriately
std::vector<Particle> autevolve(std::vector<Particle> disk, double dt, std::ofstream & tofile, double T = 2e5) {

  tofile << "# 1 m    2 rx  3 ry   4 rz   5 vx     6  vy  7  vz     8 u       9 ye       10 e      11 rho       12 time" << std::endl;

  double t = 0;
  int i = 0;

  std::cout << "Number of total steps ~ " << T/dt << std::endl;
  //BEGIN time loop
  while (t < T) {
    
    //std::cout <<"t"<< i << "\t";
    //disk[0].preport();           //report position of just the first particle



    if (i % 1000 == 0)
    {
      std::cout << "step " << i << std::endl;
      std::cout << t/T*100 << " percent complete." << std::endl;
    }

      //then record the time

    //begin looping over all particles

    for (size_t n = 0; n < disk.size(); n++)
    {

        //try checking only every 10 steps, so that I can update accurately, but not enter the for loop as many times
        //this won't give me as accurate of an answer for the turning point, but it will be faster, as for loops are slow
     //if (i % 1 == 0) {

      //enter if CURRENT radial velocity is positive
      if (check(disk[n].r, disk[n].v) == 1)
      {

        //enter this case only if velocity ever used to be negative
        if (disk[n].vcheck == 0)
        {
          save(disk[n], tofile, t);                   //save particle data to file, along with time halted
          disk.erase(disk.begin() + n);   //when a particle is deleted
          n = n - 1;               //subtract 1 from the index, so the next particle isn't skipped
          //std::cout << "p" << n+1 << " was deleted" << std::endl;
        }

        //enter if within range but moving outwards, keep updating
        else
        {
          disk[n] = RKupdate(disk[n], dt);
          //std::cout << "p" << n << " was updated, still moving outwards" << std::endl;
        }

      }


      //ENTER IF particle is MOVING INWARDS
      else {
        disk[n].vcheck = 0;       //affirms that the particle _has_ in the past had negative radial velocity
        disk[n] = RKupdate(disk[n], dt);    //update for the case of falling particle
        //std::cout << "p" << n << " is moving inwards"<< std::endl;
      }

     //} //end the stopping condition loop

     //else {
     //  disk[n] = RKupdate(disk[n], dt); 
     //}

    }   //end looping over particles

    t = t + dt;
    i = i + 1;
    //std::cout << i << std::endl;

  }     //finish time whileloop



  //AFTER TIME T, we must record the rest of the particles
  for (size_t n = 0; n < disk.size(); n++)
  {
    save(disk[n], tofile, t);
    //disk.erase(disk.begin() + n);   //when a particle is deleted
    //n = n - 1;               //subtract 1 from the index, so the next particle isn't skipped
  }

  return disk;

}
/////////////////////////////


void solidangle(std::vector<Particle> disk, std::ofstream & outf) {

  outf << "# phi    costheta    mass     ye   rx    ry     rz   vx      vy      vz      u       e       rho" << std::endl;
  
  for (size_t n = 0; n < disk.size(); n++)
  {
    if (disk[n].vcheck == 1)
    {
      double atanphi;
      double tan1;
      double tan0;
      double costheta;
      tan1 = disk[n].v[1];
      tan0 = disk[n].v[0];
      atanphi = atan2(tan1, tan0);
      costheta = disk[n].v[2]/disk[n].v.norm();

      outf << atanphi << " " << costheta << " " << disk[n].m << " " << disk[n].ye << " "
      <<disk[n].r[0]<< " " << disk[n].r[1] << " " << disk[n].r[2]<<" " << disk[n].v[0]
      <<" " << disk[n].v[1]<<" " << disk[n].v[2]<<" " << disk[n].u << " " << disk[n].e << " " << disk[n].rho << std::endl;
    }
  }

}

/////////////////////////////////////////////



void manevolve(std::vector<Particle> disk, std::ofstream & tofile, double dt, double T = 2e5){
//This will evolve all particles for a predetermined time, N*dt
  double t = 0;
  int    i = 0;
  std::cout <<"total time: "<< T << "\t"<< std::endl;
  
  //let's check to see how many have v >1, BEFORE updating position//////////////////
  int illegalvcount = 0;
  for (size_t n = 0; n < disk.size(); n++)
  {
    if (std::abs(disk[n].v.norm()) > 1)
    {
      ++illegalvcount;
    }
  }
  std::cout << "# of particles with |v| > 1 AFTER v correction, BEFORE update: "<< illegalvcount<< std::endl;
  /////////////////////////////////////////////////////////


  //BEGIN LOOPING OVER TIME
  while (t < T)
  {

    //std::cout <<"t"<< i << "\t";
    //disk[0].preport();              //report just the first particle at each time

    //begin looping over all particles for given time
    for (size_t n = 0; n < disk.size(); n++)
    {
      disk[n] = RKupdate(disk[n], dt);    //update particle
    }


    std::cout <<"fraction complete: "<< t/T << "\t"<< std::endl;
    t = t + dt;
    i = i + 1;
  }



 //Now let's check to see how many have v >1, AFTER updating position//////////////////
  illegalvcount = 0;
  for (size_t n = 0; n < disk.size(); n++)
  {
    if (std::abs(disk[n].v.norm()) > 1)
    {
      ++illegalvcount;
    }
  }
  std::cout << "# of particles with |v| > 1 AFTER update: "<< illegalvcount<< std::endl;
  /////////////////////////////////////////////////////////



  for (size_t n = 0; n < disk.size(); n++)
    {
        save(disk[n], tofile, t);    //save particle
    }

}

//RODRIGO's TIME: 3231.3 = ~16.03ms
//M14_M7_S8
//we want avg density over r, z
// plotting z/r as a function of log(r)

// plot the following:
// mass weighted average of vphi, vr , vz, rho, ye, entropy

//output into file, the radius, z value, and the averages stated above



//////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////



int main() {

  std::string dt_;
  std::string MASS_;
  std::string particledat_;
  std::string domain_restricted_;
  double dt;
  std::ifstream inputfile("bhsim.input");
  
  while (!inputfile.eof()) {         //while the end of file is not reached,
                                    //peek at first character.   
    if (inputfile.peek() == '#') inputfile.ignore(256, '\n');  //does the line start with "#"? If so, ignore until new line
      
    else
    { //reading bhsim.input and assigning each of its entries to appropriate objects
      std::getline(inputfile, dt_);
      dt = std::stod(dt_);
      std::getline(inputfile, MASS_);
      MASS = std::stod(MASS_);
      std::getline(inputfile, particledat_);
      std::getline(inputfile, domain_restricted_);
      domain_restricted = std::stoi(domain_restricted_);
    }
  }

  std::ifstream particledat(particledat_);
  std::ofstream fallbackdat("fallback.dat");
  std::ofstream solidangledat("solidangle.dat");
  std::ofstream finaldat("final.dat");

  //USER INPUT: If 0, then domain unchanged. If 1, it's restricted according to read(). 

  if (domain_restricted==1) std::cout << "Confirm that domain is restricted as desired." << std::endl;
  //

  std::vector<Particle> disk;
  disk  = read(particledat);

///////////////////

  //std::cout << disk[0].v << ", " << disk[0].r << "velocity" <<
  //disk[0].vcheck << std::endl;
  //std::cout << "dot"<< disk[0].r.dot(disk[0].v)<< std::endl;
  


  //if (i==58455) std::cout << "INITIAL VELOCITY: "<<disk[i].v <<std::endl;

  std::vector<Particle> newdisk;
  newdisk = autevolve(disk, dt, fallbackdat);
  //manevolve(disk, finaldat, dt, 3231.3);
  //std::cout << newdisk[0].m << ", " << newdisk[0].r << ", " << newdisk[0].v << std::endl;
  solidangle(newdisk, solidangledat);

  return 0;
}


//END OF PROGRAM            /////////////////////////////////////////

//TO DO NEXT:

//////////////////////
