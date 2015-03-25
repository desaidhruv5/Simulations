#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include "Eigen/Eigen/Dense"

//////////////////////////////////////////////////////////////////*


/*
That's about 1.5 million points, because I included every single grid point for which
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





//checks direction of particle velocity
//returns 1 if radial velocity is positive, else returns 0
bool check (Eigen::Vector3f r, Eigen::Vector3f v) {


  bool b;
  if (r.dot(v) > 0) b =1;
  else b = 0;

  return b;
}

///////////////



class Particle {    //this is a vector of particles
public:
  float m;
  Eigen::Vector3f r;
  Eigen::Vector3f v;
  float u;
  float ye;
  bool vcheck;



  Particle (float m_, Eigen::Vector3f r_, Eigen::Vector3f v_, float u_, float ye_) :    //constructor for object of Particle class
            m(m_),    r(r_),                          v(v_),     u(u_),     ye(ye_){}



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
//DEFINING function 'Read' takes in file, outputs a vector of particles, a.k.a. 'disk'
std::vector <Particle> read(std::ifstream & infile) {

//STEPS:
//reads each line, 1: creates object of type particle, 2:labels it, 3: assigns position, 4: assigns velocity, 5: ...

  //here we define what each element on each line represents
  float m, x, y, z, vx, vy, vz, u, ye, e, rho;

  Eigen::Vector3f r;
  Eigen::Vector3f v;
  std::vector<Particle> disk;
  Particle p(m, r, v, u, ye);

  int i = 0;
 

  while (!infile.eof()) {
        //peek at first character.   
    if (infile.peek() == '#') infile.ignore(256, '\n');  //does the line start with "#"? If so, ignore until new line
      
    else {
      infile >> m >> x >> y >> z >> vx >> vy >> vz >> u >> ye >> e >> rho;
      disk.push_back(p);              //creating entry for ith particle in the disk
      disk[i].m = m;                  //mass of ith particle in the disk
      disk[i].r[0] = x;          //ith particle of the disk, position vector, create x component
      disk[i].r[1] = y;          //...
      disk[i].r[2] = z; 
      disk[i].v[0] = vx;         //ith particle of the disk, position vector, create x component
      disk[i].v[1] = vy;         //...
      disk[i].v[2] = vz;
      disk[i].u    = u;
      disk[i].ye = ye;                //disk[i] is ith particle, ye is composition of particle
      disk[i].vcheck = check (disk[i].r, disk[i].v); //checking if intial radial velocity is positive/negative



      //now correct the velocities based on energy ( E = KE + PE )
      float correctv = (   disk[i].u - 1 + 7.942/disk[i].r.norm()   )*2;
      correctv = sqrt(correctv);
      std::cout << "correct: " << correctv << std::endl;
      std::cout << "actual: " << disk[i].v.norm() << std::endl;
      std::cout << (correctv - disk[i].v.norm() )/correctv << std::endl; 




      ++i;
    }
    
  }


  return disk;
}




//COMPUTING FIELD at position r away from black hole, which is at origin
Eigen::Vector3f fieldon(Particle p) {
  const float G = 1;
  const float M = 7.942;    //mass of black hole, only source of gravitational field

  float R = p.r.norm();

  //know mass, position
  return p.r*(-1*G*M/(R*R*R));
}




Particle RKupdate(Particle p1, float dt) { //This updates 1 particle via RK 4th order method
                                              //after a time step of 'dt'

  Eigen::Vector3f dv;
  Eigen::Vector3f dr;

  Particle p  = p1;
  //std::vector<float> j1 = p.v;                 //velocity at initial
  //Field k1    = fieldon(p);          //field at initial

  p.r = p1.r + .5*dt*p.v;                       //approximates field at the 
  Eigen::Vector3f j2 = p.v + dt*.5*fieldon(p1);    //velocity at midpt using field initial
  Eigen::Vector3f k2 = fieldon(p);     //field at midpt using velocity initial

  p.r = p1.r + .5*dt*j2;
  Eigen::Vector3f j3 = p.v + dt*.5*k2;    //velocity at midpt. using field 
  Eigen::Vector3f k3 = fieldon(p);

  p.r = p1.r + dt*j3;
  Eigen::Vector3f j4 = p.v + dt*k3;
  //std::vector<float> k4    = fieldon(p);


  dv = dt/6.* ( fieldon(p1) + 2*k2  + 2*k3  +  fieldon(p)  );       //find dv
  dr = dt/6.* ( p.v         + 2*j2 +  2*j3  +  j4   );       //find dr


  //update this info back into particle info:
  p.r = p1.r + dr; 
  p.v = p1.v + dv;


  return p;
}


void save(Particle p, std::ofstream & outf, float t) {

  outf << p.m << " " << p.r[0] << " " <<
  p.r[1] << " " << p.r[2] << " " <<
  p.v[0] << " " << p.v[1] << " " <<
  p.v[2] << " " << p.ye << t << std::endl;
}


void autevolve(std::vector<Particle> disk, float dt, std::ofstream & tofile, float T = 2e5) { //This will auto-stop
 
  float t = 0;
  int i = 0;

  std::cout << "Number of total steps ~ " << T/dt << std::endl;
  //BEGIN time loop
  while (t < T) {
    
    //std::cout <<"t"<< i << "\t";
    //disk[0].preport();           //report position of just the first particle



    if (i % 10000) {
      std::cout << "step " << i << std::endl;
      std::cout << t/T*100 << " percent complete." << std::endl;

    }

      //then record the time

    //begin looping over all particles

    for (size_t n = 0; n < disk.size(); n++) {

        //try checking only every 10 steps, so that I can update accurately, but not enter the for loop as many times
        //this won't give me as accurate of an answer for the turning point, but it will be faster, as for loops are slow
      if (i % 1 == 0) {
        //enter if CURRENT radial velocity is positive
        if (check(disk[n].r, disk[n].v) == 1) {

          //enter this case only if velocity ever used to be negative
          if (disk[n].vcheck == 0) {
            save(disk[n], tofile, t);                   //save particle data to file
            tofile << "t" << t << std::endl;
            disk.erase(disk.begin() + n);   //when a particle is deleted
            n = n - 1;               //subtract 1 from the index, so the next particle isn't skipped
            //std::cout << "p" << n+1 << " was deleted" << std::endl;
          }

          //enter if within range but moving outwards, keep updating
          else {
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

      } //end the stopping condition loop



      else {
        disk[n] = RKupdate(disk[n], dt); 
      }

    }   //end looping over particles

    t = t + dt;
    i = i + 1;
    //std::cout << i << std::endl;

  }     //finish time whileloop



  //AFTER TIME T, we must record the rest of the particles
  for (size_t n = 0; n < disk.size(); n++) {
    save(disk[n], tofile, t);
    disk.erase(disk.begin() + n);   //when a particle is deleted
    n = n - 1;               //subtract 1 from the index, so the next particle isn't skipped

  }



}
/////////////////////////////

void manevolve(std::vector<Particle> disk, float dt, float T = 2e5){
//This will evolve all particles for a predetermined time, N*dt
  float t = 0;
  int   i = 0;
  while (t < T) {    //begin time loop

    //std::cout <<"t"<< i << "\t";
    //disk[0].preport();              //report just the first particle at each time


//begin looping over all particles for given time
    for (size_t n = 0; n < disk.size(); n++) {
 
      disk[n] = RKupdate(disk[n], dt);    //update particle

    }

  t = t + dt;
  i = i + 1;
  }
}




//Evolve should take in dt, disk vectorlist
//RKupdate should just take in one particle, dt


//////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////



int main() {

//USER INPUT


  float dt;
  //std::cout << "dt: ";
  //std::cin >> dt;
  //Total time is 2e5
  dt = 1;

  //we want 1% error. So find the position it converges to by decreasing delta t sufficiently,

  std::vector<Particle> disk;
  std::ifstream infile("numbers.dat");

  disk  = read(infile);

///////////////////
  std::ofstream outf("output.dat");
  outf << "# 1 m    2 rx  3 ry   4 rz   5 vx     6  vy  7  vz     8 ye " << std::endl;
  

  //std::cout << disk[0].v << ", " << disk[0].r << "velocity" <<
  //disk[0].vcheck << std::endl;
  //std::cout << "dot"<< disk[0].r.dot(disk[0].v)<< std::endl;
  


  //autevolve(disk, dt, outf);
  //manevolve(disk, dt);


  return 0;
}


//END OF PROGRAM            /////////////////////////////////////////

//TO DO NEXT:
//RKupdate
//Find a way to estimate error/accuracy
//turn into array


//////////////////////





////////////////////////////////////////////
