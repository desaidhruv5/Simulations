#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include "Eigen/Eigen/Dense"

//////////////////////////////////////////////////////////////////*


//checks direction of particle velocity
//returns 1 if radial velocity is positive, else returns 0
bool check (Eigen::Vector3f r, Eigen::Vector3f v) {


  bool b;
  if (r.dot(v) > 0) b =1;
  else b =0;

  return b;
}
///////////////



class Particle {    //this is a vector of particles
public:
  float m;
  Eigen::Vector3f r;
  Eigen::Vector3f v;
  float ye;
  bool vcheck;



  Particle (float m_, Eigen::Vector3f r_, Eigen::Vector3f v_, float ye_) :    //constructor for object of Particle class
            m(m_),    r(r_),                 v(v_),                 ye(ye_){}



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
  Particle p(m, r, v, ye);

  int i = 0;
 

  while (!infile.eof()) {
        //peek at first character.   
    if (infile.peek() == '#') infile.ignore(256, '\n');    //is this "#"? If so, then ignore until next line
      
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
      disk[i].ye = ye;                //disk[i] is ith particle, ye is composition of particle
      disk[i].vcheck = check (disk[i].r, disk[i].v); //returns 0 if radial velocity is negative, else returns 1

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

  p.r = p1.r + .5*dt*p.v;
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


//void save(Particle p, **changetoutputstream**std::ifstream & tofile) {
  
//}

/*
void autevolve(std::vector<Particle> disk, float dt, float T = 2e5) { //This will auto-stop
 
  float t = 0;
  int i = 0;
  while (t < T) { //begin time loop

    //std::cout <<"t"<< i << "\t";
    disk[0].preport();           //report position of just the first particle

    Eigen::Vector3f dv;          //create update vectors
    Eigen::Vector3f dr;


      //need if statements -- if particle falls closer than <100km,
      //then record the time, and stop updating that particle

      //to find periapse, need to save closest position. if r_i^2 > r_(i+1)^2 then peri =r_(i+1). else, peri= r_i

    for (size_t n = 0; n < disk.size(); n++) { //looping over all particles

      if (check(disk[n].r, disk[n].v) == 1) {   //enter if CURRENT radial velocity is positive

        if (disk[n].vcheck == 0) {   //enter this case only if velocity ever used to be negative
          //save particle data, delete from list, for the case of a now positive vr, which was previously negative
          //when a particle is deleted, subtract 1 from the index, so the next particle isn't skipped
          //n = n - 1;
        }


        else {
          RKupdate(disk[n], dt);    //if within range but moving outwards, keep updating
        }

      }


      else {
        disk[n].vcheck = 0;       //affirms that the particle _has_ in the past had negative radial velocity
        RKupdate(disk[n], dt);    //update for the case of falling particle
      }

      t = t + dt;
      i = i + 1;
    }//end looping over particles

  } //end time whileloop

  //now we must record the rest of the particles
  for (size_t n = 0; n < disk.size(); n++) {
    //save(disk[n], tofile);
    //n = n - 1;
  
  }
}
*/////////////////////////////


void manevolve(std::vector<Particle> disk, float dt, float T = 2e5){
//This will evolve all particles for a predetermined time, N*dt
  float t = 0;
  int   i = 0;
  while (t < T) {    //begin time loop

    std::cout <<"t"<< i << "\t";
    disk[0].preport();              //report just the first particle at each time

    for (size_t n = 0; n < disk.size(); n++) {   //begin looping over all particles for given time
 
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
  std::cout << "dt: ";
  std::cin >> dt;

  //we want 1% error. So find the position it converges to by decreasing delta t sufficiently,
  //then 

  
  std::vector<Particle> disk;
  std::ifstream infile("numbers.dat");
  disk  = read(infile);

///////////////////

  manevolve(disk, dt);

  return 0;
}





//END OF PROGRAM            /////////////////////////////////////////

//TO DO NEXT:
//change to:  save to file, initial position and velocity
//RKupdate
//Find a way to estimate error/accuracy
//turn into array





//////////////////////





////////////////////////////////////////////
