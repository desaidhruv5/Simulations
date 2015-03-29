#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include "Eigen/Eigen/Dense"

//CONVERGENCE TESTER



class Particle {    //this is a vector of particles
public:
  float m;
  Eigen::Vector3f r;
  Eigen::Vector3f v;
  float u;
  float ye;
  bool vcheck;



  Particle (float m_, Eigen::Vector3f r_, Eigen::Vector3f v_, float u_, float ye_) :    //constructor for object of Particle class
            m(m_),    r(r_),                           v(v_),     u(u_),     ye(ye_){}



  void
  preport() {                              //report particle info
    std::cout << "Mass:\t" << m
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


//COMPUTING FIELD at position r away from black hole, which is at origin
Eigen::Vector3f fieldon(Particle p) {
  const float G = 1;
  const float M = 1;    //mass of black hole, only source of gravitational field

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


Particle manevolve(std::vector<Particle> disk, float N, float T = 2*3.14159265358979323846264) {
//This will evolve all particles for a predetermined time, N*dt

  float t = 0;
  int   i = 1;
  float dt = T/N;
  disk[0].preport();
  Particle p = disk[0];

  while (t <= T) {    //begin time loop

  	p = RKupdate(p, dt);    //update particle

    //std::cout <<"t"<< i << "\t";
    //disk[0].preport();              //report just the first particle at each time

  	t = t + dt;
  	i = i + 1;
  }

  return p;

}


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

      ++i;
    }
    
  }


  return disk;
}




int main(){
//USER INPUT


  float N;
  //std::cout << "dt: ";
  //std::cin >> dt;
  //Total time is 2e5
  N = 9e3;
  float T = 2*3.14159265358979323846264;

  //we want 1% error. So find the position it converges to by decreasing delta t sufficiently

  std::vector<Particle> disk;
  std::ifstream infile("test.dat");

  disk  = read(infile);

///////////////////
  

  //std::cout << disk[0].v << ", " << disk[0].r << "velocity" <<
  //disk[0].vcheck << std::endl;
  //std::cout << "dot"<< disk[0].r.dot(disk[0].v)<< std::endl;


  std::cout << "final position: " << manevolve(disk, N).r << std::endl;

  Eigen::Vector3f p = manevolve(disk, N).r;

  Eigen::Vector3f r(1, 0 ,0);
  p = p - r;
  float d = p.norm(); //d = error, since d/radius = d/1 = d

  std::cout << "N: " << N << std::endl;
  std::cout << "dt: " << T/N << std::endl;

  std::cout << "error: " << d << std::endl;

	return 0;
}






