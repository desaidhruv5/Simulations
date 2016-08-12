#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include "Eigen_lib/Eigen/Dense"
#include <cstdio>
#include <time.h>
//CONVERGENCE TESTER

double MASS = 6.0552018;
double G = 1;



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
  bool vcheck = 0;        //1 if particle is moving outwards, otherwise 0
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
  std::cout << "There are " << tail.size()<<" particles." << std::endl;
  return tail;
}





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



std::vector<Particle> manevolve(std::vector<Particle> tail, double tol, double T = 2e5)
{

  double t = 0;
  int i=0;
  while (t < T)
  {    //begin time loop
    ++i;

    //  REPLACE T WITH T[i], being the respective density stops

    //std::cout <<"t"<< i << "\t";
    //tail[0].preport();              //report just the first particle at each time


//begin looping over all particles for given time

    t=T;
    for (size_t n = 0; n < tail.size(); n++)                          //for each particle,
    {
      //pass statement, for particles ready to compute density or those that have reached end of simulation
      //adjust time step before or after RK4?
      if (tail[n].t<T)
      {
        tail[n].dt = tol*2*M_PI*pow(tail[n].r.norm(), 1.5)/sqrt(MASS);  //calculate minimum time step for given tolerance
        if (tail[n].dt + tail[n].t > T)                                 //if time step + current time exceeds final time,
        {
          tail[n].dt = T- tail[n].t;                                    //force time step to be such that particle ends at final time
          //std::cout << "tail[n].dt: "<< tail[n].dt << std::endl;
        }
        else
        {
          t = 0;
          i = i-1;                                                        //else, still more evolution to go.
        }
        RKupdate(tail[n]);    //update particle

      }
      
    }


  }
  return tail;
}



void converge_data(double tol, double error, std::ofstream & outf)
{
  //std::cout << "dt: "<<dt <<", r_diff: " << error << std::endl;
  outf << tol<<" "<< error<< std::endl;
}


int main()
{
//USER INPUT

  std::ofstream convergeplot("convergeplot.dat");

  double pindex = 20;    //index of particle for which convergence test is run

  //we want 1% error. So find the position it converges to by decreasing delta t sufficiently

  std::vector<Particle> tail;
  std::vector<std::vector<Particle>> final_tail;



  std::ifstream infile("converge.dat");

  tail  = read(infile);

///////////////////
  std::vector<double> tol;
  std::vector<double> r_diff;
  for (int j = 1; j < 11; ++j)
  {
    tol.push_back(0.0001*j);
  }

  //std::cout << "pow(.5,2) " << pow(.5,2) << std::endl;
  //tail[0].vcheck << std::endl;
  double T = 2.e5;
  final_tail.push_back(manevolve(tail, tol[0], T) );
  for (int j = 1; j < tol.size(); ++j)
  {
    final_tail.push_back(manevolve(tail, tol[j], T) );
    std::cout << "final dt of 0th particle in list: "<< final_tail[j][pindex].dt << std::endl;

    double abs_diff =(final_tail[j][pindex].r - final_tail[0][pindex].r).norm(); 
    r_diff.push_back(abs_diff/final_tail[0][pindex].r.norm());
    converge_data(tol[j-1], r_diff[j-1], convergeplot);
  }


  double NUM = 100;
  std::vector<double> dtime;
  double base = pow(T, 1/NUM);
  int exp=0;
  while (pow(base, exp) < T)
  {

    dtime.push_back(pow(base,exp));
    std::cout << "dtime["<< exp <<"] = "<< dtime[exp]<< std::endl;
    ++exp;
  }
  std::cout << "T = "<< T<< std::endl;


	return 0;
}


///////////////////////////
//plot log error vs. log dt




