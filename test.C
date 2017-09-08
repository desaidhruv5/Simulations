#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include "Eigen_lib/Eigen/Dense"
#include <cstdio>
#include <omp.h>





class Particle {      //this is a vector of particles
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
  bool dflag = 0;     //set to 1 for particles that we track densities for


  Particle (double m_, Eigen::Vector3d r_, Eigen::Vector3d v_) :    //constructor for object of Particle class
               m(m_),               r(r_),              v(v_){}



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
//--------------------------------------------------------------------------------------



void set_dflags(std::vector<Particle> & disk){

  for (int i = 0; i < int(disk.size() ); ++i)
  {
    disk[i].dflag=1;
  }

}



bool compdist (Particle p1, Particle p2, Particle center)
{ //comparing r1-r to r2-r
  if ( (p1.id != center.id) && (p2.id != center.id) )
  {
    return ( (p1.r-center.r).norm() < (p2.r-center.r).norm() );  //This returns 1 if r1-r < r2-r
  }

  else {
    return 1;
  }
}

std::vector<Particle> v;

void evolve (std::vector<Particle> & v)
{


  double m0 = 0;
  double m1 = m0;
  double m2 = m0;
  Eigen::Vector3d r0(0,0,0);
  Eigen::Vector3d r1(1,0,0);
  Eigen::Vector3d r2(2,0,0);

  Eigen::Vector3d v0(0,0,0);
  Eigen::Vector3d v1(0,0,0);
  Eigen::Vector3d v2(0,0,0);

  Particle p0(m0,r0,v0);
  Particle p1(m1,r1,v1);
  Particle p2(m2,r2,v2);

  p0.id = 0;
  p1.id = 1;
  p2.id = 2;



  v.push_back(p0);
  v.push_back(p2);
  v.push_back(p1);
  
  auto compare = [p0] (Particle part1, Particle part2)
  {
    return compdist(part1, part2, p0);
  };

  std::cout<< "Before: center.r=" << p0.r <<std::endl;


  std::partial_sort(v.begin(), v.begin()+3, v.end(), compare);


  std::cout<< "After: v0=" << v[0].r << ", v1=" << v[1].r <<", v2=" << v[2].r <<std::endl;

  std::cout<< "After: center.r=" << p0.r <<std::endl;
  std::cout<< "dflag=" << v[0].dflag <<std::endl;

  set_dflags(v);


  std::cout<< "dflag=" << v[0].dflag <<std::endl;

  std::vector<Particle> a = v;
  std::cout<< "Before: a=" << a[0].r <<std::endl;

  v[0].r[0]=10000;
  std::cout<< "After: a=" << a[0].r <<std::endl;




}



#define sq(x) x*x

#define cu(x) x*x*x

#define qu(x) x*x*x*x




std::vector<std::string> split(std::string str, char delimiter)
{
  std::vector<std::string> internal;
  std::stringstream ss(str); // Turn the string into a stream.
  std::string tok;
  
  while(getline(ss, tok, delimiter)) {
    internal.push_back(tok);
  }
  
  return internal;
}



int main ()
{
  
  Eigen::Vector3d r1(1.10000000000001,0,0);

  //std::cout << r1[0] << std::endl;
  //std::cout << r1[0]-1.1 << std::endl;
  std::ofstream fallbackdat("test.dat");
  fallbackdat << "# m << x << y << z << vx << vy << vz << u << rho << temp << ye << e << t" << std::endl;
  fallbackdat << "# m << x << y << z << vx << vy << vz << u << rho << temp << ye << e << t" << std::endl;
  fallbackdat << "# m << x << y << z << vx << vy << vz << u << rho << temp << ye << e << t" << std::endl;
  std::ifstream inputfile("test.dat");


  std::string line;
  std::getline(inputfile, line);

  std::vector<std::string> v = split(line, ' ');

  int i=0;
  while(i<int(v.size()) )
  {
    std::cout << v[i] <<std::endl;
    ++i;
  }

  int a =4;
  goto AA;
  std::cout.precision(17);
  std::cout << sq(4)<<" "<<cu(3)<<" " << fmax(4.0010212123,4.0010212124) <<std::endl;
  AA: ;

  return 0;
}


















