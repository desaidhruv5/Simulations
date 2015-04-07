#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include "Eigen/Eigen/Dense"



class Particle {    //Defining object 'particle'
public:
float m;
std::vector<float> r;
std::vector<float> v;
float ye;


  Particle (float m_, std::vector<float> r_, std::vector<float> v_, float ye_) :    //constructor for object of Particle class
            m(m_),    r(r_),                 v(v_),                 ye(ye_){}

};





std::vector<Particle> read(std::ifstream & infile) {

//STEPS:
//reads each line, 1: creates object of type particle, 2:labels it, 3: assigns position, 4: assigns velocity, 5:

  //here we define what each element on each line represents
  float m, x, y, z, vx, vy, vz, u, ye, e, rho;

  std::vector<float> r;
  std::vector<float> v;
  std::vector<Particle> disk;
  Particle p(m, r, v, ye);




  int i = 0;



    while (!infile.eof()) {
        //peek at first character.   
      if (infile.peek() == '#') infile.ignore(256, '\n');    //is this "#"? If so, then ignore until next line
      
      else {
      infile >> m >> x >> y >> z >> vx >> vy >> vz >> u >> ye >> e >> rho;
      disk.push_back(p);               //creating entry for ith particle in the disk
      disk[i].m = m;                  //mass of ith particle in the disk
      disk[i].r.push_back(x);         //ith particle of the disk, position vector, create x component
      disk[i].r.push_back(y);         //...
      disk[i].r.push_back(z);
      disk[i].v.push_back(vx);
      disk[i].v.push_back(vy);
      disk[i].v.push_back(vz);
      disk[i].ye = ye;                //disk[i] is ith particle, ye is composition of particle
      //std::cout << disk[i].m <<std::endl;
      ++i;
      }
    
    }

  return disk;
}



int count(std::ifstream & infile) {
	int i = 0;
	std::string unused;
	
	while ( std::getline(infile, unused) )
   	++i;

  	return i;		//so there are i lines
}








int main() {

/*
	std::ifstream infile("numbers.txt");

	float m, x, y, z, vx, vy, vz, ye;

	while (infile >> m >> x >> y >> z >> vx >> vy >> vz >> ye) {
		float c(m);
		float d(x);
		float e(y);
		std::cout << c << ", " << d << ", " << e << std::endl;
	}
*/
	std::vector<Particle> disk;
	std::ifstream infile("numbers.dat");
	disk  = read(infile);
	std::cout<< infile.is_open() << std::endl;
	infile.close();
	std::ifstream infile1("numbers.dat");
	std::cout << count(infile1) << std::endl;


	std::cout << disk[0].m << ", " << disk[1].m << ", " << disk[2].m <<  ", " << disk[0].ye << std::endl;

	return 0;
}












