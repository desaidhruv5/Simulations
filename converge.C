#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include "Eigen/Eigen/Dense"

//CONVERGENCE TESTER

//should loop through all particles, checking convergence of each one.
//should get a ratio for each particle


int main(){
    //these below three vectors are the end positions of the 1st particle in the "numbers.C" list,
    //with a time step of 1/2, 1, and 2 respectively, after evolving for time = 2e5.
	Eigen::Vector3f v12(58.5, 40.2, 5.54); 
	Eigen::Vector3f  v1(59,   39.8, 5.35); 
	Eigen::Vector3f  v2(63.4, 38.3, 4.68);

	float h = (v1 - v12).norm();
	float g = (v2 - v12).norm();



	std::cout << h/g <<std::endl;
    //for 4th order, the ratio should be 1/17, for 3rd order, ratio should be 1/9
	std::cout << "Should be: " << 1/17. << std::endl;


	return 0;
}