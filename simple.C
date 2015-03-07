#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include "Eigen/Eigen/Dense"




float check (Eigen::Vector3f v, int a = 6) {

	bool b;

	v[0] = v[0] + a;
	return v[0];
}


int main(){

	Eigen::Vector3f v12(58.5, 40.2, 5.54); 
	Eigen::Vector3f  v1(59,   39.8, 5.35); 
	Eigen::Vector3f  v2(63.4, 38.3, 4.68);

	float h = (v1 - v12).norm();
	float g = (v2 - v12).norm();



	std::cout << h/g <<std::endl;
	std::cout << "...should be: " << 1/9. << std::endl;


	return 0;
}