#include <iostream>
#include <cmath>
using namespace std;
//3 inputs
// a x^2 +bx =c
// return value of x

//declare class called KS
// no member data, one method called ksolveit



/*
double * ksolveit (double a, double b, double c) {
   double * x = new double[2];
	x[0] = (-b + sqrt(b*b + 4*a*c))/(2*a);
    *(x+1) = (-b - sqrt(b*b + 4*a*c))/(2*a);
    return x;
}
*/


class KS {
public:
	double * solveit (double a, double c=0)
	{
    double b = 4;
  	double * x = new double[2];
	x[0] = (-b + sqrt(b*b + 4*a*c))/(2*a);
    *(x+1) = (-b - sqrt(b*b + 4*a*c))/(2*a);
    return x;
	}


};



double * ksolveit(double a, double b, double c) {
	KS ks;
    double * xks = ks.solveit(a,c);
    return xks;
	}


int main() {
	double a;
 	a        = 1;
 	double b = 0;
 	double c = -4;
	double * kx = ksolveit(a,b,c);
	std::cout << "When a = " << a ;
	std::cout << ", b = " << b ;
	std::cout << ", and c = " << c << std::endl;
    std::cout << "kx = " << kx[0] << " " << *(kx+1) << std::endl;
	delete kx;
	
    KS	ks;
    double * xks = ks.solveit(a, c);
    std::cout << "xks = " << xks[0] << " " << *(xks+1) << std::endl;
	delete xks;

  	double * xks2 = ksolveit(a,b,c);
    std::cout << "xks = " << xks2[0] << " " << *(xks2+1) << std::endl;
    delete xks2;
	

	return 0;
}