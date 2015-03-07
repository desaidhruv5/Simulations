#include <iostream>
#include <cmath>
#include <vector>
using namespace std;


//given x, y, z, string, find the values and print to screen (write a function)
//in test program, call the function in main 



void report(string ID2, float p2, float x2, float y2, float z2) {
		cout << "pID: " << ID2
		<< ", pm = " << p2
		<< ", px = " << x2
		<< ", py = " << y2
		<< ", pz = " << z2 << endl;
}



class Force {		//defining objects of type 'Force', which has 3 components (force directions)
	float fx;
	float fy;
	float fz;
	public:
		void report() {		//
		cout << "Force = (" << fx << ", " << fy << ", " << fz << ")" << endl;
	}
	
	Force(float a, float b, float c) {
		fx = a;
		fy = b;
		fz = c;
	};

	friend Force operator+(Force f1, Force f2);

};


Force operator+(Force f1, Force f2) {
	
	float x = f1.fx + f2.fx;
	float y = f1.fy + f2.fy;
	float z = f1.fz + f2.fz;
	Force f (x, y, z);
	
	return f;
}



class Particle {

	private:
	string pID;
	float px;
	float py;
	float pz;
	float pm;


	public:

//can we bring private defs inside public?? should work but doesn't?...



	Particle (string id, float p, float x1, float y1, float z1) {
		pID = id;
		px  = x1;
		py  = y1;
		pz  = z1;
		pm = p;
	}

	void report() {
		::report(pID, pm, px, py, pz);
	}

	Force forceat(float tx, float ty, float tz){
		//compute force in x direction of this particle at pos. x = x3, y3, z3
		const int G = 1;
		float rx = px - tx;
		float ry = py - ty;
		float rz = pz - tz;
		float r = sqrt(rx*rx + ry*ry + rz*rz);
		float M = G*pm/(r*r*r);
		//know mass, position
		float Mx = M*rx;
		float My = M*ry;
		float Mz = M*rz;
		Force f(Mx, My, Mz);
		cout << "The force at (" << tx << ", "<< ty << ", " << tz << ") is:" << endl;
		return f;
	}

};



Force forceat(float tx, float ty, float tz, vector<Particle> vp) {		//finds force at tx,ty,tz due to 2 particles
	Force f1 = vp[0].forceat(tx, ty, tz);
	Force f2 = vp[1].forceat(tx, ty, tz);
	Force f = f1 + f2;
	return f;

}




int main() {

	string ID1 = "Earth";
	float ma1 = 5;		//info for particle 1
	float i = 2;
	float j = 3;
	float k = 4;

	string ID2 = "Venus";
	float ma2 = 3;		//info for particle 2
	float l1 = -2;
	float l2 = -2;
	float l3 = 2;


	float h1 = 3;		//location at which force is measured
	float h2 = 3;
	float h3 = 3;

	Particle earth(ID1, ma1, i, j, k); //create object named earth of class Particle 
	earth.report();						//calling function 'report' inside class Particle



	Force f1 = earth.forceat (h1, h2, h3);		//declaring force object 
	f1.report();
					//calling function 'report' inside class Particle

	Particle venus(ID2, ma2, l1, l2, l3); //create object named earth of class Particle 
	venus.report();						//calling function 'report' inside class Particle




vector<Particle> vp;		//declaring vector named vp with elements of type Particle
vp.push_back (earth);		//1st element called earth
vp.push_back (venus);

	for (int n = 0; n < vp.size(); n++)
        vp[n].report();
    forceat(h1, h2, h3, vp).report();
	return 0;
}

