#include <iostream>
#include <cmath>
#include <vector>
#include <string>

//////////////////////////////////////////////////////////////////

class Time {                            //Defining Time class

public:

  int t;
  Time(float t_) : t(t_) {}

};

class Triple {

public:

  float x;
  float y;
  float z;
  Triple(float x_, float y_, float z_) :
          x(x_), y(y_), z(z_) {};

  //Triple(const Triple & r) : x(r.x), y(r.y), z(r.z) {};

  Triple                //Declare method here
  operator+ (const Triple & );

  // kmd: Displacement operator-(Position)
  Triple
  operator- ( const Triple & );

  Triple
  operator* (const float &);

  double
  operator* (const Triple &);

  //Triple
  //operator* (const Time &);

  //declare methods, but don't define them

};



class Displace : public Triple {          //Defining Displace class
public:
  Displace(Triple fooey) : Triple(fooey) {};

};


class Position : public Displace {        //Defining Position class

public:
  Position(Triple a) :    //constructor for object of type position
          Displace(a) {};         // order of initializers must match member data

  void
  report() const;

  Displace
  operator-(const Position& ) const;

};


class Velocity : public Triple {        //Defining Velocity class
public:

  Velocity(Triple v_) :
          Triple(v_) {};

  void
  report() const;

  Displace
  operator* (const Time &); //DEFINE

  Velocity
  operator* (const float &);

  Velocity
  operator+ (const Velocity & f);


};

//defining objects of type 'Field', which has 3 components (it's a Triple)
class Field : public Triple {           //Defining Field class
public:
  Field(Triple f_) :
        Triple(f_) {}; 

  void
  report() const;


  Velocity
  operator* (const Time &);  //DEFINE

  Field
  operator* (const float &);

  Field
  operator+ (const Field & f);


};


//////////////////////////////////////////////////////////////
Triple
Triple::operator+ (const Triple & e) {
    return Triple(x+e.x, y+e.y, z+e.z);
}

double
Triple::operator* (const Triple & e) {
    return x * e.x + (*this).y * e.y + this->z * e.z;
}


Triple
Triple::operator* (const float & m) {
    return Triple(x*m, y*m, z*m);
}


Triple
Triple::operator- (const Triple & d) {
  return Triple(x-d.x, y-d.y, z-d.z);
}


Displace
Position::operator- (const Position& e) const {
  return (Triple) *this - (Triple) e;
}


Velocity
Field::operator* (const Time& m) {
  return Velocity( Triple(x*m.t, y*m.t, z*m.t));
}


Displace
Velocity::operator* (const Time& m) {
  return Displace( Triple(x*m.t, y*m.t, z*m.t));
}


//Triple
//Triple::operator* (const Time& m){
//  return Triple(x*m.t, y*m.t, z*m.t);
//}

Field
Field::operator* (const float & m) {
  return Field(Triple(m*x, m*y, m*z));
}

Field
Field::operator+ (const Field & f) {
  return Field(Triple(f.x*x, f.y*y, f.z*z));
}



Velocity
Velocity::operator* (const float & m) {
  return Velocity(Triple(m*x, m*y, m*z));
}

Velocity
Velocity::operator+ (const Velocity & f) {
  return Velocity(Triple(f.x*x, f.y*y, f.z*z));
}

//////////////////////////////////
void
Position::report() const {
  std::cout << "The Position is: (" <<
  x << ", " << y << ", " << z
  << ")." << std::endl;
}

//////////////////////////////////

void
Velocity::report() const {
  std::cout << "The Velocity is: (" <<
  x << ", " << y << ", " << z
  << ")." << std::endl;
}

//////////////////////////////////

void
Field::report() const {
  std::cout << "Field = (" <<
  x << ", " << y << ", " << z
  << ")" << std::endl;
}
////////////////////////



////////////////////////
class Particle {

  public:
    std::string pID;  //ID of particle
    float pm;         //mass of particle
    Position r;       //position of particle
    Velocity v;       //velocity of particle
    
    Particle (std::string pID_, float pm_, Position r_, Velocity v_) :    //constructor for object of Particle class
      pID(pID_), pm(pm_), r(r_), v(v_) {} 

    Position                //defining  position
    position() const {return r;}

    Velocity
    velocity() const {return v;}  //defining particle's velocity

    bool
    operator== ( Particle ) const;


    void
    report() {                              //report particle info
      std::cout << "\t"<< pID
      << " x = " << r.x
      << " y = " << r.y
      << "\tz = " << r.z << std::endl;
    }

    Field
    fieldat(const Position & h) const{
      //compute field due to this particle at pos. h = x3, y3, z3
      
      const float G = 6.67e-11;
      
      ////CAN DELETE THIS LATER, SINCE R SHOULD NEVER EQUAL H WHEN WE EVALUATE THE FIELD


      //else {
        Displace d = r - h;
        float R = sqrt( (r-h) * (r-h));
        float M = G*pm/(R*R*R);
        //know mass, position
        Field f(d*M);
        if (0) 
          std::cout << "   field calc due to "
        << pID << " : pidx = " 
        << r.x << " atx = " << h.x << " r-h= " << d.x
        << " R = " << R << " pm = " << pm << " G = " << G
        << " M = " << M << " fx = " << f.x
        << std::endl;
        return f;
      }
    //}
};



//Defining particles being equal
bool
Particle::operator== ( Particle d) const {
  return bool(pID == d.pID);
}


///////////////////////

// Defining fieldat(h, vp) function, finds total field due to n particles at position p

Field
netfield(const Particle & p, std::vector<Particle> vp) {

  Field f(Triple(0,0,0));


  for (size_t n = 0; n < vp.size(); n++) {
    
    if (p == vp[n]) {
    }

    else {
      f = f + vp[n].fieldat(p.r);
    }

    if (0)
      std::cout << " at " << p.pID
    << ", field due to " << vp[n].pID << ". dfx = " << vp[n].fieldat(p.r).x
    << " fx = " << f.x
    << std::endl;
  }
  return f;

}




//////////////
//compute using Euler Method
void compute(std::vector<Particle> p, Time dt, int T) {
 
  for (int t = 0; t <= T; ++t) {    //at each time t,
    std::cout <<"t"<< t;

     // create a dr, dv vector to find the next state vector
    std::vector<Velocity> dv; 
    std::vector<Displace> dr;

    for (size_t n = 0; n < p.size(); n++) {   //for each particle,
      
      p[n].report();  // report the positions and velocities

      Field f = netfield( p[n] , p); //find field at p[n] due to all others
      dv.push_back( f * dt);                   //find dv
      dr.push_back( p[n].v * dt + (f*dt)*dt*0.5 );    //find dr
      
      if (0) std::cout << "  At time " << t << ", p" << n
                << " f,dv,dr = " 
                << " " << f.x
                << " " << dv[n].x 
                << " " << dr[n].x
                << std::endl;
    }
    //update this info back into particle info:
    for (size_t n = 0; n < p.size(); n++)  {
      p[n].r = p[n].r + dr[n]; 
      p[n].v = p[n].v + dv[n];
    }
  
  }
}


void RKcompute(std::vector<Particle> p, const Time & dt, const int & T) {

  for (int t = 0; t <= T; ++t) {
    std::cout <<"t"<< t;

    // compute the field at each particle position due to all the other particles
    std::vector<Velocity> dv;
    std::vector<Displace> dr;

    for (size_t n = 0; n < p.size(); n++) { 
      // report the positions and velocities of each particle
      p[n].report();

      Particle q  = p[n];
      Velocity j1 = q.v;               //velocity at initial
      Field k1    = netfield(q , p);               //field at initial
      
      q.r = p[n].r + (j1*dt)*.5;
      Velocity j2 = q.v + (k1*dt)*.5;    //velocity at midpt using field initial
      Field k2    = netfield(q , p);     //field at midpt using velocity initial

      q.r = p[n].r + (j2*dt)*.5;
      Velocity j3 = q.v + (k2*dt)*.5;    //velocity at midpt. using field 
      Field k3    = netfield(q , p);
      
      q.r = p[n].r + j3*dt;
      Velocity j4 = q.v + k3*dt;
      Field k4    = netfield( q, p);



      dv.push_back( ((k1 + k2*2 + k3*2 + k4)*dt) *(1/6.)  );       //find dv
      dr.push_back( ((j1 + j2*2 + j3*2 + j4)*dt) *(1/6.)  );       //find dr

      

      if (0) std::cout << "  At time " << t << ", p" << n
                << " f,dv,dr = " 
                << " " << k1.x
                << " " << dv[n].y 
                << " " << dr[n].y
                << std::endl;
    }
    //update this info back into particle info:
    for (size_t n = 0; n < p.size(); n++)  {
      p[n].r = p[n].r + dr[n]; 
      p[n].v = p[n].v + dv[n];
    }
  
  }
}


//////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
int main() {

//USER INPUT

  std::vector<Particle> p;    //declaring vector named vp with elements of type Particle
  
  //now, to fill the particle vector with particles:
  p.push_back (Particle("Sun", 2e30,
          Position(Triple(0, 0, 0)), Velocity(Triple(0, 0, 0)) ) );        //1st element called earth
  p.push_back (Particle("Earth", 6e24,
          Position(Triple( 1.5e11, 0, 0)), Velocity(Triple(0, 3e4, 0))  ) );
//  p.push_back (Particle("Moon", 7.3e22,
//          Position(Triple( 1.5e11, 3.84e8, 0)), Velocity(Triple(-1e3, 3e4, 0))  ) );

  Time dt = 3600*24;       //step size for computation
  int T = 365;             //number of total steps

///////////////////

  RKcompute(p, dt, T);
//compute(p, dt, T);




//END OF PROGRAM            /////////////////////////////////////////

//TO DO NEXT:
//change to:  save to file, initial position and velocity
//RKcompute
//Find a way to estimate error/accuracy





//////////////////////










////////////////////////////////////////////
  return 0;
}