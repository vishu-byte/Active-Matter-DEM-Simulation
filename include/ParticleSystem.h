//
// Created by Vishu Saini on 31/08/23
//

#ifndef PARTICLE_SIMULATION_PARTICLE_SYSTEM_H
#define PARTICLE_SIMULATION_PARTICLE_SYSTEM_H

#include <cmath>
#include <math.h>
#include <random>
#include <vector>
#include <list>

namespace ParSim { // for particle simulation

class Particle {

public:
  // cartesian coordinates
  double x;
  double y;
  double vx;
  double vy;
  double alpha;
  double omega;
  double radius;
  double vx_activity;
  double vy_activity;
  double theta;
  double omega_activity;
  double force_radial[2];      //radial forces fx fy
  double force_tangential[2];   //tangential forces fx fy
  double torque{0};                           // torque

  //virial pressure
  double sigma[4];

  // For storing forces of one step before
  double position_prev[2];   // (x, y)
  double velocity_prev[2];   // (vx, vy)
  double alpha_prev;
  double omega_prev;

  double force_radial_prev[2];    // radial force: fx and fy
  double force_tangential_prev[2]; 
  double torque_prev{0};                           // torque
  // std:: vector <double> IC{0,0,0,0,0};

  std:: list<int> Neighbours;   // Neighbour list of zero intial length



  Particle();                    // default constructor
  Particle(int, double, double); // parameterized constructor: int N, double phi
  Particle(double x_cor, double y_cor, double speed, double direction,
           double orientation);                 // parameterized constructor
  virtual ~Particle(){};                        // virtual destructor
  void random_initialize(int, double, double);  // randomly initializes particle
  void Lattice_initialize(int, double, double); // lattice initializes particle
  void circ_initialize(int, double, double);    // lattice initializes in circle
};

class ParticleSystem {
public:
  int no_of_particles;
  int n; // no. of particles in a row
  double L;
  double phi;

  std::vector<double> lattice_grid{0, 0};
  Particle *particle_array{nullptr}; // creating particle array on heap

  ParticleSystem(int, double, double);      // constructor
  ParticleSystem(int, int, double, double); // parameterized constructor
  virtual ~ParticleSystem();                // destructor
  Particle *const get_particles(); // constant pointer, can not change address
                                   // of memory block to which it points
  double distance(Particle &par1, Particle &par2);

  double min_sep(double &x1, double &x2);

  double nearest_img_dist(Particle &par1, Particle &par2);

  double nearest_img_dist_wall_y(Particle &par1, Particle &par2);

  double dist_from_origin(Particle &par);
};

class Boxresize {
public:
  double L;
  double phi;

public:
  Boxresize(){};
};

} // namespace ParSim

#endif // PARTICLE_SIMULATION_PARTICLE_AND_SYSTEM_H
