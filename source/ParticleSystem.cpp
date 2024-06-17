/*MIT License

* Copyright (c) [2024] [Vishu Saini]

* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:

* The above copyright notice and this permission notice shall be included in all
* copies or substantial portions of the Software.

* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*/

#include "../include/ParticleSystem.h"
#include <cmath>
#include <iostream>
#include <iterator>
#include <math.h>
#include <random>
#include <stdlib.h>

/*Class Particle definitions ------------*/
ParSim::Particle::Particle() { // default constructor
  x = 0;
  y = 0;

  // forces
  force_radial[0] = 0;
  force_radial[1] = 0;
  force_tangential[0] = 0;
  force_tangential[1] = 0;

  // pressure
  sigma[0] = 0.0;
  sigma[1] = 0.0;
  sigma[2] = 0.0;
  sigma[3] = 0.0;

  // previous step information
  position_prev[0] = 0;
  position_prev[1] = 0;
  velocity_prev[0] = 0;
  velocity_prev[1] = 0;

  force_radial_prev[0] = 0;
  force_radial_prev[1] = 0;
  force_tangential_prev[0] = 0;
  force_tangential_prev[1] = 0;
}

ParSim::Particle::Particle(int N, double omega, double L) {
  Lattice_initialize(N, omega, L);
}

ParSim::Particle::Particle(double x_cor, double y_cor, double v_x, double v_y,
                           double orientation) {
  x = x_cor;
  y = y_cor;
  vx = v_x;
  vy = v_y;
  alpha = orientation;
}

void ParSim::Particle::random_initialize(int N, double phi, double L) {

  // std:: cout << "Random init called with: (N,phi,L)" << N<< " "<<phi <<" " <<
  // L << std::endl;
  std::random_device rd;
  std::uniform_real_distribution<double> x_coordinate(-1, 1);
  std::uniform_real_distribution<double> y_coordinate(-1, 1);
  std::uniform_real_distribution<double> vx_dist(-1, 1);
  std::uniform_real_distribution<double> vy_dist(-1, 1);
  std::uniform_real_distribution<double> alpha_dist(-1, 1);
  std::uniform_real_distribution<double> theta_dist(-1, 1);
  std::uniform_real_distribution<double> omega_dist(-1, 1);

  x = 1 * (L / 2) * x_coordinate(rd);
  y = 1 * (L / 2) * y_coordinate(rd); // random distribution

  // Generate random particle speed. Speed is squared causing
  // particle distribution to be exponential instead of linear.
  vx = 0.4 * vx_dist(rd);
  vy = 0.4 * vy_dist(rd);

  // Generate random particle orientation (0 to 2pi) and omegas
  alpha = alpha_dist(rd);
  omega = 0 * M_PI * omega_dist(rd);

  // Generate random V0
  vx_activity = 0 * vx_dist(rd);
  vy_activity = 0 * vy_dist(rd);

  // Generate random theta
  theta = 0 * M_PI * theta_dist(rd);

  // Generatoe random omega
  omega_activity = 0 * M_PI * omega_dist(rd);
}

void ParSim::Particle::Lattice_initialize(int N, double omega_act, double L) {

  // std:: cout << "Random init called with: (N,phi,L)" << N<< " "<<phi <<" " <<
  // L << std::endl;

  double L_inside_box = L - 2.2; // length of initialization box; 2.2 should be
                                 // replaced with 1.1xsigma
  double spacing = (L_inside_box) / std::sqrt(N); //
  static double x_cor = 0;                        // executed only once
  static double y_cor = 0;
  std::random_device rd;
  std::uniform_real_distribution<double> vx_dist(-1, 1);
  std::uniform_real_distribution<double> vy_dist(-1, 1);
  std::uniform_real_distribution<double> alpha_dist(-1, 1);
  std::uniform_real_distribution<double> theta_dist(-1, 1);
  std::uniform_real_distribution<double> omega_dist(-1, 1);

  // std:: cout << lattice_grid_ << ' ' << lattice_grid[1] << std::endl;
  if (x_cor > sqrt(N) - 1) {
    x_cor = 0;
    y_cor += 1;
  }

  // std::cout << x_cor << ' ' << y_cor << std::endl;

  // lattice grid distribution
  x = (spacing)*x_cor - 0.5 * (L - 2.2);
  y = (spacing)*y_cor - 0.5 * (L - 2.2);

  // Generate random particle speed.
  vx = 0.4 * vx_dist(rd);
  vy = 0.4 * vy_dist(rd);

  // Generate random particle orientation (0 to 2pi) and omegas
  alpha = 2 * M_PI * alpha_dist(rd);
  omega = 0 * M_PI * omega_dist(rd);

  // Generate random V0
  vx_activity = 0 * vx_dist(rd);
  vy_activity = 0 * vy_dist(rd);

  // Generate random theta -- active velocity director
  theta = 2 * M_PI * theta_dist(rd);

  // Generate random omega0
  omega_activity = omega_act; // 0.01--anticlockwise rotation

  // std:: cout << omega_activity << std:: endl;

  x_cor = (x_cor + 1);
}

void ParSim::Particle::circ_initialize(int N, double phi, double L) {

  // std:: cout << "Random init called with: (N,phi,L)" << N<< " "<<phi <<" " <<
  // L << std::endl;

  double spacing =
      (L - 2.2) / std::sqrt(N); // 2.2 should be replaced with 1.1xsigma
  static double x_cor = 0;      // executed only once
  static double y_cor = 0;
  std::random_device rd;
  std::uniform_real_distribution<double> vx_dist(-1, 1);
  std::uniform_real_distribution<double> vy_dist(-1, 1);
  std::uniform_real_distribution<double> alpha_dist(-1, 1);
  std::uniform_real_distribution<double> theta_dist(-1, 1);
  std::uniform_real_distribution<double> omega_dist(-1, 1);

  // std:: cout << lattice_grid_ << ' ' << lattice_grid[1] << std::endl;
  if (x_cor > sqrt(N) - 1) {
    x_cor = 0;
    y_cor += 1;
  }

  // std::cout << x_cor << ' ' << y_cor << std::endl;

  // lattice grid distribution
  x = (spacing)*x_cor - 0.5 * (L - 2.2);
  y = (spacing)*y_cor - 0.5 * (L - 2.2);

  // Generate random particle speed.
  vx = 10 * vx_dist(rd);
  vy = 10 * vy_dist(rd);

  // Generate random particle orientation (0 to 2pi) and omegas
  alpha = 2 * M_PI * alpha_dist(rd);
  omega = 0 * M_PI * omega_dist(rd);

  // Generate random V0
  vx_activity = 0 * vx_dist(rd);
  vy_activity = 0 * vy_dist(rd);

  // Generate random theta -- active velocity director
  theta = 2 * M_PI * theta_dist(rd);

  // Generate random omega0
  omega_activity = 100; // 0.1--anticlockwise rotation

  x_cor = (x_cor + 1);
}

/*Class Particle System definitions----------------*/
ParSim::ParticleSystem::ParticleSystem(int N, double omega, double dim) {
  this->no_of_particles = N;
  this->particle_array = new Particle[no_of_particles];
  this->L = dim;

  for (int k = 0; k < N; k++) {

    Particle temp(N, omega, L);
    particle_array[k] = temp;

    // change position of kth particle here
  }
}

ParSim::ParticleSystem::~ParticleSystem() { delete[] particle_array; }

namespace ParSim {
Particle *const ParticleSystem::get_particles() { return particle_array; }
} // namespace ParSim

double ParSim::ParticleSystem::distance(Particle &par1, Particle &par2) {
  double distance = sqrt((par1.x - par2.x) * (par1.x - par2.x) +
                         (par1.y - par2.y) * (par1.y - par2.y));
  return distance;
}

double ParSim::ParticleSystem::min_sep(double &x1, double &x2) {
  double dx = x1 - x2;
  if (dx < (-L / 2)) {
    dx += L;
  } else if (dx > (L / 2)) {
    dx -= L;
  }

  return dx;
}

double ParSim::ParticleSystem::nearest_img_dist(Particle &par1,
                                                Particle &par2) {

  double dist;
  double a = min_sep(par1.x, par2.x);
  double b = min_sep(par1.y, par2.y);
  dist = sqrt(a * a + b * b);

  return dist;
}

double ParSim::ParticleSystem::nearest_img_dist_wall_y(Particle &par1,
                                                       Particle &par2) {

  double dist;
  double a = min_sep(par1.x, par2.x);
  dist = sqrt(a * a + (par1.y - par2.y) * (par1.y - par2.y));

  return dist;
}

double ParSim::ParticleSystem::dist_from_origin(Particle &par) {
  double distance;
  distance = sqrt(par.x * par.x + par.y * par.y);
  return distance;
}
