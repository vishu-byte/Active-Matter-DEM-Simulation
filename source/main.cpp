#include "../include/Physics.h"
#include <bits/types/time_t.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <stdlib.h>
#include <time.h>

void state_before_simulation(std::ofstream &log, ParSim::ParticleSystem &parsym,
                             ParSim ::Physics &physics, int steps,
                             double dimension, double phi, double omega,
                             double Pecr);
void state_after_simulation(std::ofstream &log, ParSim::ParticleSystem &parsym,
                            ParSim ::Physics &physics);

int main() {

  // 1)Creating and initializaing particle system

  /*Parameters*/
  /*Try to stick to S.I units to make sense out of numbers*/
  int Number_of_particles = 25;
  int Number_of_time_steps = 100000;

  // Mips parameters
  double omega = 0.1;                              // torque  = 6
  double phi = 0.65;                            // packing fraction
  double L = sqrt(Number_of_particles / (phi)); // Radius of boundary

  ParSim::ParticleSystem parsym(
      Number_of_particles, omega,
      L); // create a simple system (using lattice init)
  ParSim::Physics Physics;

  ParSim::Particle *const particle =
      parsym.get_particles(); // get access to particles

  /*Setting physics parameters -- all game to be played here */
  Physics.parameters[8] = 0.001; // time step
  Physics.parameters[1] = 2.0;   // interaction_diameter sigma
  Physics.parameters[2] = 10.0;  // mass
  Physics.parameters[13] =
      Physics.parameters[2] * (Physics.parameters[1]) *
      (Physics.parameters[1]) / 8.0; // I = (1/8)m*(simga)^2



  Physics.parameters[5] = 1000.0;    // gamma_t
  Physics.parameters[10] = 500; // gamma_r

  Physics.parameters[4] = 0.8; // mu

  // Langevin_Dyanmics.parameters[12] = 0.0; // KbT

  Physics.WCA_parameters[0] = 1.0; // WCA epsilon
  Physics.WCA_parameters[1] =
      Physics.parameters[1]; // WCA sigma

  Physics.circular_wall_parameters[0] = 1.0; // WCA epsilon
  Physics.circular_wall_parameters[1] =
      Physics.parameters[1]; // WCA sigma

  Physics.NL_skin_depth =
      2 * Physics.parameters[1]; // 2*simga, NL_skin_depth, DeltaR
  Physics.NL_tau = 200;          // NL_Tau, tau

  /*---obsolete params---*/
  Physics.parameters[0] = 1; // k
  Physics.parameters[3] = 1; // radius

  Physics.parameters[11] =
      Physics.parameters[10]; // Dr  -- rotational diffusion

  Physics.parameters[6] = 0.00000001; // epsilon1  -- softening length
  Physics.parameters[7] =
      M_PI / 10000000; // epsilon2 -- softening omega
  Physics.parameters[9] =
      0.5 * Physics.parameters[5] /
      pow((Physics.parameters[2] * Physics.parameters[0]),
          0.5); // zeta

  double Pecr = 50; // rotational Peclet number

  /* 2)Reading initial state positions*/

  // std::ifstream input_state("init_state.txt");
  // int i = 0;
  // double temp;

  // if (!input_state) { // file couldn't be opened
  //   std::cerr << "\033[31mError: Init_state file could not be opened"
  //             << std::endl;
  //   exit(1);
  // }

  // while (input_state >> particle[i].x >>
  //        particle[i].y) { // read particle x and y and z in temp
  //   i++;
  // }

  // std::cout << "\033[32mInput state read for " << i << " particles ..."
  //           << std ::endl;

  // input_state.close();

  // 3)Creating a data file for storage and log-----------

  std::ofstream data_output;
  std::ofstream log;
  std::ofstream data_all;

  data_output.open("data.xyz");
  log.open("log.txt");
  data_all.open("data_all.txt");

  // Print the state before the simulation in log
  state_before_simulation(log, parsym, Physics, Number_of_time_steps,
                          L, phi, omega, Physics.parameters[10]);

  log << "-x-x-x-x-x-Simulation initiated-x-x-x-x-x- " << std::endl;
  std::cout << "\033[35m-x-x-x-x-x-Simulation initiated-x-x-x-x-x- "
            << std::endl;
  std::cout << "\033[36mNo. of time steps: " << Number_of_time_steps
            << std::endl;

  time_t start = time(&start); // for measuring total runtime

  unsigned limit = 200000;

  // 4) Main simulation loop--------------
  for (int step = 0; step < Number_of_time_steps; step++) {

    // writing data of this state to file (will be used for rendering the
    // system in ovito), write every nth state

    if (step % 200 == 0) {
      data_output << Number_of_particles << std::endl;
      data_output << "Lattice="
                  << "\"10.0 0.0 0.0 0.0 10.0 0.0 0.0 0.0 0.0\"" << std::endl;
      // first store current configuration
      for (int i = 0; i < parsym.no_of_particles; ++i) {
        data_output << particle[i].x << ' ' << particle[i].y << ' ' << 0 << ' '
                    << cos(particle[i].alpha) << ' ' << sin(particle[i].alpha)
                    << ' ' << 0 << ' ' << particle[i].alpha << ' '
                    << particle[i].vx << ' ' << particle[i].vy << ' '
                    << particle[i].omega << ' ' << std::endl;
      }
    }

    if (step % 200 == 0) {
      for (int i = 0; i < parsym.no_of_particles; ++i) {
        data_all << step << ' ' << i << ' ' << particle[i].x << ' '
                 << particle[i].y << ' ' << particle[i].vx << ' '
                 << particle[i].vy << ' ' << particle[i].alpha << ' '
                 << particle[i].omega << ' ' << particle[i].force_radial[0]
                 << ' ' << particle[i].force_radial[1] << ' '
                 << particle[i].force_tangential[0] << ' '
                 << particle[i].force_tangential[1] << ' ' << particle[i].torque
                 << ' ' << particle[i].sigma[0] << ' ' << particle[i].sigma[1]
                 << ' ' << particle[i].sigma[2] << ' ' << particle[i].sigma[3]
                 << std::endl;
      }
    }

    if (step % 100 == 0) {
      std ::cout << "\033[37m----------Step count: " << step << " \033[32m["
                 << step / double(Number_of_time_steps) * 100 << " %]"
                 << std::endl;
      log << "----------Step count: " << step << std::endl;
    }

    // Manipulate particle positions for next iteration.
    Physics.evolve_system_ERM_NL(parsym, step);
  }

  time_t end = time(&end);

  data_output.close();
  data_all.close();

  std::cout << "-x-x-x-x-x-Simulation ended-x-x-x-x-x-" << std::endl;
  log << "\033[35m-x-x-x-x-x-Simulation ended-x-x-x-x-x-" << std::endl;
  log << "Runtime: " << end - start << " seconds" << std::endl;

  /*----------------------------------*/

  // Print the state before the simulation in log
  state_after_simulation(log, parsym, Physics);

  return 0;
}

void state_before_simulation(std::ofstream &log, ParSim::ParticleSystem &parsym,
                             ParSim ::Physics &physics, int steps,
                             double dimension, double phi, double omega,
                             double gamma_r) {

  ParSim::Particle *const particle =
      parsym.get_particles(); // get access to paticles

  log << "-------Parameters and state------" << std::endl;
  log << "Number of particles: " << parsym.no_of_particles << std::endl
      << "Time step: " << physics.parameters[8] << std::endl
      << "Number of time steps: " << steps << std::endl
      << "Dimension: " << dimension << std::endl
      << "phi: " << phi << std::endl
      << "k: " << physics.parameters[0] << std::endl
      << "Interaction diameter (sigma): " << physics.parameters[1] << std::endl
      << "Mass (m): " << physics.parameters[2] << std::endl
      << "Omega: " << omega << std::endl
      << "gamma_t: " << physics.parameters[5] << std::endl
      << "gamma_r: " << gamma_r << std::endl
      << "mu: " << physics.parameters[4] << std::endl;

  log << "Energy-momentum before the collision: " << std::endl;
  log << "Total Energy: "
      << physics.EnergyMomentum(parsym)[0] + physics.EnergyMomentum(parsym)[1]
      << std::endl;
  log << "Translational K.Energy: " << physics.EnergyMomentum(parsym)[0]
      << std::endl;
  log << "Rotational K.Energy: " << physics.EnergyMomentum(parsym)[1]
      << std::endl;
  log << "Momentum: "
      << "(" << physics.EnergyMomentum(parsym)[2] << ", "
      << physics.EnergyMomentum(parsym)[3] << ")" << std::endl;

  log << "-------Initial conditions------" << std::endl;

  for (int i = 0; i < parsym.no_of_particles; ++i) {
    log << "Particle: " << i << " ---------" << std::endl;
    log << "x, y = " << particle[i].x << ", " << particle[i].y << std::endl;
    log << "V = " << particle[i].vx << ", " << particle[i].vy << std::endl;
    log << "Omega = " << particle[i].omega << std::endl;
    log << "V0 = " << particle[i].vx_activity << ", " << particle[1].vy_activity
        << std::endl;
    log << "Omega0 = " << particle[i].omega_activity << std::endl;
  }
};

void state_after_simulation(std::ofstream &log, ParSim::ParticleSystem &parsym,
                            ParSim::Physics &physics) {
  log << "Energy-momentum After the collision: " << std::endl;
  log << "Total Energy: "
      << physics.EnergyMomentum(parsym)[0] + physics.EnergyMomentum(parsym)[1]
      << std::endl;
  log << "Translational K.Energy: " << physics.EnergyMomentum(parsym)[0]
      << std::endl;
  log << "Rotational K.Energy: " << physics.EnergyMomentum(parsym)[1]
      << std::endl;
  log << "Momentum: "
      << "(" << physics.EnergyMomentum(parsym)[2] << ", "
      << physics.EnergyMomentum(parsym)[3] << ")" << std::endl;
}
