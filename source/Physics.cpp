//
// Created by Vishu Saini.
//

#include "../include/Physics.h"
#include "../include/ParticleSystem.h"
#include <cmath>
#include <fstream>
#include <math.h>
#include <random>
#include <stdlib.h>
#include <vector>

// #######################################################//
//----------------------------------------
//
//  Class: Physics definitions
//
//--------------------------------------
// #######################################################//

// void ParSim::Physics::default_physics(){
//  move(elapsed);
// }

ParSim::Physics::Physics()
    : mt((std::random_device())()) { // initialize seed of mt using random
                                     // device

  parameters[0] = 0;
  parameters[1] = 0;
  parameters[2] = 0;
  parameters[3] = 0;
  parameters[4] = 0;
  parameters[5] = 0;
  parameters[6] = 0;
  parameters[7] = 0;
  parameters[8] = 0;
  parameters[9] = 0;
  parameters[10] = 0;
  parameters[11] = 0;
  parameters[12] = 0;

  WCA_parameters[0] = 0; // WCA epsilon
  WCA_parameters[1] = 0; // WCA sigma

  circular_wall_parameters[0] = 0.0; // WCA epsilon
  circular_wall_parameters[1] = 0.0; // WCA sigma

  NL_skin_depth = 0.0; // DeltaR
  NL_tau = 0;          // tau
}

/*Force linker + integrators-- */

// 1) Particle-Particle Force linkers ------------
void ParSim::Physics::Force_PP(ParSim::ParticleSystem &ps) {

  // Loop1: through all particles
  for (int i = 0; i < ps.no_of_particles; ++i) {

    // log << "1st loop -- " << std::endl;

    // Store previous step forces
    ps.particle_array[i].force_radial_prev[0] =
        ps.particle_array[i].force_radial[0];
    ps.particle_array[i].force_radial_prev[1] =
        ps.particle_array[i].force_radial[1];
    ps.particle_array[i].force_tangential_prev[0] =
        ps.particle_array[i].force_tangential[0];
    ps.particle_array[i].force_tangential_prev[1] =
        ps.particle_array[i].force_tangential[1];
    ps.particle_array[i].torque_prev = ps.particle_array[i].torque;

    // Reset the current forces
    ps.particle_array[i].force_radial[0] = 0;
    ps.particle_array[i].force_radial[1] = 0;
    ps.particle_array[i].force_tangential[0] = 0;
    ps.particle_array[i].force_tangential[1] = 0;
    ps.particle_array[i].torque = 0;

    // Unary force of damping. Always there. Translational and Rotational
    // activities added.
    ps.particle_array[i].force_radial[0] +=
        -1 * (this->parameters[5]) * ps.particle_array[i].vx +
        ps.particle_array[i].vx_activity;
    ps.particle_array[i].force_radial[1] +=
        -1 * (this->parameters[5]) * ps.particle_array[i].vy +
        ps.particle_array[i].vy_activity;
    ps.particle_array[i].torque +=
        -1 * (this->parameters[5]) * ps.particle_array[i].omega +
        ps.particle_array[i].omega_activity;

    // Knary force calculation --- Loop2: through all particles
    for (int j = 0; j < ps.no_of_particles; ++j) {
      if (j == i) { // no self coupling
        continue;
      }

      double d = ps.distance(ps.particle_array[i], ps.particle_array[j]);

      // U

      if (d <= (this->parameters[1])) {
        // mag. of radial force
        double N = (this->parameters[0]) *
                   ((this->parameters[1]) -
                    d); // magnitude of radial force used as normal reaction
        // radial interaction force
        ps.particle_array[i].force_radial[0] +=
            N * (ps.particle_array[i].x - ps.particle_array[j].x) /
            (d + (this->parameters[6]));

        ps.particle_array[i].force_radial[1] +=
            N * (ps.particle_array[i].y - ps.particle_array[j].y) /
            (d + (this->parameters[6]));

        // tangential friction force
        double omega_sum =
            (ps.particle_array[i].omega + ps.particle_array[j].omega);

        ps.particle_array[i].force_tangential[0] +=
            -(this->parameters[4]) *
            ((this->parameters[0]) * ((this->parameters[1]) - d)) *
            (omega_sum / (abs(omega_sum + (this->parameters[7])))) *
            ((ps.particle_array[i].y - ps.particle_array[j].y) /
             (d + (this->parameters[6])));

        ps.particle_array[i].force_tangential[1] +=
            -(this->parameters[4]) *
            ((this->parameters[0]) * ((this->parameters[1]) - d)) *
            (omega_sum / (abs(omega_sum + (this->parameters[7])))) *
            (-(ps.particle_array[i].x - ps.particle_array[j].x) /
             (d + (this->parameters[6])));

        // torque on particle

        if (omega_sum != 0) {
          ps.particle_array[i].torque +=
              -(this->parameters[4]) *
              ((this->parameters[0]) * ((this->parameters[1]) - d)) *
              (omega_sum / (abs(omega_sum + (this->parameters[7])))) * d;
        }

        if (omega_sum == 0) {
          ps.particle_array[i].torque +=
              -(this->parameters[4]) *
              ((this->parameters[0]) * ((this->parameters[1]) - d)) *
              (omega_sum / (abs(omega_sum + (this->parameters[7])))) * d;
        }
      }
    }
  }
}

void ParSim::Physics::Force_PP2(ParSim::ParticleSystem &ps) {

  // Loop1: through all particles
  for (int i = 0; i < ps.no_of_particles; ++i) {

    // log << "1st loop -- " << std::endl;

    // Store previous step forces
    ps.particle_array[i].force_radial_prev[0] =
        ps.particle_array[i].force_radial[0];
    ps.particle_array[i].force_radial_prev[1] =
        ps.particle_array[i].force_radial[1];
    ps.particle_array[i].force_tangential_prev[0] =
        ps.particle_array[i].force_tangential[0];
    ps.particle_array[i].force_tangential_prev[1] =
        ps.particle_array[i].force_tangential[1];
    ps.particle_array[i].torque_prev = ps.particle_array[i].torque;

    // Reset the current forces
    ps.particle_array[i].force_radial[0] = 0;
    ps.particle_array[i].force_radial[1] = 0;
    ps.particle_array[i].force_tangential[0] = 0;
    ps.particle_array[i].force_tangential[1] = 0;
    ps.particle_array[i].torque = 0;

    // Unary force of damping. Always there. Translational and Rotational
    // activities added.
    ps.particle_array[i].force_radial[0] +=
        -1 * (this->parameters[5]) * ps.particle_array[i].vx +
        ps.particle_array[i].vx_activity;
    ps.particle_array[i].force_radial[1] +=
        -1 * (this->parameters[5]) * ps.particle_array[i].vy +
        ps.particle_array[i].vy_activity;
    ps.particle_array[i].torque +=
        -1 * (this->parameters[5]) * ps.particle_array[i].omega +
        ps.particle_array[i].omega_activity;
  }
  // Loop2: N(N-1)/2 loop for pair force calculation
  for (int i = 0; i < (ps.no_of_particles - 1); ++i) {

    for (int j = i + 1; j < ps.no_of_particles; ++j) {

      double d = ps.distance(ps.particle_array[i], ps.particle_array[j]);

      // U

      if (d <= (this->parameters[1])) {
        // mag. of radial force
        double N = (this->parameters[0]) *
                   ((this->parameters[1]) -
                    d); // magnitude of radial force used as normal reaction
        // radial interaction force
        double frx = N * (ps.particle_array[i].x - ps.particle_array[j].x) /
                     (d + (this->parameters[6]));

        double fry = N * (ps.particle_array[i].y - ps.particle_array[j].y) /
                     (d + (this->parameters[6]));

        ps.particle_array[i].force_radial[0] += frx;
        ps.particle_array[i].force_radial[1] += fry;

        ps.particle_array[j].force_radial[0] += -frx;
        ps.particle_array[j].force_radial[1] += -fry;

        // tangential friction force
        double omega_sum =
            (ps.particle_array[i].omega + ps.particle_array[j].omega);

        double ftx = -(this->parameters[4]) *
                     ((this->parameters[0]) * ((this->parameters[1]) - d)) *
                     (omega_sum / (abs(omega_sum + (this->parameters[7])))) *
                     ((ps.particle_array[i].y - ps.particle_array[j].y) /
                      (d + (this->parameters[6])));

        double fty = -(this->parameters[4]) *
                     ((this->parameters[0]) * ((this->parameters[1]) - d)) *
                     (omega_sum / (abs(omega_sum + (this->parameters[7])))) *
                     (-(ps.particle_array[i].x - ps.particle_array[j].x) /
                      (d + (this->parameters[6])));

        ps.particle_array[i].force_tangential[0] += ftx;
        ps.particle_array[i].force_tangential[1] += fty;

        ps.particle_array[j].force_tangential[0] += -ftx;
        ps.particle_array[j].force_tangential[1] += -fty;

        // torque on particle

        if (omega_sum != 0) {
          ps.particle_array[i].torque +=
              -(this->parameters[4]) *
              ((this->parameters[0]) * ((this->parameters[1]) - d)) *
              (omega_sum / (abs(omega_sum + (this->parameters[7])))) * d;

          ps.particle_array[j].torque +=
              -(this->parameters[4]) *
              ((this->parameters[0]) * ((this->parameters[1]) - d)) *
              (omega_sum / (abs(omega_sum + (this->parameters[7])))) * d;
        }

        if (omega_sum == 0) {
          ps.particle_array[i].torque +=
              -(this->parameters[4]) *
              ((this->parameters[0]) * ((this->parameters[1]) - d)) *
              (omega_sum / (abs(omega_sum + (this->parameters[7])))) * d;
          ps.particle_array[j].torque +=
              -(this->parameters[4]) *
              ((this->parameters[0]) * ((this->parameters[1]) - d)) *
              (omega_sum / (abs(omega_sum + (this->parameters[7])))) * d;
        }
      }
    }
  }
}

void ParSim::Physics::Force_PP_PBC(ParSim::ParticleSystem &ps) {

  Particle *const particle = ps.get_particles();

  // Loop1: through all particles
  for (int i = 0; i < ps.no_of_particles; ++i) {

    // Store previous step forces
    particle[i].force_radial_prev[0] = particle[i].force_radial[0];
    particle[i].force_radial_prev[1] = particle[i].force_radial[1];
    particle[i].force_tangential_prev[0] = particle[i].force_tangential[0];
    particle[i].force_tangential_prev[1] = particle[i].force_tangential[1];
    particle[i].torque_prev = particle[i].torque;

    // Reset the current forces and pressures
    particle[i].force_radial[0] = 0;
    particle[i].force_radial[1] = 0;
    particle[i].force_tangential[0] = 0;
    particle[i].force_tangential[1] = 0;
    particle[i].torque = 0;

    particle[i].sigma[0] = 0.0;
    particle[i].sigma[1] = 0.0;
    particle[i].sigma[2] = 0.0;
    particle[i].sigma[3] = 0.0;

    // wall forces and torque
    // double force_wall_y = 0.0;
    // double force_wall_fric_x = 0.0;
    // double torque_wall = 0.0;
    // double yp = 0.0;

    // if (ps.particle_array[i].y > 0) {
    //   yp = ps.particle_array[i].y + (this->parameters[1] / 2);

    // } else {
    //   yp = ps.particle_array[i].y - (this->parameters[1] / 2);
    // }

    // if (abs(yp) > ps.L / 2) {
    //   double norm = 2 * (abs(yp) - (ps.L / 2)); // kwall = 2
    //   force_wall_y =
    //       -norm * ps.particle_array[i].y / (abs(ps.particle_array[i].y));

    //   force_wall_fric_x =
    //       -(this->parameters[4] / 2.0) * abs(norm) *
    //       (ps.particle_array[i].omega /
    //        (abs(ps.particle_array[i].omega + (this->parameters[7])))) *
    //       (-(ps.particle_array[i].y) / abs(ps.particle_array[i].y));

    //   if (ps.particle_array[i].y > 0) {
    //     torque_wall =
    //         -(this->parameters[4] / 2.0) * abs(norm) *
    //         (ps.particle_array[i].omega /
    //          (abs(ps.particle_array[i].omega + (this->parameters[7])))) *
    //         (abs((ps.L / 2) - (ps.particle_array[i].y)));

    //   } else {
    //     torque_wall =
    //         -(this->parameters[4] / 2.0) * abs(norm) *
    //         (ps.particle_array[i].omega /
    //          (abs(ps.particle_array[i].omega + (this->parameters[7])))) *
    //         (abs((ps.L / 2) + (ps.particle_array[i].y)));
    //   }

    // } else {
    //   force_wall_y = 0.0;
    //   force_wall_fric_x = 0.0;
    //   torque_wall = 0.0;
    // }

    // Unary force of damping. Always there. Translational and Rotational
    // activities added (note this is not radial !!!).

    particle[i].force_radial[0] +=
        0; //-1 * (this->parameters[5]) * particle[i].vx;

    particle[i].force_radial[1] +=
        0; //-1 * (this->parameters[5]) * particle[i].vy;

    particle[i].torque += -1 * (this->parameters[10]) * particle[i].omega +
                          (particle[i].omega_activity * (this->parameters[10]));

    // Knary force calculation --- Loop2: through all particles
    for (int j = 0; j < ps.no_of_particles; ++j) {
      if (j == i) { // no self coupling
        continue;
      }

      // Distance from the nearest image of jth particle

      // i)calculate min separation dx and dy
      double dx = particle[i].x - particle[j].x;

      if (dx < (-ps.L / 2)) {
        dx += ps.L;
      } else if (dx > (ps.L / 2)) {
        dx -= ps.L;
      }

      double dy = particle[i].y - particle[j].y;

      if (dy < (-ps.L / 2)) {
        dy += ps.L;
      } else if (dy > (ps.L / 2)) {
        dy -= ps.L;
      }
      // ii) Nearest image distance
      double d = std::sqrt(dx * dx + dy * dy);

      // U

      if (d <= (this->parameters[1])) {

        // relative velocities calculation (ith particle at rest)
        double vrel_x = particle[j].vx - particle[i].vx;
        double vrel_y = particle[j].vy - particle[i].vy;

        // mag. of radial force
        double N = (this->parameters[0]) *
                   ((this->parameters[1]) -
                    d); // magnitude of radial force used as normal reaction

        // components of unit separation vector r

        double r1 = dx / (d + (this->parameters[6]));
        double r2 = dy / (d + (this->parameters[6]));

        double vreldotr = (r1 * vrel_x + r2 * vrel_y);
        // Radial interaction force with jth particle
        double frx = N * r1 + (this->parameters[5]) * vreldotr * r1;
        double fry = N * r2 + (this->parameters[5]) * vreldotr * r2;

        particle[i].force_radial[0] += frx;
        particle[i].force_radial[1] += fry;
        // Tangential friction force with particle j
        double omega_sum = (particle[i].omega + particle[j].omega);

        double ftx = -(this->parameters[4]) *
                         ((this->parameters[0]) * ((this->parameters[1]) - d)) *
                         (omega_sum / 2.0) * r2 +
                     (this->parameters[5]) * (vrel_x - vreldotr * r1);

        double fty = -(this->parameters[4]) *
                         ((this->parameters[0]) * ((this->parameters[1]) - d)) *
                         (omega_sum / 2.0) * (-r1) +
                     (this->parameters[5]) * (vrel_y - vreldotr * r2);

        particle[i].force_tangential[0] += ftx;
        particle[i].force_tangential[1] += fty;

        // torque on particle (generic torque, should be written as cross
        // product instead )

        if (omega_sum != 0) {
          particle[i].torque +=
              -(this->parameters[4]) *
              ((this->parameters[0]) * ((this->parameters[1]) - d)) *
              (omega_sum / 2.0) * d;
        }

        if (omega_sum == 0) {
          particle[i].torque +=
              -(this->parameters[4]) *
              ((this->parameters[0]) * ((this->parameters[1]) - d)) *
              (omega_sum / 2.0) * d;
        }

        // pressure term calculation

        particle[i].sigma[0] += (frx + ftx) * (particle[i].x - particle[j].x);
        particle[i].sigma[1] += (frx + ftx) * (particle[i].y - particle[j].y);
        particle[i].sigma[2] += (fry + fty) * (particle[i].x - particle[j].x);
        particle[i].sigma[3] += (fry + fty) * (particle[i].y - particle[j].y);
      }
    }
  }
}

void ParSim::Physics::Force_PP_PBC2(ParSim::ParticleSystem &ps) {

  Particle *const particle = ps.get_particles();

  // Loop1: through all particles
  for (int i = 0; i < ps.no_of_particles; ++i) {

    // Store previous step forces
    particle[i].force_radial_prev[0] = particle[i].force_radial[0];
    particle[i].force_radial_prev[1] = particle[i].force_radial[1];
    particle[i].force_tangential_prev[0] = particle[i].force_tangential[0];
    particle[i].force_tangential_prev[1] = particle[i].force_tangential[1];
    particle[i].torque_prev = particle[i].torque;

    // Reset the current forces and pressures
    particle[i].force_radial[0] = 0;
    particle[i].force_radial[1] = 0;
    particle[i].force_tangential[0] = 0;
    particle[i].force_tangential[1] = 0;
    particle[i].torque = 0;

    particle[i].sigma[0] = 0.0;
    particle[i].sigma[1] = 0.0;
    particle[i].sigma[2] = 0.0;
    particle[i].sigma[3] = 0.0;

    // Unary force of damping. Always there. Translational and Rotational
    // activities added (note this is not radial !!!).

    particle[i].force_radial[0] +=
        0; //-1 * (this->parameters[5]) * particle[i].vx;

    particle[i].force_radial[1] +=
        0; //-1 * (this->parameters[5]) * particle[i].vy;

    particle[i].torque += -1 * (this->parameters[10]) * particle[i].omega +
                          (particle[i].omega_activity * (this->parameters[10]));
  }
  // Loop2: N(N-1)/2 loop for pair force calculation
  for (int i = 0; i < (ps.no_of_particles - 1); ++i) {

    for (int j = i + 1; j < ps.no_of_particles; ++j) {

      // Distance from the nearest image of jth particle

      // i)calculate min separation dx and dy
      double dx = particle[i].x - particle[j].x;

      if (dx < (-ps.L / 2)) {
        dx += ps.L;
      } else if (dx > (ps.L / 2)) {
        dx -= ps.L;
      }

      double dy = particle[i].y - particle[j].y;

      if (dy < (-ps.L / 2)) {
        dy += ps.L;
      } else if (dy > (ps.L / 2)) {
        dy -= ps.L;
      }
      // ii) Nearest image distance
      double d = std::sqrt(dx * dx + dy * dy);

      // U

      if (d <= (this->parameters[1])) {

        // relative velocities calculation (ith particle at rest)
        double vrel_x = particle[j].vx - particle[i].vx;
        double vrel_y = particle[j].vy - particle[i].vy;

        // mag. of radial force
        double N = (this->parameters[0]) *
                   ((this->parameters[1]) -
                    d); // magnitude of radial force used as normal reaction

        // components of unit separation vector r

        double r1 = dx / (d + (this->parameters[6]));
        double r2 = dy / (d + (this->parameters[6]));

        double vreldotr = (r1 * vrel_x + r2 * vrel_y);
        vreldotr = 0;
        // Radial interaction force with jth particle
        double frx = N * r1 + (this->parameters[5]) * vreldotr * r1;
        double fry = N * r2 + (this->parameters[5]) * vreldotr * r2;

        particle[i].force_radial[0] += frx;
        particle[i].force_radial[1] += fry;

        particle[j].force_radial[0] += -frx;
        particle[j].force_radial[1] += -fry;
        // Tangential friction force with particle j
        double omega_sum = (particle[i].omega + particle[j].omega);

        double ftx = -(this->parameters[4]) *
                         ((this->parameters[0]) * ((this->parameters[1]) - d)) *
                         (omega_sum / 2.0) * r2 +
                     (this->parameters[5]) * (vrel_x - vreldotr * r1);

        double fty = -(this->parameters[4]) *
                         ((this->parameters[0]) * ((this->parameters[1]) - d)) *
                         (omega_sum / 2.0) * (-r1) +
                     (this->parameters[5]) * (vrel_y - vreldotr * r2);

        particle[i].force_tangential[0] += ftx;
        particle[i].force_tangential[1] += fty;

        particle[j].force_tangential[0] += -ftx;
        particle[j].force_tangential[1] += -fty;

        // torque on particle (generic torque, should be written as cross
        // product instead )

        if (omega_sum != 0) {
          particle[i].torque +=
              -(this->parameters[4]) *
              ((this->parameters[0]) * ((this->parameters[1]) - d)) *
              (omega_sum / 2.0) * d;
          particle[j].torque +=
              -(this->parameters[4]) *
              ((this->parameters[0]) * ((this->parameters[1]) - d)) *
              (omega_sum / 2.0) * d;
        }

        if (omega_sum == 0) {
          particle[i].torque +=
              -(this->parameters[4]) *
              ((this->parameters[0]) * ((this->parameters[1]) - d)) *
              (omega_sum / 2.0) * d;
          particle[j].torque +=
              -(this->parameters[4]) *
              ((this->parameters[0]) * ((this->parameters[1]) - d)) *
              (omega_sum / 2.0) * d;
        }

        // pressure term calculation

        particle[i].sigma[0] += (frx + ftx) * (particle[i].x - particle[j].x);
        particle[i].sigma[1] += (frx + ftx) * (particle[i].y - particle[j].y);
        particle[i].sigma[2] += (fry + fty) * (particle[i].x - particle[j].x);
        particle[i].sigma[3] += (fry + fty) * (particle[i].y - particle[j].y);

        particle[j].sigma[0] += -(frx + ftx) * (particle[i].x - particle[j].x);
        particle[j].sigma[1] += -(frx + ftx) * (particle[i].y - particle[j].y);
        particle[j].sigma[2] += -(fry + fty) * (particle[i].x - particle[j].x);
        particle[j].sigma[3] += -(fry + fty) * (particle[i].y - particle[j].y);
      }
    }
  }
}

void ParSim::Physics::Force_PP_CRB(ParSim::ParticleSystem &ps) {

  Particle *const particle = ps.get_particles();

  // Loop1: through all particles
  for (int i = 0; i < ps.no_of_particles; ++i) {

    // Store previous step forces
    particle[i].force_radial_prev[0] = particle[i].force_radial[0];
    particle[i].force_radial_prev[1] = particle[i].force_radial[1];
    particle[i].force_tangential_prev[0] = particle[i].force_tangential[0];
    particle[i].force_tangential_prev[1] = particle[i].force_tangential[1];
    particle[i].torque_prev = particle[i].torque;

    // Reset the current forces and pressures
    particle[i].force_radial[0] = 0;
    particle[i].force_radial[1] = 0;
    particle[i].force_tangential[0] = 0;
    particle[i].force_tangential[1] = 0;
    particle[i].torque = 0;

    particle[i].sigma[0] = 0.0;
    particle[i].sigma[1] = 0.0;
    particle[i].sigma[2] = 0.0;
    particle[i].sigma[3] = 0.0;

    // wall forces -----------------
    double force_wall_x = 0.0;
    double force_wall_y = 0.0;

    double force_wall_fric_x = 0.0;
    double force_wall_fric_y = 0.0;

    double torque_wall = 0.0;

    double rp = 0.0;

    double magnitude = 0.0;

    rp = std::sqrt(particle[i].x * particle[i].x +
                   particle[i].y * particle[i].y) +
         (this->parameters[1] /
          2); // calculates position vector magnitude + (sigma/2)

    double radius = ps.L;

    if (rp > radius) {

      double flag = ((this->parameters[1]) -
                     abs(rp - radius)) / // proxy for (r/sigma) in WCA
                    this->circular_wall_parameters[1];
      double flag5 = flag * flag * flag * flag * flag;

      double norm =
          4 *
          (this->circular_wall_parameters[0] /
           this->circular_wall_parameters[1]) *
          (1.0 /
           (flag5 * flag5 * flag5 * flag5 * flag5 + (this->parameters[6]))) *
          ((48.0 / (flag5 * flag5 * flag5 * flag5 * flag * flag * flag * flag +
                    (this->parameters[6]))) -
           (24.0)); // WCA force

      force_wall_x = -norm * particle[i].x / (rp - this->parameters[1]);

      force_wall_y = -norm * particle[i].y / (rp - this->parameters[1]);

      force_wall_fric_x = -(this->parameters[4] / 2.0) * abs(norm) *
                          (particle[i].omega /
                           (abs(particle[i].omega + (this->parameters[7])))) *
                          (-(particle[i].y) / (rp - this->parameters[1]));

      force_wall_fric_y = -(this->parameters[4] / 2.0) * abs(norm) *
                          (particle[i].omega /
                           (abs(particle[i].omega + (this->parameters[7])))) *
                          ((particle[i].x) / (rp - this->parameters[1]));

    } else {
      force_wall_x = 0.0;
      force_wall_fric_x = 0.0;

      force_wall_y = 0.0;
      force_wall_fric_y = 0.0;
    }

    // wall torque

    // if (rp > radius) {
    //   torque_wall =
    //       -(this->parameters[4]/2.0) * abs(norm) *
    //       (particle[i].omega /
    //        (abs(particle[i].omega + (this->parameters[7])))) *
    //       ((ps.L / 2) - (rp - this->parameters[1]));

    // } else {
    //   torque_wall = 0.0;
    // }

    // wall friction reduced to null

    torque_wall = 0.0;
    force_wall_fric_x = 0.0;
    force_wall_fric_y = 0.0;

    // Unary force of damping. Always there. Translational and Rotational
    // activities added. (note these are not purely radial and tangential)

    particle[i].force_radial[0] += -1 * (this->parameters[5]) * particle[i].vx +
                                   force_wall_x + force_wall_fric_x;

    particle[i].force_radial[1] += -1 * (this->parameters[5]) * particle[i].vy +
                                   force_wall_y + force_wall_fric_y;

    particle[i].torque +=
        -1 * (this->parameters[10]) * particle[i].omega +
        (particle[i].omega_activity * (this->parameters[10])) + torque_wall;

    // Knary force calculation --- Loop2: through all particles
    for (int j = 0; j < ps.no_of_particles; ++j) {
      if (j == i) { // no self coupling
        continue;
      }

      // double d = ps.distance(par[i], par[j]);

      double d = sqrt(
          (particle[i].x - particle[j].x) * (particle[i].x - particle[j].x) +
          (particle[i].y - particle[j].y) * (particle[i].y - particle[j].y));

      // U

      if (d <= (this->parameters[1])) {

        // relative velocities calculation (ith particle at rest)
        double vrel_x = 0; // particle[j].vx - particle[i].vx;
        double vrel_y = 0; // particle[j].vy - particle[i].vy;

        // mag. of radial force
        double N = (this->parameters[0]) *
                   ((this->parameters[1]) -
                    d); // magnitude of radial force used as normal reaction

        // components of unit separation vector r
        double r1 =
            (particle[i].x - particle[j].x) / (d + (this->parameters[6])); // rx
        double r2 =
            (particle[i].y - particle[j].y) / (d + (this->parameters[6])); // ry

        double vreldotr = (r1 * vrel_x + r2 * vrel_y);
        // Radial interaction force with jth particle
        double frx = N * r1 + (this->parameters[5]) * vreldotr * r1;
        double fry = N * r2 + (this->parameters[5]) * vreldotr * r2;

        particle[i].force_radial[0] += frx;

        particle[i].force_radial[1] += fry;

        // Tangential friction force with particle j
        double omega_sum = (particle[i].omega + particle[j].omega);

        double ftx = -(this->parameters[4]) *
                         ((this->parameters[0]) * ((this->parameters[1]) - d)) *
                         (omega_sum / (2.0 + (this->parameters[7]))) * r2 +
                     (this->parameters[5]) * (vrel_x - vreldotr * r1);

        double fty = -(this->parameters[4]) *
                         ((this->parameters[0]) * ((this->parameters[1]) - d)) *
                         (omega_sum / (2.0 + (this->parameters[7]))) * (-r1) +
                     (this->parameters[5]) * (vrel_y - vreldotr * r2);

        particle[i].force_tangential[0] += ftx;

        particle[i].force_tangential[1] += fty;

        // torque on particle

        if (omega_sum != 0) {
          particle[i].torque +=
              -(this->parameters[4]) *
              ((this->parameters[0]) * ((this->parameters[1]) - d)) *
              (omega_sum / (2.0 + (this->parameters[7]))) * (d / 2);
        }

        if (omega_sum == 0) {
          particle[i].torque +=
              -(this->parameters[4]) *
              ((this->parameters[0]) * ((this->parameters[1]) - d)) *
              (omega_sum / (2.0 + (this->parameters[7]))) * (d / 2);
        }

        // pressure term calculation

        particle[i].sigma[0] += (frx + ftx) * (particle[i].x - particle[j].x);
        particle[i].sigma[1] += (frx + ftx) * (particle[i].y - particle[j].y);
        particle[i].sigma[2] += (fry + fty) * (particle[i].x - particle[j].x);
        particle[i].sigma[3] += (fry + fty) * (particle[i].y - particle[j].y);
      }
    }
  }
}

void ParSim::Physics::Force_PP_CRB2(ParSim::ParticleSystem &ps) {

  Particle *const particle = ps.get_particles();

  // Loop1: through all particles
  for (int i = 0; i < ps.no_of_particles; ++i) {

    // Store previous step forces
    particle[i].force_radial_prev[0] = particle[i].force_radial[0];
    particle[i].force_radial_prev[1] = particle[i].force_radial[1];
    particle[i].force_tangential_prev[0] = particle[i].force_tangential[0];
    particle[i].force_tangential_prev[1] = particle[i].force_tangential[1];
    particle[i].torque_prev = particle[i].torque;

    // Reset the current forces and pressures
    particle[i].force_radial[0] = 0;
    particle[i].force_radial[1] = 0;
    particle[i].force_tangential[0] = 0;
    particle[i].force_tangential[1] = 0;
    particle[i].torque = 0;

    particle[i].sigma[0] = 0.0;
    particle[i].sigma[1] = 0.0;
    particle[i].sigma[2] = 0.0;
    particle[i].sigma[3] = 0.0;

    // wall forces -----------------
    double force_wall_x = 0.0;
    double force_wall_y = 0.0;

    double force_wall_fric_x = 0.0;
    double force_wall_fric_y = 0.0;

    double torque_wall = 0.0;

    double rp = 0.0;

    double magnitude = 0.0;

    rp = std::sqrt(particle[i].x * particle[i].x +
                   particle[i].y * particle[i].y) +
         (this->parameters[1] /
          2); // calculates position vector magnitude + (sigma/2)

    double radius = ps.L;

    if (rp > radius) {

      double flag = ((this->parameters[1]) -
                     abs(rp - radius)) / // proxy for (r/sigma) in WCA
                    this->circular_wall_parameters[1];
      double flag5 = flag * flag * flag * flag * flag;

      double norm =
          4 *
          (this->circular_wall_parameters[0] /
           this->circular_wall_parameters[1]) *
          (1.0 /
           (flag5 * flag5 * flag5 * flag5 * flag5 + (this->parameters[6]))) *
          ((48.0 / (flag5 * flag5 * flag5 * flag5 * flag * flag * flag * flag +
                    (this->parameters[6]))) -
           (24.0)); // WCA force

      force_wall_x = -norm * particle[i].x / (rp - this->parameters[1]);

      force_wall_y = -norm * particle[i].y / (rp - this->parameters[1]);

      force_wall_fric_x = -(this->parameters[4] / 2.0) * abs(norm) *
                          (particle[i].omega /
                           (abs(particle[i].omega + (this->parameters[7])))) *
                          (-(particle[i].y) / (rp - this->parameters[1]));

      force_wall_fric_y = -(this->parameters[4] / 2.0) * abs(norm) *
                          (particle[i].omega /
                           (abs(particle[i].omega + (this->parameters[7])))) *
                          ((particle[i].x) / (rp - this->parameters[1]));

    } else {
      force_wall_x = 0.0;
      force_wall_fric_x = 0.0;

      force_wall_y = 0.0;
      force_wall_fric_y = 0.0;
    }

    // wall torque

    // if (rp > radius) {
    //   torque_wall =
    //       -(this->parameters[4]/2.0) * abs(norm) *
    //       (particle[i].omega /
    //        (abs(particle[i].omega + (this->parameters[7])))) *
    //       ((ps.L / 2) - (rp - this->parameters[1]));

    // } else {
    //   torque_wall = 0.0;
    // }

    // wall friction reduced to null

    torque_wall = 0.0;
    force_wall_fric_x = 0.0;
    force_wall_fric_y = 0.0;

    // Unary force of damping. Always there. Translational and Rotational
    // activities added. (note these are not purely radial and tangential)

    particle[i].force_radial[0] += -1 * (this->parameters[5]) * particle[i].vx +
                                   force_wall_x + force_wall_fric_x;

    particle[i].force_radial[1] += -1 * (this->parameters[5]) * particle[i].vy +
                                   force_wall_y + force_wall_fric_y;

    particle[i].torque +=
        -1 * (this->parameters[10]) * particle[i].omega +
        (particle[i].omega_activity * (this->parameters[10])) + torque_wall;
  }
  // Loop 2: N(N-1)/2 Loop for pair force calculation ------
  for (int i = 0; i < (ps.no_of_particles - 1); ++i) {

    for (int j = i + 1; j < ps.no_of_particles; ++j) {

      // double d = ps.distance(par[i], par[j]);

      double d = sqrt(
          (particle[i].x - particle[j].x) * (particle[i].x - particle[j].x) +
          (particle[i].y - particle[j].y) * (particle[i].y - particle[j].y));

      // U

      if (d <= (this->parameters[1])) {

        // relative velocities calculation (ith particle at rest)
        double vrel_x = 0; // particle[j].vx - particle[i].vx;
        double vrel_y = 0; // particle[j].vy - particle[i].vy;

        // mag. of radial force
        double N = (this->parameters[0]) *
                   ((this->parameters[1]) -
                    d); // magnitude of radial force used as normal reaction

        // components of unit separation vector r
        double r1 =
            (particle[i].x - particle[j].x) / (d + (this->parameters[6])); // rx
        double r2 =
            (particle[i].y - particle[j].y) / (d + (this->parameters[6])); // ry

        double vreldotr = (r1 * vrel_x + r2 * vrel_y);
        // Radial interaction force with jth particle
        double frx = N * r1 + (this->parameters[5]) * vreldotr * r1;
        double fry = N * r2 + (this->parameters[5]) * vreldotr * r2;

        particle[i].force_radial[0] += frx;
        particle[i].force_radial[1] += fry;

        particle[j].force_radial[0] += -frx;
        particle[j].force_radial[1] += -fry;

        // Tangential friction force with particle j
        double omega_sum = (particle[i].omega + particle[j].omega);

        double ftx = -(this->parameters[4]) *
                         ((this->parameters[0]) * ((this->parameters[1]) - d)) *
                         (omega_sum / (2.0 + (this->parameters[7]))) * r2 +
                     (this->parameters[5]) * (vrel_x - vreldotr * r1);

        double fty = -(this->parameters[4]) *
                         ((this->parameters[0]) * ((this->parameters[1]) - d)) *
                         (omega_sum / (2.0 + (this->parameters[7]))) * (-r1) +
                     (this->parameters[5]) * (vrel_y - vreldotr * r2);

        particle[i].force_tangential[0] += ftx;
        particle[i].force_tangential[1] += fty;

        particle[j].force_tangential[0] += -ftx;
        particle[j].force_tangential[1] += -fty;

        // torque on particle

        if (omega_sum != 0) {
          particle[i].torque +=
              -(this->parameters[4]) *
              ((this->parameters[0]) * ((this->parameters[1]) - d)) *
              (omega_sum / (2.0 + (this->parameters[7]))) * (d / 2);
          particle[j].torque +=
              -(this->parameters[4]) *
              ((this->parameters[0]) * ((this->parameters[1]) - d)) *
              (omega_sum / (2.0 + (this->parameters[7]))) * (d / 2);
        }

        if (omega_sum == 0) {
          particle[i].torque +=
              -(this->parameters[4]) *
              ((this->parameters[0]) * ((this->parameters[1]) - d)) *
              (omega_sum / (2.0 + (this->parameters[7]))) * (d / 2);
          particle[j].torque +=
              -(this->parameters[4]) *
              ((this->parameters[0]) * ((this->parameters[1]) - d)) *
              (omega_sum / (2.0 + (this->parameters[7]))) * (d / 2);
        }

        // pressure term calculation

        particle[i].sigma[0] += (frx + ftx) * (particle[i].x - particle[j].x);
        particle[i].sigma[1] += (frx + ftx) * (particle[i].y - particle[j].y);
        particle[i].sigma[2] += (fry + fty) * (particle[i].x - particle[j].x);
        particle[i].sigma[3] += (fry + fty) * (particle[i].y - particle[j].y);

        particle[j].sigma[0] += -(frx + ftx) * (particle[i].x - particle[j].x);
        particle[j].sigma[1] += -(frx + ftx) * (particle[i].y - particle[j].y);
        particle[j].sigma[2] += -(fry + fty) * (particle[i].x - particle[j].x);
        particle[j].sigma[3] += -(fry + fty) * (particle[i].y - particle[j].y);
      }
    }
  }
}
void ParSim::Physics::Force_PP_CRB_WCA(ParSim::ParticleSystem &ps) {

  Particle *const particle = ps.get_particles();

  // Loop1: through all particles
  for (int i = 0; i < ps.no_of_particles; ++i) {

    // Store previous step forces
    particle[i].force_radial_prev[0] = particle[i].force_radial[0];
    particle[i].force_radial_prev[1] = particle[i].force_radial[1];
    particle[i].force_tangential_prev[0] = particle[i].force_tangential[0];
    particle[i].force_tangential_prev[1] = particle[i].force_tangential[1];
    particle[i].torque_prev = particle[i].torque;

    // Reset the current forces and pressures
    particle[i].force_radial[0] = 0;
    particle[i].force_radial[1] = 0;
    particle[i].force_tangential[0] = 0;
    particle[i].force_tangential[1] = 0;
    particle[i].torque = 0;

    particle[i].sigma[0] = 0.0;
    particle[i].sigma[1] = 0.0;
    particle[i].sigma[2] = 0.0;
    particle[i].sigma[3] = 0.0;

    // wall forces -----------------
    double force_wall_x = 0.0;
    double force_wall_y = 0.0;

    double force_wall_fric_x = 0.0;
    double force_wall_fric_y = 0.0;

    double torque_wall = 0.0;

    double rp = 0.0;

    double magnitude = 0.0;

    rp = std::sqrt(particle[i].x * particle[i].x +
                   particle[i].y * particle[i].y) +
         (this->parameters[1] /
          2); // calculates position vector magnitude + (sigma/2)

    double radius = ps.L;

    if (rp > radius) {

      double flag = ((this->parameters[1]) -
                     abs(rp - radius)) / // proxy for (r/sigma) in WCA
                    this->circular_wall_parameters[1];
      double flag5 = flag * flag * flag * flag * flag;

      double norm =
          4 *
          (this->circular_wall_parameters[0] /
           this->circular_wall_parameters[1]) *
          (1.0 / (flag5 * flag * flag + (this->parameters[6]))) *
          ((12.0 / (flag5 + (this->parameters[6]))) - (6.0)); // WCA force

      force_wall_x = -norm * particle[i].x / (rp - this->parameters[1]);

      force_wall_y = -norm * particle[i].y / (rp - this->parameters[1]);

      force_wall_fric_x = -(this->parameters[4] / 2.0) * abs(norm) *
                          (particle[i].omega /
                           (abs(particle[i].omega + (this->parameters[7])))) *
                          (-(particle[i].y) / (rp - this->parameters[1]));

      force_wall_fric_y = -(this->parameters[4] / 2.0) * abs(norm) *
                          (particle[i].omega /
                           (abs(particle[i].omega + (this->parameters[7])))) *
                          ((particle[i].x) / (rp - this->parameters[1]));

    } else {
      force_wall_x = 0.0;
      force_wall_fric_x = 0.0;

      force_wall_y = 0.0;
      force_wall_fric_y = 0.0;
    }

    // wall torque

    // if (rp > radius) {
    //   torque_wall =
    //       -(this->parameters[4]/2.0) * abs(norm) *
    //       (particle[i].omega /
    //        (abs(particle[i].omega + (this->parameters[7])))) *
    //       ((ps.L / 2) - (rp - this->parameters[1]));

    // } else {
    //   torque_wall = 0.0;
    // }

    // wall friction reduced to null

    torque_wall = 0.0;
    force_wall_fric_x = 0.0;
    force_wall_fric_y = 0.0;

    // Unary force of damping. Always there. Translational and Rotational
    // activities added. (note these are not purely radial and tangential)

    particle[i].force_radial[0] += -1 * (this->parameters[5]) * particle[i].vx +
                                   force_wall_x + force_wall_fric_x;

    particle[i].force_radial[1] += -1 * (this->parameters[5]) * particle[i].vy +
                                   force_wall_y + force_wall_fric_y;

    particle[i].torque +=
        -1 * (this->parameters[10]) * particle[i].omega +
        (particle[i].omega_activity * (this->parameters[10])) + torque_wall;

    // Knary force calculation --- Loop2: through all particles
    for (int j = 0; j < ps.no_of_particles; ++j) {
      if (j == i) { // no self coupling
        continue;
      }

      // double d = ps.distance(par[i], par[j]);

      double d = sqrt(
          (particle[i].x - particle[j].x) * (particle[i].x - particle[j].x) +
          (particle[i].y - particle[j].y) * (particle[i].y - particle[j].y));

      // U

      if (d <= 1.12 * (this->parameters[1])) {

        // relative velocities calculation (ith particle at rest)
        double vrel_x = 0; // particle[j].vx - particle[i].vx;
        double vrel_y = 0; // particle[j].vy - particle[i].vy;

        // mag. of radial force  (positive)
        double r_by_sigma = d / this->WCA_parameters[1]; // r/sigma in WCA
        double r_by_sigma5 =
            r_by_sigma * r_by_sigma * r_by_sigma * r_by_sigma * r_by_sigma;

        double N =
            4 * (this->WCA_parameters[0] / this->WCA_parameters[1]) *
            (1.0 /
             (r_by_sigma5 * r_by_sigma5 * r_by_sigma * r_by_sigma * r_by_sigma +
              (this->parameters[6]))) *
            ((24.0 / (r_by_sigma5 * r_by_sigma5 * r_by_sigma * r_by_sigma +
                      (this->parameters[6]))) -
             (12.0)); // WCA force magnitude (positive)

        // components of unit separation vector r
        double r1 =
            (particle[i].x - particle[j].x) / (d + (this->parameters[6])); // rx
        double r2 =
            (particle[i].y - particle[j].y) / (d + (this->parameters[6])); // ry

        double vreldotr = (r1 * vrel_x + r2 * vrel_y);
        // Radial interaction force with jth particle
        double frx = N * r1 + (this->parameters[5]) * vreldotr * r1;
        double fry = N * r2 + (this->parameters[5]) * vreldotr * r2;

        particle[i].force_radial[0] += frx;

        particle[i].force_radial[1] += fry;

        // Tangential friction force with particle j
        double omega_sum = (particle[i].omega + particle[j].omega);

        double ftx = -(this->parameters[4]) * 1 * (omega_sum / 2) * r2 +
                     (this->parameters[5]) * (vrel_x - vreldotr * r1);

        double fty = -(this->parameters[4]) * 1 * (omega_sum / 2) * (-r1) +
                     (this->parameters[5]) * (vrel_y - vreldotr * r2);

        particle[i].force_tangential[0] += ftx;

        particle[i].force_tangential[1] += fty;

        // torque on particle

        if (omega_sum != 0) {
          particle[i].torque +=
              -(this->parameters[4]) * 1 *
              (omega_sum / (abs(omega_sum) + (this->parameters[7]))) * d;
        }

        if (omega_sum == 0) {
          particle[i].torque +=
              -(this->parameters[4]) * 1 *
              (omega_sum / (abs(omega_sum) + (this->parameters[7]))) * d;
        }

        // pressure term calculation

        particle[i].sigma[0] += (frx + ftx) * (particle[i].x - particle[j].x);
        particle[i].sigma[1] += (frx + ftx) * (particle[i].y - particle[j].y);
        particle[i].sigma[2] += (fry + fty) * (particle[i].x - particle[j].x);
        particle[i].sigma[3] += (fry + fty) * (particle[i].y - particle[j].y);
      }
    }
  }
}

void ParSim::Physics::Force_PP_CRB_WCA2(ParSim::ParticleSystem &ps) {

  Particle *const particle = ps.get_particles();

  // Loop1: through all particles----------------
  for (int i = 0; i < ps.no_of_particles; ++i) {

    // Store previous step forces
    particle[i].force_radial_prev[0] = particle[i].force_radial[0];
    particle[i].force_radial_prev[1] = particle[i].force_radial[1];
    particle[i].force_tangential_prev[0] = particle[i].force_tangential[0];
    particle[i].force_tangential_prev[1] = particle[i].force_tangential[1];
    particle[i].torque_prev = particle[i].torque;

    // Reset the current forces and pressures
    particle[i].force_radial[0] = 0;
    particle[i].force_radial[1] = 0;
    particle[i].force_tangential[0] = 0;
    particle[i].force_tangential[1] = 0;
    particle[i].torque = 0;

    particle[i].sigma[0] = 0.0;
    particle[i].sigma[1] = 0.0;
    particle[i].sigma[2] = 0.0;
    particle[i].sigma[3] = 0.0;

    // wall forces -----------------
    double force_wall_x = 0.0;
    double force_wall_y = 0.0;

    double force_wall_fric_x = 0.0;
    double force_wall_fric_y = 0.0;

    double torque_wall = 0.0;

    double rp = 0.0;

    double magnitude = 0.0;

    rp = std::sqrt(particle[i].x * particle[i].x +
                   particle[i].y * particle[i].y) +
         (this->parameters[1] /
          2); // calculates position vector magnitude + (sigma/2)

    double radius = ps.L;

    if (rp > radius) {

      double flag = ((this->parameters[1]) -
                     abs(rp - radius)) / // proxy for (r/sigma) in WCA
                    this->circular_wall_parameters[1];
      double flag5 = flag * flag * flag * flag * flag;

      double norm =
          4 *
          (this->circular_wall_parameters[0] /
           this->circular_wall_parameters[1]) *
          (1.0 / (flag5 * flag * flag + (this->parameters[6]))) *
          ((12.0 / (flag5 + (this->parameters[6]))) - (6.0)); // WCA force

      force_wall_x = -norm * particle[i].x / (rp - this->parameters[1]);

      force_wall_y = -norm * particle[i].y / (rp - this->parameters[1]);

      force_wall_fric_x = -(this->parameters[4] / 2.0) * abs(norm) *
                          (particle[i].omega /
                           (abs(particle[i].omega + (this->parameters[7])))) *
                          (-(particle[i].y) / (rp - this->parameters[1]));

      force_wall_fric_y = -(this->parameters[4] / 2.0) * abs(norm) *
                          (particle[i].omega /
                           (abs(particle[i].omega + (this->parameters[7])))) *
                          ((particle[i].x) / (rp - this->parameters[1]));

    } else {
      force_wall_x = 0.0;
      force_wall_fric_x = 0.0;

      force_wall_y = 0.0;
      force_wall_fric_y = 0.0;
    }

    // wall torque

    // if (rp > radius) {
    //   torque_wall =
    //       -(this->parameters[4]/2.0) * abs(norm) *
    //       (particle[i].omega /
    //        (abs(particle[i].omega + (this->parameters[7])))) *
    //       ((ps.L / 2) - (rp - this->parameters[1]));

    // } else {
    //   torque_wall = 0.0;
    // }

    // wall friction reduced to null

    torque_wall = 0.0;
    force_wall_fric_x = 0.0;
    force_wall_fric_y = 0.0;

    // Unary force of damping. Always there. Translational and Rotational
    // activities added. (note these are not purely radial and tangential)

    particle[i].force_radial[0] += -1 * (this->parameters[5]) * particle[i].vx +
                                   force_wall_x + force_wall_fric_x;

    particle[i].force_radial[1] += -1 * (this->parameters[5]) * particle[i].vy +
                                   force_wall_y + force_wall_fric_y;

    particle[i].torque +=
        -1 * (this->parameters[10]) * particle[i].omega +
        (particle[i].omega_activity * (this->parameters[10])) + torque_wall;
  }

  // Loop 2: N(N-1)/2 Loop for pair force calculation ------
  for (int i = 0; i < (ps.no_of_particles - 1); ++i) {

    for (int j = i + 1; j < ps.no_of_particles; ++j) {

      // double d = ps.distance(par[i], par[j]);

      double d = sqrt(
          (particle[i].x - particle[j].x) * (particle[i].x - particle[j].x) +
          (particle[i].y - particle[j].y) * (particle[i].y - particle[j].y));

      // U

      if (d <= 1.12 * (this->parameters[1])) {

        // relative velocities calculation (ith particle at rest)
        double vrel_x = 0; // particle[j].vx - particle[i].vx;
        double vrel_y = 0; // particle[j].vy - particle[i].vy;

        // mag. of radial force  (positive)
        double r_by_sigma = d / this->WCA_parameters[1]; // r/sigma in WCA
        double r_by_sigma5 =
            r_by_sigma * r_by_sigma * r_by_sigma * r_by_sigma * r_by_sigma;

        double N =
            4 * (this->WCA_parameters[0] / this->WCA_parameters[1]) *
            (1.0 /
             (r_by_sigma5 * r_by_sigma5 * r_by_sigma * r_by_sigma * r_by_sigma +
              (this->parameters[6]))) *
            ((24.0 / (r_by_sigma5 * r_by_sigma5 * r_by_sigma * r_by_sigma +
                      (this->parameters[6]))) -
             (12.0)); // WCA force magnitude (positive)

        // components of unit separation vector r
        double r1 =
            (particle[i].x - particle[j].x) / (d + (this->parameters[6])); // rx
        double r2 =
            (particle[i].y - particle[j].y) / (d + (this->parameters[6])); // ry

        double vreldotr = (r1 * vrel_x + r2 * vrel_y);
        // Radial interaction force with jth particle
        double frx = N * r1 + (this->parameters[5]) * vreldotr * r1;
        double fry = N * r2 + (this->parameters[5]) * vreldotr * r2;

        particle[i].force_radial[0] += frx;
        particle[i].force_radial[1] += fry;

        particle[j].force_radial[0] += -frx;
        particle[j].force_radial[1] += -fry;

        // Tangential friction force with particle j
        double omega_sum = (particle[i].omega + particle[j].omega);

        double ftx = -(this->parameters[4]) * 1 * (omega_sum / 2) * r2 +
                     (this->parameters[5]) * (vrel_x - vreldotr * r1);

        double fty = -(this->parameters[4]) * 1 * (omega_sum / 2) * (-r1) +
                     (this->parameters[5]) * (vrel_y - vreldotr * r2);

        particle[i].force_tangential[0] += ftx;
        particle[i].force_tangential[1] += fty;

        particle[j].force_tangential[0] += -ftx;
        particle[j].force_tangential[1] += -fty;

        // torque on particle

        if (omega_sum != 0) {
          particle[i].torque +=
              -(this->parameters[4]) * 1 *
              (omega_sum / (abs(omega_sum) + (this->parameters[7]))) * d;
          particle[j].torque +=
              -(this->parameters[4]) *
              ((this->parameters[0]) * ((this->parameters[1]) - d)) *
              (omega_sum / 2.0) * d;
        }

        if (omega_sum == 0) {
          particle[i].torque +=
              -(this->parameters[4]) * 1 *
              (omega_sum / (abs(omega_sum) + (this->parameters[7]))) * d;
          particle[j].torque +=
              -(this->parameters[4]) *
              ((this->parameters[0]) * ((this->parameters[1]) - d)) *
              (omega_sum / 2.0) * d;
        }

        // pressure term calculation

        particle[i].sigma[0] += -(frx + ftx) * (particle[i].x - particle[j].x);
        particle[i].sigma[1] += -(frx + ftx) * (particle[i].y - particle[j].y);
        particle[i].sigma[2] += -(fry + fty) * (particle[i].x - particle[j].x);
        particle[i].sigma[3] += -(fry + fty) * (particle[i].y - particle[j].y);

        particle[j].sigma[0] += -(frx + ftx) * (particle[i].x - particle[j].x);
        particle[j].sigma[1] += -(frx + ftx) * (particle[i].y - particle[j].y);
        particle[j].sigma[2] += -(fry + fty) * (particle[i].x - particle[j].x);
        particle[j].sigma[3] += -(fry + fty) * (particle[i].y - particle[j].y);
      }
    }
  }
}
// 2) Neighbour list force linkers ----------------

void ParSim::Physics::Neighbours_search_PBC(ParticleSystem &ps) {

  Particle *const particle = ps.get_particles();

  // Loop1: through all particles
  for (int i = 0; i < ps.no_of_particles; ++i) {

    // Clear the present list
    particle[i].Neighbours.clear(); // this function takes time

    for (int j = 0; j < ps.no_of_particles; ++j) {
      // Distance from the nearest image of jth particle

      // i)calculate min separation dx and dy
      double dx = particle[i].x - particle[j].x;

      if (dx < (-ps.L / 2)) {
        dx += ps.L;
      } else if (dx > (ps.L / 2)) {
        dx -= ps.L;
      }

      double dy = particle[i].y - particle[j].y;

      if (dy < (-ps.L / 2)) {
        dy += ps.L;
      } else if (dy > (ps.L / 2)) {
        dy -= ps.L;
      }
      // ii) Nearest image distance
      double d = std::sqrt(dx * dx + dy * dy);

      if (d <= (this->parameters[1]) + this->NL_skin_depth) { // rs = d + 2.0
        particle[i].Neighbours.push_back(
            j); // j is a neighbour of i and also vice versa
      }
    }
  }
}

void ParSim::Physics::Neighbours_search(ParticleSystem &ps) {

  Particle *const particle = ps.get_particles();

  // Loop1: through all particles
  for (int i = 0; i < ps.no_of_particles; ++i) {

    // Clear the present list
    particle[i].Neighbours.clear(); // this function takes time

    for (int j = 0; j < ps.no_of_particles; ++j) {

      // i) Calculate separation dx and dy
      double dx = particle[i].x - particle[j].x;
      double dy = particle[i].y - particle[j].y;

      // ii) Particle distance
      double d = std::sqrt(dx * dx + dy * dy);

      if (d <= (this->parameters[1]) + (this->NL_skin_depth)) { // rs = d + 2.0
        particle[i].Neighbours.push_back(
            j); // j is a neighbour of i and also vice versa
      }
    }
  }
}

void ParSim::Physics::Force_NL_PBC(ParSim::ParticleSystem &ps) {

  Particle *const particle = ps.get_particles();

  // Loop1: through all particles
  for (int i = 0; i < ps.no_of_particles; ++i) {

    // Store previous step forces
    particle[i].force_radial_prev[0] = particle[i].force_radial[0];
    particle[i].force_radial_prev[1] = particle[i].force_radial[1];
    particle[i].force_tangential_prev[0] = particle[i].force_tangential[0];
    particle[i].force_tangential_prev[1] = particle[i].force_tangential[1];
    particle[i].torque_prev = particle[i].torque;

    // Reset the current forces and pressures
    particle[i].force_radial[0] = 0;
    particle[i].force_radial[1] = 0;
    particle[i].force_tangential[0] = 0;
    particle[i].force_tangential[1] = 0;
    particle[i].torque = 0;

    particle[i].sigma[0] = 0.0;
    particle[i].sigma[1] = 0.0;
    particle[i].sigma[2] = 0.0;
    particle[i].sigma[3] = 0.0;

    // Unary force of damping. Always there. Translational and Rotational
    // activities added (note this is not radial !!!).

    particle[i].force_radial[0] +=
        0; //-1 * (this->parameters[5]) * particle[i].vx;

    particle[i].force_radial[1] +=
        0; //-1 * (this->parameters[5]) * particle[i].vy;

    particle[i].torque += -1 * (this->parameters[10]) * particle[i].omega +
                          (particle[i].omega_activity * (this->parameters[10]));

    // Force calculation for Neighbour list particles
    std::list<int>::iterator j;

    for (j = particle[i].Neighbours.begin(); j != particle[i].Neighbours.end();
         j++) {

      // Distance from the nearest image of jth particle

      // i)calculate min separation dx and dy
      double dx = particle[i].x - particle[*j].x;

      if (dx < (-ps.L / 2)) {
        dx += ps.L;
      } else if (dx > (ps.L / 2)) {
        dx -= ps.L;
      }

      double dy = particle[i].y - particle[*j].y;

      if (dy < (-ps.L / 2)) {
        dy += ps.L;
      } else if (dy > (ps.L / 2)) {
        dy -= ps.L;
      }
      // ii) Nearest image distance
      double d = std::sqrt(dx * dx + dy * dy);

      // U

      if (d <= (this->parameters[1])) {

        // relative velocities calculation (ith particle at rest)
        double vrel_x = particle[*j].vx - particle[i].vx;
        double vrel_y = particle[*j].vy - particle[i].vy;

        // mag. of radial force
        double N = (this->parameters[0]) *
                   ((this->parameters[1]) -
                    d); // magnitude of radial force used as normal reaction

        // components of unit separation vector r

        double r1 = dx / (d + (this->parameters[6]));
        double r2 = dy / (d + (this->parameters[6]));

        double vreldotr = (r1 * vrel_x + r2 * vrel_y);
        // Radial interaction force with jth particle
        double frx = N * r1 + (this->parameters[5]) * vreldotr * r1;
        double fry = N * r2 + (this->parameters[5]) * vreldotr * r2;

        particle[i].force_radial[0] += frx;
        particle[i].force_radial[1] += fry;
        // Tangential friction force with particle j
        double omega_sum = (particle[i].omega + particle[*j].omega);

        double ftx = -(this->parameters[4]) *
                         ((this->parameters[0]) * ((this->parameters[1]) - d)) *
                         (omega_sum / 2.0) * r2 +
                     (this->parameters[5]) * (vrel_x - vreldotr * r1);

        double fty = -(this->parameters[4]) *
                         ((this->parameters[0]) * ((this->parameters[1]) - d)) *
                         (omega_sum / 2.0) * (-r1) +
                     (this->parameters[5]) * (vrel_y - vreldotr * r2);

        particle[i].force_tangential[0] += ftx;
        particle[i].force_tangential[1] += fty;

        // torque on particle (generic torque, should be written as cross
        // product instead )

        if (omega_sum != 0) {
          particle[i].torque +=
              -(this->parameters[4]) *
              ((this->parameters[0]) * ((this->parameters[1]) - d)) *
              (omega_sum / 2.0) * d;
        }

        if (omega_sum == 0) {
          particle[i].torque +=
              -(this->parameters[4]) *
              ((this->parameters[0]) * ((this->parameters[1]) - d)) *
              (omega_sum / 2.0) * d;
        }

        // pressure term calculation

        particle[i].sigma[0] += (frx + ftx) * (particle[i].x - particle[*j].x);
        particle[i].sigma[1] += (frx + ftx) * (particle[i].y - particle[*j].y);
        particle[i].sigma[2] += (fry + fty) * (particle[i].x - particle[*j].x);
        particle[i].sigma[3] += (fry + fty) * (particle[i].y - particle[*j].y);
      }
    }
  }
}

void ParSim::Physics::Force_NL_CRB_WCA(ParSim::ParticleSystem &ps) {

  Particle *const particle = ps.get_particles();

  // Loop1: through all particles----------------
  for (int i = 0; i < ps.no_of_particles; ++i) {

    // Store previous step forces
    particle[i].force_radial_prev[0] = particle[i].force_radial[0];
    particle[i].force_radial_prev[1] = particle[i].force_radial[1];
    particle[i].force_tangential_prev[0] = particle[i].force_tangential[0];
    particle[i].force_tangential_prev[1] = particle[i].force_tangential[1];
    particle[i].torque_prev = particle[i].torque;

    // Reset the current forces and pressures
    particle[i].force_radial[0] = 0;
    particle[i].force_radial[1] = 0;
    particle[i].force_tangential[0] = 0;
    particle[i].force_tangential[1] = 0;
    particle[i].torque = 0;

    particle[i].sigma[0] = 0.0;
    particle[i].sigma[1] = 0.0;
    particle[i].sigma[2] = 0.0;
    particle[i].sigma[3] = 0.0;

    // wall forces -----------------
    double force_wall_x = 0.0;
    double force_wall_y = 0.0;

    double force_wall_fric_x = 0.0;
    double force_wall_fric_y = 0.0;

    double torque_wall = 0.0;

    double rp = 0.0;

    double magnitude = 0.0;

    rp = std::sqrt(particle[i].x * particle[i].x +
                   particle[i].y * particle[i].y) +
         (this->parameters[1] /
          2); // calculates position vector magnitude + (sigma/2)

    double radius = ps.L;

    if (rp > radius) {

      double flag = ((this->parameters[1]) -
                     abs(rp - radius)) / // proxy for (r/sigma) in WCA
                    this->circular_wall_parameters[1];
      double flag5 = flag * flag * flag * flag * flag;

      double norm =
          4 *
          (this->circular_wall_parameters[0] /
           this->circular_wall_parameters[1]) *
          (1.0 / (flag5 * flag * flag + (this->parameters[6]))) *
          ((12.0 / (flag5 + (this->parameters[6]))) - (6.0)); // WCA force

      force_wall_x = -norm * particle[i].x / (rp - this->parameters[1]);

      force_wall_y = -norm * particle[i].y / (rp - this->parameters[1]);

      force_wall_fric_x = -(this->parameters[4] / 2.0) * abs(norm) *
                          (particle[i].omega /
                           (abs(particle[i].omega + (this->parameters[7])))) *
                          (-(particle[i].y) / (rp - this->parameters[1]));

      force_wall_fric_y = -(this->parameters[4] / 2.0) * abs(norm) *
                          (particle[i].omega /
                           (abs(particle[i].omega + (this->parameters[7])))) *
                          ((particle[i].x) / (rp - this->parameters[1]));

    } else {
      force_wall_x = 0.0;
      force_wall_fric_x = 0.0;

      force_wall_y = 0.0;
      force_wall_fric_y = 0.0;
    }

    // wall torque

    // if (rp > radius) {
    //   torque_wall =
    //       -(this->parameters[4]/2.0) * abs(norm) *
    //       (particle[i].omega /
    //        (abs(particle[i].omega + (this->parameters[7])))) *
    //       ((ps.L / 2) - (rp - this->parameters[1]));

    // } else {
    //   torque_wall = 0.0;
    // }

    // wall friction reduced to null

    torque_wall = 0.0;
    force_wall_fric_x = 0.0;
    force_wall_fric_y = 0.0;

    // Unary force of damping. Always there. Translational and Rotational
    // activities added. (note these are not purely radial and tangential)

    particle[i].force_radial[0] += -1 * (this->parameters[5]) * particle[i].vx +
                                   force_wall_x + force_wall_fric_x;

    particle[i].force_radial[1] += -1 * (this->parameters[5]) * particle[i].vy +
                                   force_wall_y + force_wall_fric_y;

    particle[i].torque +=
        -1 * (this->parameters[10]) * particle[i].omega +
        (particle[i].omega_activity * (this->parameters[10])) + torque_wall;

    // Force calculation for Neighbour list particles
    std::list<int>::iterator j;

    for (j = particle[i].Neighbours.begin(); j != particle[i].Neighbours.end();
         j++) {

      // double d = ps.distance(par[i], par[j]);

      double d = sqrt(
          (particle[i].x - particle[*j].x) * (particle[i].x - particle[*j].x) +
          (particle[i].y - particle[*j].y) * (particle[i].y - particle[*j].y));

      // U

      if (d <= 1.12 * (this->parameters[1])) {

        // relative velocities calculation (ith particle at rest)
        double vrel_x = 0; // particle[j].vx - particle[i].vx;
        double vrel_y = 0; // particle[j].vy - particle[i].vy;

        // mag. of radial force  (positive)
        double r_by_sigma = d / this->WCA_parameters[1]; // r/sigma in WCA
        double r_by_sigma5 =
            r_by_sigma * r_by_sigma * r_by_sigma * r_by_sigma * r_by_sigma;

        double N =
            4 * (this->WCA_parameters[0] / this->WCA_parameters[1]) *
            (1.0 /
             (r_by_sigma5 * r_by_sigma5 * r_by_sigma * r_by_sigma * r_by_sigma +
              (this->parameters[6]))) *
            ((24.0 / (r_by_sigma5 * r_by_sigma5 * r_by_sigma * r_by_sigma +
                      (this->parameters[6]))) -
             (12.0)); // WCA force magnitude (positive)

        // components of unit separation vector r
        double r1 = (particle[i].x - particle[*j].x) /
                    (d + (this->parameters[6])); // rx
        double r2 = (particle[i].y - particle[*j].y) /
                    (d + (this->parameters[6])); // ry

        double vreldotr = (r1 * vrel_x + r2 * vrel_y);
        // Radial interaction force with jth particle
        double frx = N * r1 + (this->parameters[5]) * vreldotr * r1;
        double fry = N * r2 + (this->parameters[5]) * vreldotr * r2;

        particle[i].force_radial[0] += frx;
        particle[i].force_radial[1] += fry;

        // Tangential friction force with particle j
        double omega_sum = (particle[i].omega + particle[*j].omega);

        double ftx = -(this->parameters[4]) * 1 * (omega_sum / 2.0) * r2 +
                     (this->parameters[5]) * (vrel_x - vreldotr * r1);

        double fty = -(this->parameters[4]) * 1 * (omega_sum / 2.0) * (-r1) +
                     (this->parameters[5]) * (vrel_y - vreldotr * r2);

        particle[i].force_tangential[0] += ftx;
        particle[i].force_tangential[1] += fty;

        // torque on particle

        if (omega_sum != 0) {
          particle[i].torque +=
              -(this->parameters[4]) * 1 * (omega_sum / 2.0) * d;
        }

        if (omega_sum == 0) {
          particle[i].torque +=
              -(this->parameters[4]) * 1 * (omega_sum / 2.0) * d;
        }

        // pressure term calculation

        particle[i].sigma[0] += -(frx + ftx) * (particle[i].x - particle[*j].x);
        particle[i].sigma[1] += -(frx + ftx) * (particle[i].y - particle[*j].y);
        particle[i].sigma[2] += -(fry + fty) * (particle[i].x - particle[*j].x);
        particle[i].sigma[3] += -(fry + fty) * (particle[i].y - particle[*j].y);
      }
    }
  }
}

// 3) Integrators ----------------------------
void ParSim::Physics::Euler_Integrator(ParSim::Particle &par, int step) {

  double m = this->parameters[2];
  double time_step = this->parameters[8];
  double Fx;
  double Fy;
  double Tau;
  double dvx;
  double dvy;
  double dw;
  // Total force and torque on this particle
  Fx = par.force_radial[0] + par.force_tangential[0];
  Fy = par.force_radial[1] + par.force_tangential[1]; // -g newton, m = 1
  // log << "Ftangential: " << par.force_tangential[1] << std::endl;
  Tau = par.torque;

  dvx = (Fx / m) * time_step;
  dvy = (Fy / m) * time_step;
  // log << "dvy: " << dvy << std::endl;
  dw = (Tau)*time_step;

  // update the attributes

  par.x += par.vx * time_step;
  par.y += par.vy * time_step;
  par.vx += dvx;
  par.vy += dvy;

  // log << "vy: " << par.vy << std::endl;

  par.alpha += par.omega * time_step;
  par.omega += dw;
}

void ParSim::Physics::Vel_Verlet_Integrator(ParSim::Particle &par, int step) {

  double m = this->parameters[2];
  double time_step = this->parameters[8];

  // update the attributes

  par.x += (time_step * par.vx) +
           (pow(time_step, 2)) *
               (par.force_radial[0] + par.force_tangential[0]) / (2 * m);
  par.y += (time_step * par.vy) +
           (pow(time_step, 2)) *
               (par.force_radial[1] + par.force_tangential[1]) / (2 * m);

  par.alpha +=
      (time_step * par.omega) + (pow(time_step, 2)) * (par.torque) / (2 * m);

  if (step != 0) {
    par.vx += time_step *
              ((par.force_radial[0] + par.force_tangential[0]) +
               (par.force_radial_prev[0] + par.force_tangential_prev[0])) /
              (2 * m);
    par.vy += time_step *
              ((par.force_radial[1] + par.force_tangential[1]) +
               (par.force_radial_prev[1] + par.force_tangential_prev[1])) /
              (2 * m);
    par.omega += time_step * (par.torque + par.torque_prev) / (2 * m);
  }

  // log << "vy: " << par.vy << std::endl;
}

void ParSim::Physics::ERM_Integrator1(ParSim::Particle &par, int step) {

  double Fx;
  double Fy;
  double Tau;
  // Total force and torque on this particle
  Fx = par.force_radial[0] + par.force_tangential[0];
  Fy = par.force_radial[1] + par.force_tangential[1]; //

  Tau = par.torque;

  // Save present attributes
  par.position_prev[0] = par.x;
  par.position_prev[1] = par.y;

  par.velocity_prev[0] = par.vx;
  par.velocity_prev[1] = par.vy;

  par.alpha_prev = par.alpha;
  par.omega_prev = par.omega;

  // update the attributes upto midpoint (x', v')

  par.x += par.vx * (this->parameters[8]) / 2; // x'
  par.y += par.vy * (this->parameters[8]) / 2;
  par.vx += (Fx / (this->parameters[2])) * (this->parameters[8]) / 2; // v'
  par.vy += (Fy / (this->parameters[2])) * (this->parameters[8]) / 2;

  par.alpha += par.omega * (this->parameters[8]) / 2;
  par.omega += (Tau / (this->parameters[2])) * (this->parameters[8]) / 2;

  // For random active force -- use Fundamental Euler integration to update by
  // half step
  std::normal_distribution<double> distribution(0.0, 1.0);

  par.theta += pow(this->parameters[11], 0.5) * distribution(this->mt) *
               (this->parameters[8] / 2); // rotational diffusion
}

void ParSim::Physics::ERM_Integrator2(ParSim::Particle &par, double L,
                                      int step) {

  double Fx;
  double Fy;
  double Tau;
  // Total force (F') and torque (Tau') on this particle
  Fx = par.force_radial[0] + par.force_tangential[0];
  Fy = par.force_radial[1] + par.force_tangential[1];

  Tau = par.torque;

  // update the attributes to next time step (x1,v1)

  par.x = par.position_prev[0] + par.vx * (this->parameters[8]); // x1
  par.y = par.position_prev[1] + par.vy * (this->parameters[8]);
  par.vx = par.velocity_prev[0] +
           (Fx / (this->parameters[2])) * (this->parameters[8]); // v1
  par.vy = par.velocity_prev[1] +
           (Fy / (this->parameters[2])) * (this->parameters[8]);

  par.alpha = par.alpha_prev + par.omega * (this->parameters[8]);
  par.omega =
      par.omega_prev + (Tau / (this->parameters[2])) * (this->parameters[8]);

  // For random active force -- use Fundamental Euler integration to update by
  // half step
  std::normal_distribution<double> distribution(0.0, 1.0);

  par.theta += pow(this->parameters[11], 0.5) * distribution(this->mt) *
               (this->parameters[8] / 2); // rotational diffusion

  // Periodic boudary condition

  // if (par.x > L / 2) {
  //   par.x -= L;
  // } else if (par.x < -L / 2) {
  //   par.x += L;
  // }

  // if (par.y > L / 2) {
  //   par.y -= L;
  // } else if (par.y < -L / 2) {
  //   par.y += L;
  // }
}

void ParSim::Physics::ERM_Integrator1_sys(ParSim::ParticleSystem &parsym,
                                          int step) {

  double Fx = 0.0;
  double Fy = 0.0;
  double Tau = 0.0;
  for (int i = 0; i < parsym.no_of_particles; ++i) {

    // reset forces for this particle
    Fx = 0.0;
    Fy = 0.0;
    Tau = 0.0;

    // Total force and torque on this ith particle
    Fx = parsym.particle_array[i].force_radial[0] +
         parsym.particle_array[i].force_tangential[0];
    Fy = parsym.particle_array[i].force_radial[1] +
         parsym.particle_array[i].force_tangential[1]; //

    Tau = parsym.particle_array[i].torque;

    // Save present attributes
    parsym.particle_array[i].position_prev[0] = parsym.particle_array[i].x;
    parsym.particle_array[i].position_prev[1] = parsym.particle_array[i].y;

    parsym.particle_array[i].velocity_prev[0] = parsym.particle_array[i].vx;
    parsym.particle_array[i].velocity_prev[1] = parsym.particle_array[i].vy;

    parsym.particle_array[i].alpha_prev = parsym.particle_array[i].alpha;
    parsym.particle_array[i].omega_prev = parsym.particle_array[i].omega;

    // update the attributes upto midpoint (x', v')

    parsym.particle_array[i].x +=
        parsym.particle_array[i].vx * (this->parameters[8]) / 2; // x'
    parsym.particle_array[i].y +=
        parsym.particle_array[i].vy * (this->parameters[8]) / 2;
    parsym.particle_array[i].vx +=
        (Fx / (this->parameters[2])) * (this->parameters[8]) / 2; // v'
    parsym.particle_array[i].vy +=
        (Fy / (this->parameters[2])) * (this->parameters[8]) / 2;

    parsym.particle_array[i].alpha +=
        parsym.particle_array[i].omega * (this->parameters[8]) / 2;
    parsym.particle_array[i].omega +=
        (2 * Tau / (this->parameters[2])) * (this->parameters[8]) / 2;

    // For random active force -- use Fundamental Euler integration to update by
    // half step
    // std::normal_distribution<double> distribution(0.0, 1.0);

    // parsym.particle_array[i].theta +=
    //     pow(this->parameters[11], 0.5) * distribution(this->mt) *
    //     (this->parameters[8] / 2); // rotational diffusion
  }
}

void ParSim::Physics::ERM_Integrator2_sys(ParSim::ParticleSystem &parsym,
                                          int step) {

  double Fx = 0.0;
  double Fy = 0.0;
  double Tau = 0.0;

  double noise_Fx = 0.0;
  double noise_Fy = 0.0;
  double noise_Tau = 0.0;

  for (int i = 0; i < parsym.no_of_particles; ++i) {
    // reset forces for this particle
    Fx = 0.0;
    Fy = 0.0;
    Tau = 0.0;

    // Total force and torque on this ith particle
    Fx = parsym.particle_array[i].force_radial[0] +
         parsym.particle_array[i].force_tangential[0];
    Fy = parsym.particle_array[i].force_radial[1] +
         parsym.particle_array[i].force_tangential[1]; //

    Tau = parsym.particle_array[i].torque;

    // update the attributes to next time step (x1,v1)

    parsym.particle_array[i].x =
        parsym.particle_array[i].position_prev[0] +
        parsym.particle_array[i].vx * (this->parameters[8]); // x1
    parsym.particle_array[i].y =
        parsym.particle_array[i].position_prev[1] +
        parsym.particle_array[i].vy * (this->parameters[8]);
    parsym.particle_array[i].vx =
        parsym.particle_array[i].velocity_prev[0] +
        (Fx / (this->parameters[2])) * (this->parameters[8]); // v1
    parsym.particle_array[i].vy =
        parsym.particle_array[i].velocity_prev[1] +
        (Fy / (this->parameters[2])) * (this->parameters[8]);

    parsym.particle_array[i].alpha =
        parsym.particle_array[i].alpha_prev +
        parsym.particle_array[i].omega * (this->parameters[8]);
    parsym.particle_array[i].omega =
        parsym.particle_array[i].omega_prev +
        (2 * Tau / (this->parameters[2])) * (this->parameters[8]);

    // For random active force -- use Fundamental Euler integration to update by
    // half step

    // Noise in the system

    std::normal_distribution<double> distribution(0.0, 1.0); // for noise

    noise_Fx = 1 *
               sqrt(2 * this->parameters[5] *
                    this->parameters[12]) * // variance = 2(gamma_t)kbT
               distribution(this->mt);

    noise_Fy = 1 *
               sqrt(2 * this->parameters[5] *
                    this->parameters[12]) * // variance = 2(gamma_t)kbT
               distribution(this->mt);
    noise_Tau = 1 *
                sqrt(2 * this->parameters[10] *
                     this->parameters[12]) * // variance = 2(gamma_r)kbT
                distribution(this->mt);

    parsym.particle_array[i].vx +=
        (noise_Fx / this->parameters[2]) * (this->parameters[8]);

    parsym.particle_array[i].vy +=
        (noise_Fy / this->parameters[2]) * (this->parameters[8]);

    parsym.particle_array->omega +=
        (2 * noise_Tau / this->parameters[2]) * (this->parameters[8]);

    // Periodic boudary condition

    // if (parsym.particle_array[i].x > parsym.L / 2) {
    //   parsym.particle_array[i].x -= parsym.L;
    // } else if (parsym.particle_array[i].x < -parsym.L / 2) {
    //   parsym.particle_array[i].x += parsym.L;
    // }
    // if (parsym.particle_array[i].y > parsym.L / 2) {
    //   parsym.particle_array[i].y -= parsym.L;
    // } else if (parsym.particle_array[i].y < -parsym.L / 2) {
    //   parsym.particle_array[i].y += parsym.L;
    // }
  }
}

// 4) Evolvers --------------------------------

void ParSim::Physics::evolve_system(ParticleSystem &parsym, int step) {

  // 1)Force-linking--------
  ParSim::Physics::Force_PP(parsym); // links forces on each object
                                     // at this stage
  // 2)Integrating----------
  for (int i = 0; i < parsym.no_of_particles; ++i) {
    // Vel_Verlet_Integrator(parsym.particle_array[i], step, log);
    ERM_Integrator1(parsym.particle_array[i], step);
    //  boundary conditions

    // if (parsym.particle_array[i].x < -1000 ||
    //     parsym.particle_array[i].x > 1000 ||
    //     parsym.particle_array[i].y < -800 || parsym.particle_array[i].y >
    //     800)
    //   parsym.particle_array[i].random_initialize();
  }
}

void ParSim::Physics::evolve_system_ERM(ParticleSystem &parsym, int step) {

  // i) Calculate force (F) from positions and velocities (x,v) --
  Force_NL_PBC(parsym); // links forces on each object

  // ii) Update x,v to x', v'  ----
  ERM_Integrator1_sys(parsym, step);

  // iii) Again calculate force (F') from x',v' ------
  Force_NL_PBC(parsym);

  //  iv) Update x',v' to xnew, vnew --------
  ERM_Integrator2_sys(parsym, step);
}

void ParSim::Physics::evolve_system_ERM_NL(ParticleSystem &parsym, int step) {

  // Update Neighbour list after every 200 steps
  if (step == 0 || step % this->NL_tau == 0) {
    Neighbours_search(parsym);
  }
  // i) Calculate force (F) from positions and velocities (x,v) --
  Force_NL_CRB_WCA(parsym); // links forces on each object

  // ii) Update x,v to x', v'  ----
  ERM_Integrator1_sys(parsym, step);

  // iii) Again calculate force (F') from x',v' ------
  Force_NL_CRB_WCA(parsym);

  //  iv) Update x',v' to xnew, vnew --------
  ERM_Integrator2_sys(parsym, step);
}

/*Conserved Quantities*/
std::vector<double> ParSim::Physics::EnergyMomentum(ParticleSystem &parsym) {
  std::vector<double> p{0, 0, 0, 0}; // kinetic, rotational, momenta
  double m = this->parameters[2];
  double I = m;
  for (int i = 0; i < parsym.no_of_particles; ++i) {
    p[0] += 0.5 * m *
            (pow(parsym.particle_array[i].vx, 2) +
             pow(parsym.particle_array[i].vy, 2)); // KE
    p[1] +=
        0.5 * I * (pow(parsym.particle_array[i].omega, 2)); // Rotational K.E
    p[2] += m * parsym.particle_array[i].vx;                // px
    p[3] += m * parsym.particle_array[i].vy;                // py
  }
  return p;
}

// #######################################################//
//----------------------------------------
//
//  Class: Langevin_Dyanmics definitions
//
//--------------------------------------
// #######################################################//

ParSim::Langevin_Dynamics::Langevin_Dynamics()
    : mt((std::random_device())()) { // initialize seed of mt using random
                                     // device

  parameters[0] = 0.0;
  parameters[1] = 0.0;
  parameters[2] = 0.0;
  parameters[3] = 0.0;
  parameters[4] = 0.0;
  parameters[5] = 0.0;
  parameters[6] = 0.0;
  parameters[7] = 0.0;
  parameters[8] = 0.0;
  parameters[9] = 0.0;
  parameters[10] = 0.0;
  parameters[11] = 0.0;
  parameters[12] = 0.0;
  parameters[13] = 0.0; // I

  // Fluctuation-Dissipation Parameters
  FD_parameters[0] = 0.0; // jeta_t
  FD_parameters[1] = 0.0; // jeta_r
  FD_parameters[2] = 0.0; // kbt
  FD_parameters[3] = 0.0; // D

  WCA_parameters[0] = 0.0; // WCA epsilon
  WCA_parameters[1] = 0.0; // WCA sigma

  circular_wall_parameters[0] = 0.0; // WCA epsilon
  circular_wall_parameters[1] = 0.0; // WCA sigma

  NL_skin_depth = 0.0; // DeltaR
  NL_tau = 0;          // tau
}

// 2) Neighbour list force linkers ----------------

void ParSim::Langevin_Dynamics::Neighbours_search_PBC(ParticleSystem &ps) {

  Particle *const particle = ps.get_particles();

  // Loop1: through all particles
  for (int i = 0; i < ps.no_of_particles; ++i) {

    // Clear the present list
    particle[i].Neighbours.clear(); // this function takes time

    for (int j = 0; j < ps.no_of_particles; ++j) {
      // Distance from the nearest image of jth particle

      // i)calculate min separation dx and dy
      double dx = particle[i].x - particle[j].x;

      if (dx < (-ps.L / 2)) {
        dx += ps.L;
      } else if (dx > (ps.L / 2)) {
        dx -= ps.L;
      }

      double dy = particle[i].y - particle[j].y;

      if (dy < (-ps.L / 2)) {
        dy += ps.L;
      } else if (dy > (ps.L / 2)) {
        dy -= ps.L;
      }
      // ii) Nearest image distance
      double d = std::sqrt(dx * dx + dy * dy);

      if (d <= (this->parameters[1]) + this->NL_skin_depth) { // rs = d + 2.0
        particle[i].Neighbours.push_back(
            j); // j is a neighbour of i and also vice versa
      }
    }
  }
}

void ParSim::Langevin_Dynamics::Neighbours_search(ParticleSystem &ps) {

  Particle *const particle = ps.get_particles();

  // Loop1: through all particles
  for (int i = 0; i < ps.no_of_particles; ++i) {

    // Clear the present list
    particle[i].Neighbours.clear(); // this function takes time

    for (int j = 0; j < ps.no_of_particles; ++j) {

      // i) Calculate separation dx and dy
      double dx = particle[i].x - particle[j].x;
      double dy = particle[i].y - particle[j].y;

      // ii) Particle distance
      double d = std::sqrt(dx * dx + dy * dy);

      if (d <= (this->parameters[1]) + (this->NL_skin_depth)) { // rs = d + 2.0
        particle[i].Neighbours.push_back(
            j); // j is a neighbour of i and also vice versa
      }
    }
  }
}

void ParSim::Langevin_Dynamics::Force_NL_PBC_WCA(ParSim::ParticleSystem &ps) {

  Particle *const particle = ps.get_particles();

  // Loop1: through all particles
  for (int i = 0; i < ps.no_of_particles; ++i) {

    // Store previous step forces
    particle[i].force_radial_prev[0] = particle[i].force_radial[0];
    particle[i].force_radial_prev[1] = particle[i].force_radial[1];
    particle[i].force_tangential_prev[0] = particle[i].force_tangential[0];
    particle[i].force_tangential_prev[1] = particle[i].force_tangential[1];
    particle[i].torque_prev = particle[i].torque;

    // Reset the current forces and pressures
    particle[i].force_radial[0] = 0;
    particle[i].force_radial[1] = 0;
    particle[i].force_tangential[0] = 0;
    particle[i].force_tangential[1] = 0;
    particle[i].torque = 0;

    particle[i].sigma[0] = 0.0;
    particle[i].sigma[1] = 0.0;
    particle[i].sigma[2] = 0.0;
    particle[i].sigma[3] = 0.0;

    // Radial and Tangential Unary forces

    particle[i].force_radial[0] += 0.0;

    particle[i].force_radial[1] += 0.0;

    particle[i].force_tangential[0] += 0.0;

    particle[i].force_tangential[1] += 0.0;

    particle[i].torque += 0.0;

    // Force calculation for Neighbour list particles
    std::list<int>::iterator j;

    for (j = particle[i].Neighbours.begin(); j != particle[i].Neighbours.end();
         j++) {

      // Distance from the nearest image of jth particle

      // i)calculate min separation dx and dy
      double dx = particle[i].x - particle[*j].x;

      // if (dx < (-ps.L / 2)) {
      //   dx += ps.L;
      // } else if (dx > (ps.L / 2)) {
      //   dx -= ps.L;
      // }

      double dy = particle[i].y - particle[*j].y;

      // if (dy < (-ps.L / 2)) {
      //   dy += ps.L;
      // } else if (dy > (ps.L / 2)) {
      //   dy -= ps.L;
      // }
      // ii) Nearest image distance
      double d = std::sqrt(dx * dx + dy * dy);

      // U

      if (d <= (this->parameters[1])) {

        // relative velocities calculation (ith particle at rest)
        double vrel_x = 0.0; // particle[*j].vx - particle[i].vx;
        double vrel_y = 0.0; // particle[*j].vy - particle[i].vy;

        // mag. of radial force  (positive)
        double r_by_sigma = d / this->WCA_parameters[1]; // r/sigma in WCA
        double r_by_sigma5 =
            r_by_sigma * r_by_sigma * r_by_sigma * r_by_sigma * r_by_sigma;

        double N =
            4 * (this->WCA_parameters[0] / this->WCA_parameters[1]) *
            (1.0 /
             (r_by_sigma5 * r_by_sigma5 * r_by_sigma * r_by_sigma * r_by_sigma +
              (this->parameters[6]))) *
            ((24.0 / (r_by_sigma5 * r_by_sigma5 * r_by_sigma * r_by_sigma +
                      (this->parameters[6]))) -
             (12.0)); // WCA force magnitude (positive)

        // components of unit separation vector r
        double r1 = (particle[i].x - particle[*j].x) /
                    (d + (this->parameters[6])); // rx
        double r2 = (particle[i].y - particle[*j].y) /
                    (d + (this->parameters[6])); // ry

        double vreldotr = (r1 * vrel_x + r2 * vrel_y);
        // Radial interaction force with jth particle
        double frx = N * r1 + (this->parameters[5]) * vreldotr * r1;
        double fry = N * r2 + (this->parameters[5]) * vreldotr * r2;

        particle[i].force_radial[0] += frx;
        particle[i].force_radial[1] += fry;

        // Tangential friction force with particle j
        double omega_sum = (particle[i].omega + particle[*j].omega);

        double ftx = -(this->parameters[4]) * 1 * (omega_sum / 2.0) * r2 +
                     (this->parameters[5]) * (vrel_x - vreldotr * r1);

        double fty = -(this->parameters[4]) * 1 * (omega_sum / 2.0) * (-r1) +
                     (this->parameters[5]) * (vrel_y - vreldotr * r2);

        particle[i].force_tangential[0] += ftx;
        particle[i].force_tangential[1] += fty;

        // torque on particle

        if (omega_sum != 0) {
          particle[i].torque +=
              -(this->parameters[4]) * 1 * (omega_sum / 2.0) * d;
        }

        if (omega_sum == 0) {
          particle[i].torque +=
              -(this->parameters[4]) * 1 * (omega_sum / 2.0) * d;
        }

        // pressure term calculation

        particle[i].sigma[0] += (frx + ftx) * (particle[i].x - particle[*j].x);
        particle[i].sigma[1] += (frx + ftx) * (particle[i].y - particle[*j].y);
        particle[i].sigma[2] += (fry + fty) * (particle[i].x - particle[*j].x);
        particle[i].sigma[3] += (fry + fty) * (particle[i].y - particle[*j].y);
      }
    }
  }
}

// 3) Integrators---------------

void ParSim::Langevin_Dynamics::LM_Intergrator1(ParSim::ParticleSystem &parsym,
                                                int step) {

  Particle *const particle = parsym.get_particles();

  double Fx = 0.0;
  double Fy = 0.0;
  double Tau = 0.0;
  for (int i = 0; i < parsym.no_of_particles; ++i) {

    // reset forces for this particle
    Fx = 0.0;
    Fy = 0.0;
    Tau = 0.0;

    // Total force and torque on this ith particle
    Fx = particle[i].force_radial[0] + particle[i].force_tangential[0];
    Fy = particle[i].force_radial[1] + particle[i].force_tangential[1]; //

    Tau = parsym.particle_array[i].torque;

    // Save present attributes
    particle[i].position_prev[0] = particle[i].x;
    particle[i].position_prev[1] = particle[i].y;

    particle[i].velocity_prev[0] = particle[i].vx;
    particle[i].velocity_prev[1] = particle[i].vy;

    particle[i].alpha_prev = particle[i].alpha;
    particle[i].omega_prev = particle[i].omega;

    // update the attributes upto midpoint

    // i) Update v(t) to v(t+(dt/2))

    particle[i].vx +=
        (Fx / (this->parameters[2])) * (this->parameters[8]) / 2; // v'
    particle[i].vy += (Fy / (this->parameters[2])) * (this->parameters[8]) / 2;

    particle[i].omega +=
        (Tau / (this->parameters[13])) * (this->parameters[8]) / 2;

    // ii) Update x(t) to x(t+(dt/2)) using v(t+(dt/2))
    particle[i].x += particle[i].vx * (this->parameters[8]) / 2; // x'
    particle[i].y += particle[i].vy * (this->parameters[8]) / 2;
    particle[i].alpha += particle[i].omega * (this->parameters[8]) / 2;

    // iii) Add Fluctuation and Dissipation changing v(t+(dt/2)) -->
    // v'(t+(dt/2))

    std::normal_distribution<double> distribution(
        0.0, 1.0); // G: noise with zero mean and unit variance

    double EXP =
        exp(-this->FD_parameters[0] * (this->parameters[2]) *
            (this->parameters[8])); // exp(-gamma_t*dt) , gamma_t = m*jeta_t
    double EXP_r =
        exp(-this->FD_parameters[1] * (this->parameters[13]) *
            (this->parameters[8])); // exp(-gamma_r*dt) , gamma_r = I*jeta_r
    double coefficient_G = (sqrt(1 - EXP * EXP) *
                            sqrt((this->parameters[2]) *
                                 this->FD_parameters[2])); // coefficient of G
    double coefficient_G_r =
        (sqrt(1 - EXP_r * EXP_r) *
         sqrt((this->parameters[13]) *
              this->FD_parameters[2])); // coefficient of G_r
    particle[i].vx =
        EXP * particle[i].vx + coefficient_G * distribution(this->mt);
    particle[i].vy =
        EXP * particle[i].vy + coefficient_G * distribution(this->mt);
    particle[i].omega =
        EXP_r * particle[i].omega + coefficient_G_r * distribution(this->mt);

    /// iv) Update  x(t+(dt/2)) to x(t + dt) using v'(t+(dt/2))
    particle[i].x += particle[i].vx * (this->parameters[8]) / 2; // x'
    particle[i].y += particle[i].vy * (this->parameters[8]) / 2;
    particle[i].theta += particle[i].omega * (this->parameters[8]) / 2;

    // For updating v'(t+(dt/2)) to v(t+dt) we need to calculate forces (that
    // are position dependent) using x(t + dt); we will now caculate forces
    // again
  }
}

void ParSim::Langevin_Dynamics::LM_Intergrator2(ParSim::ParticleSystem &parsym,
                                                int step) {

  Particle *const particle = parsym.get_particles();

  double Fx = 0.0;
  double Fy = 0.0;
  double Tau = 0.0;
  for (int i = 0; i < parsym.no_of_particles; ++i) {

    // reset forces for this particle
    Fx = 0.0;
    Fy = 0.0;
    Tau = 0.0;

    // Total force and torque on this ith particle
    Fx = particle[i].force_radial[0] + particle[i].force_tangential[0];
    Fy = particle[i].force_radial[1] + particle[i].force_tangential[1]; //

    Tau = parsym.particle_array[i].torque;

    // Save present attributes
    particle[i].position_prev[0] = particle[i].x;
    particle[i].position_prev[1] = particle[i].y;

    particle[i].velocity_prev[0] = particle[i].vx;
    particle[i].velocity_prev[1] = particle[i].vy;

    particle[i].alpha_prev = particle[i].alpha;
    particle[i].omega_prev = particle[i].omega;

    // update the attributes upto midpoint

    // v) Update v'(t+(dt/2)) to v(t+dt) using new forces f(t+dt)

    particle[i].vx += (Fx / (this->parameters[2])) * (this->parameters[8]) / 2;
    particle[i].vy += (Fy / (this->parameters[2])) * (this->parameters[8]) / 2;

    particle[i].omega +=
        (Tau / (this->parameters[13])) * (this->parameters[8]) / 2;

    // We don't need to update anything else:
  }
}

// 4) Evolvers -----------

void ParSim::Langevin_Dynamics::evolve_system_LM_NL(ParticleSystem &parsym,
                                                    int step) {

  // Update Neighbour list after every NL_tau steps
  if (step == 0 || step % this->NL_tau == 0) {
    Neighbours_search(parsym);
  }
  // LM integration: (x(t),v(t)) -- > (x(t + dt),v(t + dt))

  // i) Calculate force f(t) using (x(t)) --
  Force_NL_PBC_WCA(parsym);

  // ii) Update (x(t), v(t))  to (x(t+dt) , v'(t+(dt/2))
  LM_Intergrator1(parsym, step);

  // iii) Again calculate force f(t+dt) from x((t+dt))
  Force_NL_PBC_WCA(parsym);

  //  iv) Update v'(t+(dt/2)) to v(t+dt) --------
  LM_Intergrator2(parsym, step);
}
/*Conserved Quantities*/
std::vector<double>
ParSim::Langevin_Dynamics::EnergyMomentum(ParticleSystem &parsym) {
  std::vector<double> p{0, 0, 0, 0}; // kinetic, rotational, momenta
  double m = this->parameters[2];
  double I = m;
  for (int i = 0; i < parsym.no_of_particles; ++i) {
    p[0] += 0.5 * m *
            (pow(parsym.particle_array[i].vx, 2) +
             pow(parsym.particle_array[i].vy, 2)); // KE
    p[1] +=
        0.5 * I * (pow(parsym.particle_array[i].omega, 2)); // Rotational K.E
    p[2] += m * parsym.particle_array[i].vx;                // px
    p[3] += m * parsym.particle_array[i].vy;                // py
  }
  return p;
}