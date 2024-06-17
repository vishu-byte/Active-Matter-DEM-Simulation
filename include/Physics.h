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


/**
 * @file Physics.h
 *
 * @brief Contains classes Physics and Langevin_Dyanmics
 *
 * @author Vishu Saini
 * Contact: vishusaini220301@gmail.com
 *
 */


#ifndef PARTICLE_SIMULATION_PHYSICS_H
#define PARTICLE_SIMULATION_PHYSICS_H

#include "ParticleSystem.h"
#include <math.h>
#include <random>

#include <fstream>
#include <iostream> 
#include <vector>

namespace ParSim {
/**
 * For managing physics part of the simulation
 *
 * Evolution of a particle system in time is handled
 * using methods of this class. Evolution equations are
 * Newton's dyanmical equations.
 *
 */
class Physics { 
public:
  double parameters[14];
  double WCA_parameters[2];
  double circular_wall_parameters[2];
  double NL_skin_depth;
  int NL_tau;
  // std::random_device generator;    //this is used in constructor

  std::mt19937 mt; // all objects should share only one copy
  // std::normal_distribution<double> distribution;

  Physics(); // constructor

  virtual ~Physics(){}; // destructor
public:
  /*Force linker + integrators-- */
  // 1) Particle-Particle force linkers
  void Force_PP(ParticleSystem &); // Base force linker
  void Force_PP2(ParticleSystem &);
  void Force_PP_PBC(ParticleSystem &);
  void Force_PP_PBC2(ParticleSystem &); // force linker with PBC
  void Force_PP_CRB(ParticleSystem &);
  void Force_PP_CRB2(ParticleSystem &); // force linker with CRB
  void Force_PP_CRB_WCA(ParticleSystem &);
  void Force_PP_CRB_WCA2(ParticleSystem &); // force linker with WCA + CRB



  // 2) Neighbour list force linkers
  void Neighbours_search_PBC(ParticleSystem &); // Neighbour search method
  void Neighbours_search(ParticleSystem &); // Neighbour search method
  void Force_NL_PBC(ParticleSystem &);      // Force linker with PBC
  void Force_NL_CRB_WCA(ParticleSystem &);      // Force linker with WCA + CRB




  // 3)Integrators
  void Euler_Integrator(Particle &, int);
  void Vel_Verlet_Integrator(Particle &, int);
  void ERM_Integrator1(Particle &, int);
  void ERM_Integrator2(Particle &, double, int);
  void ERM_Integrator1_sys(ParticleSystem &, int); // updated intergrators
  void ERM_Integrator2_sys(ParticleSystem &, int); // updated intergrator

  // 4)General evolvers
  void evolve_system(ParticleSystem &, int); // takes a particle system and
                                             // moves it forward in time
  void evolve_system_ERM(ParticleSystem &, int);
  void evolve_system_ERM_NL(ParticleSystem &, int); // Neighbour list evolver

  /*Conserved quantities*/
  std::vector<double> EnergyMomentum(ParticleSystem &);
};

/**
 * For systems whose evolution is governed by Langevin equations
 *
 *
 */
class Langevin_Dynamics{
public:
  double parameters[14];
  double FD_parameters[4];
  double WCA_parameters[2];
  double circular_wall_parameters[2];
  double NL_skin_depth;
  int NL_tau;
  // std::random_device generator;    //this is used in constructor

  std::mt19937 mt; // all objects should share only one copy
  // std::normal_distribution<double> distribution;

  Langevin_Dynamics(); // constructor

  virtual ~Langevin_Dynamics(){}; // destructor
public:
  /*Force linker + integrators-- */
  // 1) Particle-Particle force linkers
  //---//

  // 2) Neighbour list force linkers
  void Neighbours_search_PBC(ParticleSystem &); // Neighbour search method
  void Neighbours_search(ParticleSystem &); // Neighbour search method
  void Force_NL_PBC_WCA(ParticleSystem &);      // Force linker with PBC
  //void Force_NL_CRB_WCA(ParticleSystem &);      // Force linker with WCA + CRB




  // 3)Integrators
  void Euler_Integrator(Particle &, int);
  void LM_Intergrator1(ParticleSystem &, int);
  void LM_Intergrator2(ParticleSystem &, int);
  

  // 4)General evolvers

  void evolve_system_LM_NL(ParticleSystem &, int); // Neighbour list evolver

  /*Conserved quantities*/
  std::vector<double> EnergyMomentum(ParticleSystem &);

};

} // namespace ParSim

#endif // PARTICLE_SIMULATION_PHYSICS_H
