//
// Created by Vishu Saini on 31/08/23
//

#ifndef PARTICLE_SIMULATION_SCREEN_H
#define PARTICLE_SIMULATION_SCREEN_H

#include "ParticleSystem.h"
#include <SDL2/SDL.h>
#include <iostream>

namespace ParSim {

class Screen {
private:
  SDL_Window *m_window{nullptr};
  SDL_Renderer *m_renderer{nullptr};
  SDL_Event m_event;
  int SCREEN_WIDTH{1000}; // parameters of window creation
  int SCREEN_HEIGHT{800};

public:
  Screen();
  virtual ~Screen();
  void update();
  void draw_particlesystem(ParticleSystem &);
  bool quit_program();

private:
  void init_SDL();
  void init_window();
  void init_renderer();
  void Draw_rect(int x, int y);
};

} // namespace ParSim

#endif // PARTICLE_SIMULATION_SCREEN_H
