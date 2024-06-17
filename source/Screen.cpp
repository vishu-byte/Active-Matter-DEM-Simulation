//
// Created by Vishu Saini on 31/08/23
//

#include "../include/Screen.h"
#include <stdlib.h>
namespace ParSim {

Screen::Screen() : m_window{nullptr}, m_renderer{nullptr} {
  // Initialize all required SDL functionality and create SDL objects.
  init_SDL();
  init_window();
  init_renderer();
}

Screen::~Screen() {
  // Destroy all SDL related objects.
  SDL_DestroyRenderer(m_renderer);
  SDL_DestroyWindow(m_window);
  SDL_Quit();
}

void Screen::init_SDL() {
  // Initialize SDL library. Returns 0 if successful.
  // Must be called before using any SDL functionality.
  if (SDL_Init(SDL_INIT_VIDEO)) { // initialized for video
    std::cout << "SDL_Init Error: \n" << SDL_GetError() << std::endl;
    exit(1);
  }
}

void Screen::init_window() {
  // Creates an SDL window struct. Returns NULL on failure.
  m_window = SDL_CreateWindow("Particle Simulation",   // Window title
                              SDL_WINDOWPOS_UNDEFINED, // Initial x position
                              SDL_WINDOWPOS_UNDEFINED, // Initial y position
                              SCREEN_WIDTH,            // Width, in pixels
                              SCREEN_HEIGHT,           // Height, in pixels
                              SDL_WINDOW_ALLOW_HIGHDPI // Increases resolution
  );

  // Check if window creation failed.
  if (!m_window) {
    std::cout << "SDL_CreateWindow Error: \n" << SDL_GetError() << std::endl;
    exit(2);
  }
}

void Screen::init_renderer() {
  // Creates an SDL renderer struct: a rendering context relative to the
  // state/information of the associated window. Textures to be rendered pass
  // through this renderer and when done rendering are sent to the window to be
  // displayed.
  m_renderer = SDL_CreateRenderer(
      m_window, // Window associated with the renderer.
      -1, // Index of rendering driver to initialize, -1 to use first available.
      SDL_RENDERER_PRESENTVSYNC | // Synchronize rendering with window refresh
                                  // rate.
          SDL_RENDERER_ACCELERATED); // Allow renderer to use hardware
                                     // acceleration.

  // Check if renderer creation failed.
  if (!m_renderer) {
    std::cout << "SDL_CreateRenderer Error: \n" << SDL_GetError() << std::endl;
    SDL_DestroyWindow(m_window);
    exit(3);
  }
}

void Screen::draw_particlesystem(ParticleSystem &parsym) {
  // Load parsym particles
  Particle *const p_particles =
      parsym.get_particles(); // returns particle array

  SDL_SetRenderDrawColor(m_renderer, 0, 0, 0, 255); // dip the brush in black
                                                    // color
  SDL_RenderClear(m_renderer);                      // clear the screen to black

  // Link parsym particles to a SDL draw object
  for (int i = 0; i < parsym.no_of_particles; ++i) {
    ParSim::Particle particle =
        p_particles[i]; // access the ith particle of particle array

    int x = static_cast<int>(particle.x + (SCREEN_WIDTH / 2));
    int y = static_cast<int>(-particle.y + (SCREEN_HEIGHT / 2));    //canvas' y axis has to be inverted first

    // // working fine till now
    // std::cout << "Particle " << i << ": " << x << std::endl;
    // std::cout << "Particle " << i << ": " << y << std::endl;

    Draw_rect(x, y);
  }

  // Loads the renderer to the SDL window.
  SDL_RenderPresent(m_renderer);
}

void Screen::Draw_rect(int x, int y) {

  // Ignore particles with coordinates outside the boundaries of the SDL window
  if (x < 0 || x >= ParSim::Screen::SCREEN_WIDTH || y < 0 ||
      y >= ParSim::Screen::SCREEN_HEIGHT)
    return;

  // Draw a rectangle at each particle position----

  SDL_SetRenderDrawColor(m_renderer, 255, 255, 255, 255);

  // Remember the SDL window uses the coordinate system with origin at top left
  // corner
  SDL_Rect rect;

  rect.w = 10;
  rect.h = 10;
  rect.y = y;
  rect.x = x; // position of rectangle

  SDL_RenderDrawRect(m_renderer, &rect);
}

bool Screen::quit_program() {
  // Check for SDL events. If the window is closed, quit the program.
  while (SDL_PollEvent(&m_event)) {
    if (m_event.type == SDL_QUIT)
      return true;
  }
  return false;
}

} // namespace ParSim
