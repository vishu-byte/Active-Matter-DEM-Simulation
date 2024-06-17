cmake_minimum_required(VERSION 3.6)
project(Particles_Simulation)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")


set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)


set(SOURCE_FILES source/main.cpp source/ParticleSystem.cpp include/ParticleSystem.h source/Physics.cpp include/Physics.h)
add_executable(Particles_Simulation ${SOURCE_FILES})

