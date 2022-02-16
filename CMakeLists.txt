cmake_minimum_required(VERSION 3.20)

# set the project name and version
project(Stochastic-Models VERSION 0.0.1)

add_library(StochasticModels STATIC lib/src/ornstein_uhlenbeck.cpp)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(APP_EXE Stochastic-Models)

# add the executable
add_executable(${APP_EXE} main.cpp)

TARGET_LINK_LIBRARIES(${APP_EXE} StochasticModels)
