cmake_minimum_required(VERSION 3.20)

# set the project name and version.
project(StochasticModels
DESCRIPTION "A package implementing fast statistical tools for stochastic calculus."
VERSION 0.0.1
)

if (NOT CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    message(STATUS "This project has a top-level one called [${CMAKE_PROJECT_NAME}]")
else()
    message(STATUS "This project is a top-level one")
endif()

add_library(${PROJECT_NAME} STATIC .)
set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)

# add subdirectory stochastic_models.
add_subdirectory(stochastic_models)

# specify the C++ standard.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# add the executable.
add_executable(run main.cpp)

# link libraries to executable.
target_link_libraries(run PUBLIC stochastic_models)
