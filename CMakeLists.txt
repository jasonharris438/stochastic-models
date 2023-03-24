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

# Add GSL package.
find_package(GSL REQUIRED)

# Add googletest installation.
include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        b796f7d44681514f58a683a3a71ff17c94edb0c1 # release-1.13.0
)
FetchContent_MakeAvailable(googletest)

# add the executable.
add_executable(stochastic-models main.cpp)

# link libraries to executable.
target_link_libraries(stochastic-models PUBLIC stochastic_models GSL::gsl)

enable_testing()

add_executable(
    tests
    stochastic_models/tests/example_test.cpp
    stochastic_models/tests/ornstein_uhlenbeck_likelihood_test.cpp
)
target_link_libraries(
    tests
    GTest::gtest_main
    stochastic_models
)

include(GoogleTest)
gtest_discover_tests(tests)
