# Stochastic Models
A C++ library with statistical tools for stochastic models.

## **Models Implemented**
### Ornstein-Uhlenbeck (Mean Reverting) Model
This assumes that the Ornstein-Uhlenbeck process takes the form
$$\mathrm{d}X(t) = \alpha(\mu-x(t))\mathrm{d}t + \sigma \mathrm{d}W(t), $$
with $\alpha$ the drift term and $\sigma$ diffusion term. The term $\mathrm{d}W(t)$ comprises the Brownian motion.

#### Implementation
I am still working on providing a clean and concise API. Useage of the parameter search in the meantime can take the form (see also [tests](stochastic_models/tests/ornstein_uhlenbeck_likelihood_test.cpp)):
```
#include <iostream>

#include "likelihood/ornstein_uhlenbeck_likelihood.h"

std::vector<double> vec; // Initialise as you wish.
// Create instance of OrnsteinUhlenbeckLikelihood model.
OrnsteinUhlenbeckLikelihood *likelihood = new OrnsteinUhlenbeckLikelihood();
/// Calculate parameters.
const std::unordered_map<std::string, const double> m =
    likelihood->calculate(vec);

/// Show contents for demonstration.
for (const auto &elem : m) {
    std::cout << elem.first << " " << elem.second << std::endl;
}
```
Useage of the simulation tool that uses the Euler–Maruyama method for the approximate numerical solution:
```
#include <iostream>

#include "sde/ornstein_uhlenbeck.h"

// New model instance with mu, alpha, sigma parameters.
OrnsteinUhlenbeckModel *model =
    new OrnsteinUhlenbeckModel(0.5, 0.01, 0.0067);
// Simulate 5,000 instances with a start point of 0.5.
const std::vector<double> vec = (*model).Simulate(0.5, 5000);
```
## **To Build**
(This is a component that I am still learning as I go with, so any suggestions on improvements or corrections are welcomed :pray:)
There is one subdirectory to be linked ->> `stochastic_models` which contains all header files for the library.

To build and install to an `install` directory, run the following:
```
$ cmake .
$ cmake --build . --target install
```
## **To Use**
It is possible to use this library statically within another application should it be
correctly linked (I've tested this).

To access in the `main` module of another application a `CMakeLists.txt` can be produced
and `#include <StochasticModels/sde/ornstein_uhlenbeck.h>` for example to access the `OrnsteinUhlenbeckModel` interface.
```
# Minimum version 3.26.
cmake_minimum_required(VERSION 3.26)
project("another-application" VERSION 0.1.0 DESCRIPTION "A project with external library")

# provide the library installation folder, so CMake could find its config
set(CMAKE_PREFIX_PATH "/install")
# the rest will be taken care of by CMake
find_package(StochasticModels CONFIG REQUIRED)

# it is an application
add_executable(${PROJECT_NAME} main.cpp)
target_include_directories(${PROJECT_NAME} PUBLIC include/StochasticModels)

target_sources(${PROJECT_NAME}
    PRIVATE
        main.cpp
)

# linking to the library
target_link_libraries(${PROJECT_NAME} PRIVATE StochasticModels)
```
