# Stochastic Models
A C++ shared library implementing the core calculations involved in stochastic models.

## **About This Project**

What started as a personal exploration into quantitative finance has evolved into a more detailed C++ library for stochastic modeling and analysis. This was just a hobby project to challenge myself to put mathematical concepts from academic papers and textbooks into code, but it has since become a more sophisticated set of tools covering stochastic differential equations, Kalman filtering, maximum likelihood estimation, and optimal trading strategies.

The codebase now features a more modular architecture, some numerical utilities, and implementations of various concepts from publicly available research including Kinetic Components Analysis and 'optimal mean reversion trading' methodologies. While it is relatively functionally and somewhat mathematically sound, there is a LOT of room for improvement in code organization, API design, performance optimization, and documentation. It will help me to keep learning too..

**I welcome contributions, feedback, and collaboration!** Whether you're interested in:
- Code reviews and architectural improvements
- Additional mathematical models or algorithms
- Performance optimizations
- Better documentation and examples
- Testing and validation
- Packaging for distribution

Please feel free to open issues, submit pull requests, or reach out with suggestions. My goal is to develop this into a robust, well-tested library that could be valuable to the broader quantitative finance and stochastic modeling community.

## **Models Implemented**
### Ornstein-Uhlenbeck (Mean Reverting) Model
This assumes that the Ornstein-Uhlenbeck process takes the form
$$\mathrm{d}X(t) = \alpha(\mu-x(t))\mathrm{d}t + \sigma \mathrm{d}W(t), $$
with $\alpha$ the drift term and $\sigma$ diffusion term. The term $\mathrm{d}W(t)$ comprises the Brownian motion.

### General Stochastic Differential Equation (SDE) Model
This assumes that the SDE process takes the form
$$\mathrm{d}X(t) = \mu(X,t)\mathrm{d}t + \sigma(X,t) \mathrm{d}W(t), $$
with $\mu(X,t)$ the drift term and $\sigma(X,t)$ diffusion term. The term $\mathrm{d}W(t)$ comprises the Brownian motion.


### Kinetic Components Analysis
A kalman-filter model adapted to capture the instantaneous first and second derivatives of a stochastic process. It is based on the Kinetic Components Analysis paper by López de Prado (2016) with the ability to produce estimates given a generic stochastic process characterised by a stochastic differential equation.

I have done my best to implement what I have read and understood from various public sources detailing the KCA method and 

**Citation:** López de Prado, Marcos and López de Prado, Marcos and Rebonato, Riccardo, Kinetic Component Analysis (June 5, 2016). Journal of Investing, Vol. 25, No. 3, 2016.

### Implementation
The library is modularised and each component serves a specific purpose in using SDEs in mathematical analysis of trading and markets.

Code modules (i.e. `hitting_times`) provide implementations of analysis outlined in the book 'Optimal Mean Reversion Trading: Mathematical Analysis and Practical Applications' by Tim Leung and Xin Li (2015).

**Citation:** Leung, T., & Li, X. (2015). Optimal Mean Reversion Trading: Mathematical Analysis and Practical Applications.

#### Distributions
Core classical statistical distributions used to derive analysis of data consumed by this library.

#### Exceptions
Custom exceptions to indicate specific error cases thrown while performing internal routines.

#### Hitting Times
Logic to calculate internal components of trading hitting time levels from Leung & Li (2015).

#### Kalman Filter
Handles all kinetic components related calculations.

#### Likelihood
Methodology to calculate maximum likelihood parameters of SDE models used in this library.

#### Numeric Utils
Contains generic numeric algorithms and 3rd party library functionality used throughout the package.

#### SDE
Handles estimating moments of and simulating SDE series.

#### Trading Levels
Controllers and types that handle orchestrating calculating optimal trading levels.

## Examples
I am still working on providing a clean and concise API. Only certain services are `PUBLIC` to enforce correct usage of the library's core modules:
```
#include <iostream>
#include "interface/stochastic_models.h"

// Declare and initialize model parameters.
const double alpha = 8;
const double mu = 0.3;
const double sigma = 0.3;
const double c = 0.02;
const double r = 0.05;

// Calculate b*.
const double value = optimalExitLevel(mu, alpha, sigma, r, c);

/// Output value.
std::cout << "Optimal exit level ~0.46683583 (" << value << ")" << std::endl;
```
**TODO:** Useage of the simulation tool that uses the Euler–Maruyama method for the approximate numerical solution:
```
#include <iostream>
#include "interface/stochastic_models.h"

// Declare and initialize model parameters.
const double alpha = 8;
const double mu = 0.3;
const double sigma = 0.3;
const int size = 5000;
const double x0 = 0.5;

// New model instance with mu, alpha, sigma parameters.
const std::vector<double> vec simulations =
    simulateOrnsteinUhlenbeck(x0, size, mu, alpha, sigma);
```
## **To Build**
There is one subdirectory to be linked ->> `stochastic_models` which contains all header files for the library.

To build and install to an `install` directory, run the following:
```
$ cmake .
$ cmake --build . --target install
```

### Dependencies

#### Boost
- Boost is included [using FetchContent](https://github.com/boostorg/cmake?tab=readme-ov-file#using-boost-with-fetchcontent)
- [`find_package`](https://cmake.org/cmake/help/latest/module/FindBoost.html#examples) method not used due to the header-only libraries being accessible without needing to build the project.

## **To Use**
It is possible to use this application as a shared library within another application should it be
correctly linked. The `.so` file created in the `src` directory is the most efficient way to use this in another project. The `../install/include` directory is also required. The linker must be able to locate these files to build the application.

To access in the `main` module of another application a `CMakeLists.txt` can be produced
and `#include <stochastic_models/sde/ornstein_uhlenbeck.h>` for example to access the `OrnsteinUhlenbeckModel` interface.
```
# Minimum version 3.30.
cmake_minimum_required(VERSION 3.30)
project("another-application" VERSION 0.1.0 DESCRIPTION "A project with external library")

# provide the library installation folder, so CMake could find its config
set(CMAKE_PREFIX_PATH "./install")
# the rest will be taken care of by CMake
find_package(stochastic_models CONFIG REQUIRED)
find_package(GSL REQUIRED) # Add GSL package.

# it is an application
add_executable(${PROJECT_NAME} main.cpp)
target_include_directories(${PROJECT_NAME} PUBLIC include/stochastic_models)

target_sources(${PROJECT_NAME}
    PRIVATE
        main.cpp
)

# linking to the library
target_link_libraries(${PROJECT_NAME} PRIVATE stochastic_models GSL::gsl)
```
Your `main.cpp` could look like:
```
#include <stochastic_models/likelihood/ornstein_uhlenbeck_likelihood.h>
#include <stochastic_models/sde/ornstein_uhlenbeck.h>

#include <iostream>
int main() {
    OrnsteinUhlenbeckModel *model =
        new OrnsteinUhlenbeckModel(0.5, 0.01, 0.0067);
    const std::vector<double> vec = (*model).Simulate(0.5, 5000);

    OrnsteinUhlenbeckLikelihood *likelihood = new OrnsteinUhlenbeckLikelihood();
    const std::unordered_map<std::string, const double> m =
        likelihood->calculate(vec);
    for (const auto &elem : m) {
        std::cout << elem.first << " " << elem.second << std::endl;
    }
    delete model;
    return 0;
}
```
