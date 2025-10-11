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

## **Procedures Implemented**

### Maximum Likelihood Estimation of Stochastic Processes
**Ornstein-Uhlenbeck (Mean Reverting) Model**
This assumes that the Ornstein-Uhlenbeck process takes the form
$$\mathrm{d}X(t) = \alpha(\mu-x(t))\mathrm{d}t + \sigma \mathrm{d}W(t), $$
with $\alpha$ the drift term and $\sigma$ diffusion term. The term $\mathrm{d}W(t)$ comprises the Brownian motion.


**General Stochastic Differential Equation (SDE) Model**
This assumes that the SDE process takes the form
$$\mathrm{d}X(t) = \mu(X,t)\mathrm{d}t + \sigma(X,t) \mathrm{d}W(t), $$
with $\mu(X,t)$ the drift term and $\sigma(X,t)$ diffusion term. The term $\mathrm{d}W(t)$ comprises the Brownian motion.

### Kinetic Components Analysis
A kalman-filter model adapted to capture the instantaneous first and second derivatives of a stochastic process. It is based on the Kinetic Components Analysis paper by López de Prado (2016) with the ability to produce estimates given a generic stochastic process characterised by a stochastic differential equation.

I have done my best to implement what I have read and understood from various public sources detailing the KCA method and have added the default fit of a stochastic differential equation as a slight modification. This should reduce the error of the estimates to some extent.

**Citation:** López de Prado, Marcos and López de Prado, Marcos and Rebonato, Riccardo, Kinetic Component Analysis (June 5, 2016). Journal of Investing, Vol. 25, No. 3, 2016.

### Optimal Trading Levels for Mean-Reverting Assets
The library is modularised and each component serves a specific purpose in using SDEs in mathematical analysis of trading and markets.

Code modules (i.e. `hitting_times`) provide implementations of analysis outlined in the book 'Optimal Mean Reversion Trading: Mathematical Analysis and Practical Applications' by Tim Leung and Xin Li (2015).

**Citation:** Leung, T., & Li, X. (2015). Optimal Mean Reversion Trading: Mathematical Analysis and Practical Applications.

## Examples
I am still working on providing a clean and concise API. Only certain services are `PUBLIC` to enforce correct usage of the library's core modules:

Outputs an optimal exit level given the parameters of an Ornstein-Uhlenbeck process:
```
#include <iostream>
#include "entrypoints/optimal_trading_levels.h"

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

Outputs the maximum likelihood estimate of an Ornstein-Uhlenbeck process parameters given a vector of data:
```
#include <iostream>
#include "entrypoints/ou_model.h"

const float tolerance = 1e-5;
// Vector of data.
const std::vector<double> test_vec{0.5, 0.25, 0.5, 0.75, 1.5, 0.5};

// Likelihood estimates.
const std::unordered_map<std::string, const double> likelihood =
    ornsteinUhlenbeckMaximumLikelihood(test_vec);

// Approx 0.58333333.
std::cout << likelihood.at("mu");
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

**_The procedures mentioned above can be accessed in a complete and safe way through the `entrypoints` module. This module provides a simple set of functions orchestrating the underlying models and calculations._**

These can be imported in your application as shown below.

It is possible to use this application as a shared library within another application should it be
correctly linked. The `.so` file created in the `src` directory is the most efficient way to use this in another project. The `../install/include` directory is also required. The linker must be able to locate these files to build the application.

To access in the `main` module of another application a `CMakeLists.txt` can be produced
and `#include <stochastic_models/sde/ornstein_uhlenbeck.h>` for example to access the `OrnsteinUhlenbeckModel` interface.
```
cmake_minimum_required(VERSION 4.0)
project(another-application VERSION 0.1.0 DESCRIPTION "A project with external library" LANGUAGES CXX)

# Add the subdirectory to build the library as a subproject to follow how stochastic-models is built.
add_subdirectory(stochastic-models)

# Require GSL for linking (ensure it's installed prior on your system).
find_package(GSL REQUIRED)

# Create the executable and link it to the subdirectory project.
add_executable(${PROJECT_NAME} main.cpp)
target_include_directories(${PROJECT_NAME} PRIVATE "${CMAKE_SOURCE_DIR}/stochastic-models/include")
target_link_libraries(${PROJECT_NAME} PRIVATE stochastic_models GSL::gsl)

# At time of writing I used the C++20 standard.
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_20)

```

Your `main.cpp` could look like:
```
#include <stochastic_models/entrypoints/ou_model.h>

#include <iostream>
int main() {
    const std::vector<double> vec = simulateOrnsteinUhlenbeck(0.5, 0.01, 0.067, 0, 50, 1);

    const std::unordered_map<std::string, const double> m =
        ornsteinUhlenbeckMaximumLikelihood(vec);

    for (const auto &elem : m) {
        std::cout << elem.first << " " << elem.second << std::endl;
    }
    return 0;
}
```

Then run 
```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```
