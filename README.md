# Stochastic Models
A C++ library with statistical tools for stochastic models.

## **Models Implemented**
### Ornstein-Uhlenbeck (Mean Reverting) Model
This assumes that the Ornstein--Uhlenbeck process takes the form
$$\mathrm{d}X(t) = -\alpha(\mu-x(t))\mathrm{d}t + \sigma \mathrm{d}W(t), $$
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
