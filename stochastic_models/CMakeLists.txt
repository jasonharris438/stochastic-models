add_library(stochastic_models
src/stochastic_model.cpp
src/ornstein_uhlenbeck.cpp
src/core.cpp
src/gaussian.cpp
src/likelihood.cpp
src/integration.cpp
src/ornstein_uhlenbeck_likelihood.cpp
src/helpers.cpp
src/hitting_time_density.cpp)

target_include_directories(stochastic_models PUBLIC include)
