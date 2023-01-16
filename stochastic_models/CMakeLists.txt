add_library(stochastic_models
src/stochastic_model.cpp
src/ornstein_uhlenbeck.cpp
src/core.cpp
src/gaussian.cpp
src/likelihood.cpp
src/ornstein_uhlenbeck_likelihood.cpp)

target_include_directories(stochastic_models PUBLIC include)
