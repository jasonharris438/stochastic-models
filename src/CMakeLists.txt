add_library(stochastic_models
adapters.cpp
stochastic_model.cpp
ornstein_uhlenbeck.cpp
core.cpp
gaussian.cpp
gsl_errors.cpp
filter.cpp
states.cpp
states_exceptions.cpp
linalg.cpp
likelihood.cpp
integration.cpp
solvers.cpp
differentiation.cpp
general_linear.cpp
general_linear_online.cpp
general_linear_likelihood.cpp
ornstein_uhlenbeck_likelihood.cpp
helpers.cpp
hitting_time_ornstein_uhlenbeck.cpp
hitting_time_density.cpp
optimal_mean_reversion.cpp
exponential_mean_reversion.cpp
trading_levels.cpp
type_conversion.cpp
interface.cpp)

target_include_directories(stochastic_models
PRIVATE
# where the library itself will look for its internal headers
"${PROJECT_SOURCE_DIR}/include"
PUBLIC
# where top-level project will look for the library's public headers
$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/stochastic_models/interface>
$<BUILD_INTERFACE:${boost_SOURCE_DIR}>
# where external projects will look for the library's public headers
$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/interface>
)

target_link_libraries(
    stochastic_models
    GSL::gsl
    nlohmann_json::nlohmann_json
)
