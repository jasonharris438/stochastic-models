add_library(stochastic_models
adapters.cpp
core.cpp
differentiation.cpp
entrypoint_general_sde.cpp
entrypoint_kca_filter.cpp
entrypoint_optimal_trading_levels.cpp
entrypoint_ou_model.cpp
exponential_mean_reversion.cpp
gaussian.cpp
general_linear.cpp
general_linear_likelihood.cpp
general_linear_online.cpp
gsl_errors.cpp
helpers.cpp
hitting_time_density.cpp
hitting_time_ornstein_uhlenbeck.cpp
integration.cpp
kca.cpp
linalg.cpp
ornstein_uhlenbeck_likelihood.cpp
ornstein_uhlenbeck_online.cpp
optimal_mean_reversion.cpp
ornstein_uhlenbeck.cpp
solvers.cpp
states.cpp
states_exceptions.cpp
stochastic_model.cpp
trading_levels.cpp
trading_levels_exponential.cpp
trading_levels_params.cpp
type_conversion.cpp
)

target_include_directories(stochastic_models
PRIVATE
# where the library itself will look for its internal headers
"${PROJECT_SOURCE_DIR}/include"
PUBLIC
# where top-level project will look for the library's public headers
$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/stochastic_models/entrypoints>
$<BUILD_INTERFACE:${boost_SOURCE_DIR}>
# where external projects will look for the library's public headers
$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/entrypoints>
)

target_link_libraries(
    stochastic_models
    GSL::gsl
    nlohmann_json::nlohmann_json
)
