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

# EXPORT_NAME pairs with the StochasticModels:: namespace, so the exported
# target is StochasticModels::StochasticModels.
set_target_properties(stochastic_models PROPERTIES
    EXPORT_NAME StochasticModels
    VERSION ${PROJECT_VERSION}
    SOVERSION 0
)

# In the export, so consumers on an older standard get one clear diagnostic.
target_compile_features(stochastic_models PUBLIC cxx_std_20)

# Include root on both sides: consumers write
# #include <stochastic_models/entrypoints/ou_model.h> in every consumption mode.
target_include_directories(stochastic_models
    PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# Boost is third-party and not warning-clean; SYSTEM keeps its headers out
# of the warning gate.
target_include_directories(stochastic_models SYSTEM PUBLIC
    $<BUILD_INTERFACE:${boost_SOURCE_DIR}>
)

target_link_libraries(
    stochastic_models
    PRIVATE
    GSL::gsl
    # BUILD_INTERFACE keeps nlohmann_json out of the installed export.
    $<BUILD_INTERFACE:nlohmann_json::nlohmann_json>
)

# Warnings apply to the library target only; fetched dependencies are not
# warning-clean. Each suppression is temporary and has an owning fix:
#   -Wno-ignored-qualifiers      const-qualified value returns (API cleanup)
#   -Wno-unused-parameter        not-implemented stubs
#   -Wno-delete-non-virtual-dtor pending virtual-destructor fix
target_compile_options(
    stochastic_models
    PRIVATE
    -Wall
    -Wextra
    -Wno-ignored-qualifiers
    -Wno-unused-parameter
    -Wno-delete-non-virtual-dtor
    $<$<BOOL:${STOCHASTIC_MODELS_WERROR}>:-Werror>
)
