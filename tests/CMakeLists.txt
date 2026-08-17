# Add googletest installation.
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  # v1.17.0
  GIT_TAG        "52eb8108c5bdec04579160ae17225d66034bd723"
)
# Off before the fetch, or googletest installs its headers and libraries
# into the package prefix alongside the real library.
set(INSTALL_GTEST OFF)
FetchContent_MakeAvailable(googletest)

enable_testing()

# Create the tests executable.
add_executable(
    unit_tests
    adapters_test.cpp
    exponential_mean_reversion_test.cpp
    filter_states_test.cpp
    filter_update_test.cpp
    gaussian_distribution_test.cpp
    gaussian_seeding_test.cpp
    general_linear_likelihood_test.cpp
    general_linear_online_test.cpp
    general_linear_test.cpp
    general_sde_model_test.cpp
    gsl_handler_guard_test.cpp
    helpers_include_hygiene_test.cpp
    hitting_time_test.cpp
    kca_filter_test.cpp
    optimal_mean_reversion_test.cpp
    optimal_trading_levels_test.cpp
    ornstein_uhlenbeck_likelihood_test.cpp
    ornstein_uhlenbeck_test.cpp
    ou_model_test.cpp
    trading_levels_test.cpp
    type_conversion_test.cpp
    utils_test.cpp)


target_include_directories(unit_tests
    PRIVATE
    # where the library itself will look for its internal headers
    "${PROJECT_SOURCE_DIR}/include"
    )
target_link_libraries(
    unit_tests
    PRIVATE
    GTest::gtest_main
    stochastic_models
    # tests include library-internal headers and call GSL directly
    nlohmann_json::nlohmann_json
    GSL::gsl
)

include(GoogleTest)
gtest_discover_tests(unit_tests)
