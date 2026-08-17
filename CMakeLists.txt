cmake_minimum_required(VERSION 3.28)

# This is a shared library project.
option(BUILD_SHARED_LIBS "Build the shared library" ON)

option(STOCHASTIC_MODELS_WERROR "Treat library warnings as errors" OFF)

# Set the project name and version.
project(StochasticModels
    DESCRIPTION "Core stochastic differential equation implementation of fast mathematical tools for stochastic calculus."
    VERSION 0.3.0
    LANGUAGES CXX)

# An empty build type compiles with no optimization and no debug info.
# Multi-config generators pick the configuration at build time instead.
if(PROJECT_IS_TOP_LEVEL AND NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Release Debug RelWithDebInfo MinSizeRel)
endif()

# specify the C++ standard.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Include the FetchContent module.
include(FetchContent)

# Fetch Boost 1.90.0
FetchContent_Declare(
    Boost
    URL      https://archives.boost.io/release/1.90.0/source/boost_1_90_0.tar.bz2
    URL_HASH SHA256=49551aff3b22cbc5c5a9ed3dbc92f0e23ea50a0f7325b0d198b705e8ee3fc305
    )
FetchContent_MakeAvailable(Boost)

# Fetch nlohmann/json 3.12.0
FetchContent_Declare(
    json
    URL         https://github.com/nlohmann/json/releases/download/v3.12.0/json.tar.xz
    URL_HASH    SHA256=42f6e95cad6ec532fd372391373363b62a14af6d771056dbfc86160e6dfff7aa
)
# Emit json headers as -isystem so they stay out of the warning gate.
set(JSON_SystemInclude ON)
# Off before the fetch, or the dependency configures with conversions on.
# Recommended by the nlohmann/json maintainers for untrusted input.
set(JSON_ImplicitConversions OFF CACHE INTERNAL "")
FetchContent_MakeAvailable(json)

# Add GSL package.
find_package(GSL REQUIRED)
include(GNUInstallDirs)
# Explicit: configure_package_config_file and write_basic_package_version_file
# must not depend on a fetched dependency including this module first.
include(CMakePackageConfigHelpers)

# CTest defines BUILD_TESTING (default ON) and calls enable_testing() at the root.
if(PROJECT_IS_TOP_LEVEL)
    include(CTest)
endif()

add_subdirectory(src)
if(PROJECT_IS_TOP_LEVEL AND BUILD_TESTING)
    add_subdirectory(tests)
endif()

# Namespaced alias: add_subdirectory consumers link the same name the package exports.
add_library(StochasticModels::StochasticModels ALIAS stochastic_models)

if(PROJECT_IS_TOP_LEVEL)
    add_subdirectory(examples/quickstart)
endif()

# create config file
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

# A consumer's install prefix belongs to the consumer.
if(PROJECT_IS_TOP_LEVEL)
    if(DEFINED CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
        message(STATUS
            "CMAKE_INSTALL_PREFIX is not set\n"
            "Default value: ${CMAKE_INSTALL_PREFIX}\n"
            "Will set it to ${CMAKE_SOURCE_DIR}/install")
        set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/install"
            CACHE PATH "Where the library will be installed to" FORCE)
    else()
        message(STATUS
            "CMAKE_INSTALL_PREFIX was already set\n"
            "Current value: ${CMAKE_INSTALL_PREFIX}")
    endif()
endif()

install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.h"
)

install(TARGETS stochastic_models
    EXPORT "${PROJECT_NAME}Targets"
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(EXPORT "${PROJECT_NAME}Targets"
    FILE "${PROJECT_NAME}Targets.cmake"
    NAMESPACE StochasticModels::
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)
