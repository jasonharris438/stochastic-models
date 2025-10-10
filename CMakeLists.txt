# Minimum version 4.0.
cmake_minimum_required(VERSION 4.0)

# This is a shared library project.
option(BUILD_SHARED_LIBS "Build the shared library" ON)

# Set the project name and version.
project(StochasticModels
    DESCRIPTION "Core stochastic differential equation implementation of fast mathematical tools for stochastic calculus."
    VERSION 0.1.5
    LANGUAGES CXX)

# specify the C++ standard.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Include the FetchContent module.
include(FetchContent)

# Fetch Boost 1.89.0
FetchContent_Declare(
    Boost
    URL      https://archives.boost.io/release/1.89.0/source/boost_1_89_0.tar.bz2
    URL_HASH SHA256=85a33fa22621b4f314f8e85e1a5e2a9363d22e4f4992925d4bb3bc631b5a0c7a
    EXCLUDE_FROM_ALL
    )
FetchContent_MakeAvailable(Boost)

# Fetch nlohmann/json 3.12.0
FetchContent_Declare(
    json
    URL         https://github.com/nlohmann/json/releases/download/v3.12.0/json.tar.xz
    URL_HASH    SHA256=42f6e95cad6ec532fd372391373363b62a14af6d771056dbfc86160e6dfff7aa
)
FetchContent_MakeAvailable(json)
option(JSON_ImplicitConversions OFF) # Recommended by nlohmann/json maintainers.

# Add GSL package.
find_package(GSL REQUIRED)
include(GNUInstallDirs)

# Add the main library.
add_library(${PROJECT_NAME} .)
set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)

# Add subdirectories.
add_subdirectory(src)
add_subdirectory(tests)

# create config file
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)

# CMake sets this variable to a ``TRUE`` value when the
# ``CMAKE_INSTALL_PREFIX`` has just been initialized to
# its default value
if(DEFINED CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
message(
        STATUS
        "CMAKE_INSTALL_PREFIX is not set\n"
        "Default value: ${CMAKE_INSTALL_PREFIX}\n"
        "Will set it to ${CMAKE_SOURCE_DIR}/install"
        )
        set(CMAKE_INSTALL_PREFIX
        "${CMAKE_SOURCE_DIR}/install"
        CACHE PATH "Where the library will be installed to" FORCE
        )
else()
    message(
        STATUS
        "CMAKE_INSTALL_PREFIX was already set\n"
        "Current value: ${CMAKE_INSTALL_PREFIX}"
        )
endif()

install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.h"
)
install(DIRECTORY "${boost_SOURCE_DIR}/boost/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/boost"
    FILES_MATCHING PATTERN "*.hpp"
)
install(DIRECTORY "${json_SOURCE_DIR}/include/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    FILES_MATCHING PATTERN "*.hpp"
)

install(TARGETS ${PROJECT_NAME}
    EXPORT "${PROJECT_NAME}Targets"
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
# generate and install export file
install(EXPORT "${PROJECT_NAME}Targets"
    FILE "${PROJECT_NAME}Targets.cmake"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)
