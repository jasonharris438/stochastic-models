cmake_minimum_required(VERSION 3.28)

project(StochasticModelsConsumer LANGUAGES CXX)

find_package(StochasticModels 0.3 REQUIRED CONFIG)

# No C++ standard is set here on purpose: the package must export it.
add_executable(consumer consumer.cpp)
target_link_libraries(consumer PRIVATE StochasticModels::StochasticModels)
