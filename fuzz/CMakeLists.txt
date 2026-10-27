# The root CMakeLists.txt requires Clang and turns on the sanitizers before it
# adds this directory. libFuzzer supplies main.
foreach(fuzzer kca_state_fuzzer kca_dimensions_fuzzer)
    add_executable(${fuzzer} ${fuzzer}.cpp)
    target_compile_options(
        ${fuzzer}
        PRIVATE
        -fsanitize=fuzzer
        -Wall
        -Wextra
        -Wno-ignored-qualifiers
        $<$<BOOL:${STOCHASTIC_MODELS_WERROR}>:-Werror>
    )
    target_link_options(${fuzzer} PRIVATE -fsanitize=fuzzer)
    target_link_libraries(
        ${fuzzer}
        PRIVATE
        stochastic_models
        nlohmann_json::nlohmann_json
    )
endforeach()

# libFuzzer writes the inputs it finds to the first corpus directory it is
# given; these keep them out of the source tree.
file(MAKE_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}/corpus/state
    ${CMAKE_CURRENT_BINARY_DIR}/corpus/dimensions
)
