# libFuzzer supplies main.
foreach(fuzzer kca_state_fuzzer kca_dimensions_fuzzer)
    add_executable(${fuzzer} ${fuzzer}.cpp)
    target_compile_options(
        ${fuzzer}
        PRIVATE
        -fsanitize=fuzzer
        ${STOCHASTIC_MODELS_WARNING_FLAGS}
    )
    target_link_options(${fuzzer} PRIVATE -fsanitize=fuzzer)
    target_link_libraries(
        ${fuzzer}
        PRIVATE
        stochastic_models
        # The support header includes library-internal headers that use json.
        nlohmann_json::nlohmann_json
    )
endforeach()

# libFuzzer writes new inputs to the first corpus directory it is given.
file(MAKE_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}/corpus/state
    ${CMAKE_CURRENT_BINARY_DIR}/corpus/dimensions
)
