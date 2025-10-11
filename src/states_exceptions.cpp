#include "stochastic_models/kalman_filter/states_exceptions.h"

filter_uninitialised::filter_uninitialised(const std::string& message)
    : std::logic_error(message) {}
filter_invalid_operation::filter_invalid_operation(const std::string& message)
    : std::logic_error(message) {}
json_parse_error::json_parse_error(const std::string& message)
    : std::runtime_error(message) {}
