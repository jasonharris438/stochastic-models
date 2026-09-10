#ifndef STOCHASTIC_MODELS_KALMAN_FILTER_STATES_EXCEPTIONS_H
#define STOCHASTIC_MODELS_KALMAN_FILTER_STATES_EXCEPTIONS_H

#include <stdexcept>
#include <string>
class filter_uninitialised : public std::logic_error {
public:
  explicit filter_uninitialised(const std::string& message);
};
class filter_invalid_operation : public std::logic_error {
public:
  explicit filter_invalid_operation(const std::string& message);
};
class json_parse_error : public std::runtime_error {
public:
  explicit json_parse_error(const std::string& message);
};
/**
 * @brief A filter dimension is out of range or the set is inconsistent.
 */
class invalid_filter_dimensions final : public std::invalid_argument {
public:
  explicit invalid_filter_dimensions(const std::string& message);
};
/**
 * @brief A source length or shape differs from the target it is copied into.
 */
class filter_shape_mismatch final : public std::invalid_argument {
public:
  explicit filter_shape_mismatch(const std::string& message);
};
#endif // STOCHASTIC_MODELS_KALMAN_FILTER_STATES_EXCEPTIONS_H
