#include <stdexcept>
#include <string>

#ifndef STOCHASTIC_MODELS_EXCEPTIONS_ERRORS_H
#define STOCHASTIC_MODELS_EXCEPTIONS_ERRORS_H

/**
 * @file
 * @brief Domain-specific exception types used across numeric helpers.
 */

/**
 * @brief Thrown when integration reaches the maximum allowed subdivisions.
 */
class IntegrationMaxIterationError final : public std::runtime_error {
public:
  IntegrationMaxIterationError(const std::string& msg)
      : std::runtime_error{msg} {}
};

/**
 * @brief Generic GSL failure mapped into a C++ exception.
 */
class GslFailedError final : public std::runtime_error {
public:
  GslFailedError(const std::string& msg) : std::runtime_error{msg} {}
};

/**
 * @brief Thrown when a root-finding routine fails to find a solution.
 */
class NoSolutionError final : public std::runtime_error {
public:
  NoSolutionError(const std::string& msg) : std::runtime_error{msg} {}
};

/**
 * @brief Division-by-zero encountered in a numeric routine.
 */
class ZeroDivError final : public std::runtime_error {
public:
  ZeroDivError(const std::string& msg) : std::runtime_error{msg} {}
};

/**
 * @brief Allocation failure mapped to std::bad_alloc.
 */
class NoMemoryError final : public std::bad_alloc {
public:
  NoMemoryError() : std::bad_alloc{} {}
};

/**
 * @brief Singular matrix / integrand encountered during integration.
 */
class SingularityError final : public std::runtime_error {
public:
  SingularityError(const std::string& msg) : std::runtime_error{msg} {}
};

/**
 * @brief Invalid number of observations supplied to a likelihood helper.
 */
class InvalidNumberObservationsError final : public std::logic_error {
public:
  InvalidNumberObservationsError(const std::string& msg)
      : std::logic_error{msg} {}
};

#endif // STOCHASTIC_MODELS_EXCEPTIONS_ERRORS_H
