#include <stdexcept>
#include <string>

#ifndef STOCHASTIC_MODELS_EXCEPTIONS_ERRORS_H
#define STOCHASTIC_MODELS_EXCEPTIONS_ERRORS_H

#include <stdexcept>
#include <string>

class IntegrationMaxIterationError final : public std::runtime_error {
public:
  IntegrationMaxIterationError(const std::string& msg)
      : std::runtime_error{msg} {}
};

class GslFailedError final : public std::runtime_error {
public:
  GslFailedError(const std::string& msg) : std::runtime_error{msg} {}
};

class NoSolutionError final : public std::runtime_error {
public:
  NoSolutionError(const std::string& msg) : std::runtime_error{msg} {}
};

class ZeroDivError final : public std::runtime_error {
public:
  ZeroDivError(const std::string& msg) : std::runtime_error{msg} {}
};

class NoMemoryError final : public std::bad_alloc {
public:
  NoMemoryError() : std::bad_alloc{} {}
};

class SingularityError final : public std::runtime_error {
public:
  SingularityError(const std::string& msg) : std::runtime_error{msg} {}
};

class InvalidNumberObservationsError final : public std::logic_error {
public:
  InvalidNumberObservationsError(const std::string& msg)
      : std::logic_error{msg} {}
};
#endif // STOCHASTIC_MODELS_EXCEPTIONS_ERRORS_H
