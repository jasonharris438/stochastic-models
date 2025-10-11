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
