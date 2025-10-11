#ifndef STATES_EXCEPTIONS_H
#define STATES_EXCEPTIONS_H

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
#endif // STATES_EXCEPTIONS_H
