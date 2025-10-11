#include "stochastic_models/exceptions/gsl_errors.h"

#include <iostream>

void custom_gsl_exception_handler(
    const char* reason, const char* file, int line, int gsl_errno
) {
  // Log the error details
  std::cerr << "GSL Error: " << reason << " in " << file << ":" << line
            << " (Error Code: " << gsl_errno << ")" << std::endl;
}
