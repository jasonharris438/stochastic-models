#include "stochastic_models/numeric_utils/differentiation.h"

#include "stochastic_models/exceptions/gsl_errors.h"
#include "stochastic_models/numeric_utils/helpers.h"

#include <gsl/gsl_deriv.h>
const double
adaptiveCentralDifferentiation(ModelFunc fn, void* model, double& x) {
  double result, error;

  gsl_function F;
  F.function = *fn;
  F.params = model;

  // Set custom error handler.
  gsl_error_handler_t* old_handler =
      gsl_set_error_handler(&custom_gsl_exception_handler);

  int status = gsl_deriv_central(&F, x, 1e-5, &result, &error);

  // No codes to ignore.
  const std::vector<int> ignore_codes = {};
  check_function_status(status, ignore_codes);

  // Restore original handler.
  gsl_set_error_handler(old_handler);

  const double value = result;

  return value;
}
