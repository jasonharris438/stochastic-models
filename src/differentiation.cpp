#include "stochastic_models/numeric_utils/differentiation.h"

#include "stochastic_models/exceptions/gsl_errors.h"
#include "stochastic_models/numeric_utils/helpers.h"

#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h>
const double
adaptiveCentralDifferentiation(ModelFunc fn, void* model, double& x) {
  double result, error;

  gsl_function F;
  F.function = *fn;
  F.params = model;

  // Set custom error handler; RAII guard restores the previous handler on
  // every exit path (including exceptions).
  GslHandlerGuard gsl_guard{&custom_gsl_exception_handler};

  int status = gsl_deriv_central(&F, x, 1e-5, &result, &error);

  // No codes to ignore.
  const std::vector<int> ignore_codes = {};
  check_function_status(status, ignore_codes);

  const double value = result;

  return value;
}
