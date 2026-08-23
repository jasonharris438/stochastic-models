#include "stochastic_models/numeric_utils/differentiation.h"

#include "stochastic_models/exceptions/gsl_errors.h"
#include "stochastic_models/numeric_utils/helpers.h"

#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h>
const double
adaptiveCentralDifferentiation(ModelFunc fn, void* model, double& x) {
  double result = 0, error = 0;

  gsl_function F;
  F.function = *fn;
  F.params = model;

  // RAII guard restores the previous GSL handler on every exit path.
  GslHandlerGuard gsl_guard{&custom_gsl_exception_handler};

  int status = gsl_deriv_central(&F, x, 1e-5, &result, &error);

  check_function_status(status, {});
  check_finite_result(result, "adaptiveCentralDifferentiation");

  return result;
}
