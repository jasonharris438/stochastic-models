#include "stochastic_models/numeric_utils/integration.h"

#include "stochastic_models/exceptions/errors.h"
#include "stochastic_models/exceptions/gsl_errors.h"
#include "stochastic_models/numeric_utils/helpers.h"

#include <algorithm>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <string>

constexpr double eps_abs = 0;
constexpr double eps_rel = 1e-7;
constexpr double tolerance_slack = 10;
constexpr double bound_floor = 1e-10;

namespace {
  // Gates results whose status was a round-off error: accept only when the
  // error estimate is inside a slack multiple of the requested tolerance.
  // The floor accepts integrands that cancel to nearly zero, where a purely
  // relative bound would collapse.
  void check_error_estimate(
      const double& error, const double& result, const char* routine
  ) {
    const double bound =
        tolerance_slack * std::max(bound_floor, eps_rel * std::abs(result));
    if (!(error <= bound)) {
      throw IntegrationToleranceError(
          std::string(routine) + " error estimate " + std::to_string(error) +
          " is greater than the accepted bound " + std::to_string(bound) + "."
      );
    }
  }
} // namespace

IntegrationState::IntegrationState(gsl_integration_workspace& w)
    : workspace(&w) {}
IntegrationState::~IntegrationState() {
  if (workspace != nullptr) {
    gsl_integration_workspace_free(workspace);
    workspace = nullptr;
  }
}

const double
adaptiveIntegration(ModelFunc fn, void* model, double& lower, double& upper) {
  IntegrationState state(*gsl_integration_workspace_alloc(1000));

  double result = 0, error = 0;

  gsl_function F;
  F.function = *fn;
  F.params = model;

  // RAII guard restores the previous GSL handler on every exit path.
  GslHandlerGuard gsl_guard{&custom_gsl_exception_handler};

  int status = gsl_integration_qags(
      &F, lower, upper, eps_abs, eps_rel, 1000, state.workspace, &result, &error
  );

  // A round-off status is tolerated here and gated by the estimate check.
  check_function_status(status, {GSL_EROUND});
  check_finite_result(result, "adaptiveIntegration");
  check_error_estimate(error, result, "adaptiveIntegration");

  return result;
}
const double
semiInfiniteIntegrationUpper(ModelFunc fn, void* model, double& lower) {
  IntegrationState state(*gsl_integration_workspace_alloc(1000));

  double result = 0, error = 0;

  gsl_function F;
  F.function = *fn;
  F.params = model;

  // RAII guard restores the previous GSL handler on every exit path.
  GslHandlerGuard gsl_guard{&custom_gsl_exception_handler};

  int status = gsl_integration_qagiu(
      &F, lower, eps_abs, eps_rel, 1000, state.workspace, &result, &error
  );

  // A round-off status is tolerated here and gated by the estimate check.
  check_function_status(status, {GSL_EROUND});
  check_finite_result(result, "semiInfiniteIntegrationUpper");
  check_error_estimate(error, result, "semiInfiniteIntegrationUpper");

  return result;
}
