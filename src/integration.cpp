#include "stochastic_models/numeric_utils/integration.h"

#include "stochastic_models/exceptions/gsl_errors.h"
#include "stochastic_models/numeric_utils/helpers.h"

#include <iostream>

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

  double result, error;

  gsl_function F;
  F.function = *fn;
  F.params = model;

  /* set custom error handler */
  gsl_error_handler_t* old_handler =
      gsl_set_error_handler(&custom_gsl_exception_handler);

  int status = gsl_integration_qags(
      &F, lower, upper, 0, 1e-7, 1000, state.workspace, &result, &error
  );

  /* restore original handler */
  gsl_set_error_handler(old_handler);

  // Check the status returned from the integration routine.
  // We are choosing to ignore a round-off error as it is not critical to the
  // current use-case.
  const std::vector<int> ignore_codes = {GSL_EROUND};
  check_function_status(status, ignore_codes);

  const double value = result;

  return value;
}
const double
semiInfiniteIntegrationUpper(ModelFunc fn, void* model, double& lower) {
  IntegrationState state(*gsl_integration_workspace_alloc(1000));

  double result, error = 0;

  gsl_function F;
  F.function = *fn;
  F.params = model;

  /* set custom error handler */
  gsl_error_handler_t* old_handler =
      gsl_set_error_handler(&custom_gsl_exception_handler);

  int status = gsl_integration_qagiu(
      &F, lower, 0, 1e-7, 1000, state.workspace, &result, &error
  );

  /* restore original handler */
  gsl_set_error_handler(old_handler);

  // Check the status returned from the integration routine.
  // We are choosing to ignore a round-off error as it is not critical to the
  // current use-case.
  const std::vector<int> ignore_codes = {GSL_EROUND};
  check_function_status(status, ignore_codes);

  const double value = result;

  return value;
}
