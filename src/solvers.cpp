#include "stochastic_models/numeric_utils/solvers.h"

#include "stochastic_models/exceptions/errors.h"
#include "stochastic_models/exceptions/gsl_errors.h"
#include "stochastic_models/numeric_utils/helpers.h"

#include <gsl/gsl_errno.h>
#include <iostream>
#include <math.h>

BrentSolverState::BrentSolverState() {
  fsolver_type = gsl_root_fsolver_brent;
  fsolver = gsl_root_fsolver_alloc(fsolver_type);
}
BrentSolverState::~BrentSolverState() {
  if (fsolver != nullptr) {
    gsl_root_fsolver_free(fsolver);
    fsolver = nullptr;
  }
}

constexpr double kEpsAbs = 1e-8;
constexpr double kEpsRel = 1e-4;

const double
brentSolver(ModelFunc fn, void* model, double& lower, double& upper) {
  if (lower >= upper) {
    throw std::invalid_argument(
        "Invalid interval: lower bound must be less than upper bound."
    );
  }

  BrentSolverState solver_state;

  if (solver_state.fsolver == nullptr) {
    std::cerr << "Error: failed to allocate memory for solver." << std::endl;
    throw NoMemoryError();
  }

  gsl_function F;
  F.function = fn;
  F.params = model;

  // RAII guard restores the previous GSL handler on every exit path.
  GslHandlerGuard gsl_guard{&custom_gsl_exception_handler};

  int status = gsl_root_fsolver_set(solver_state.fsolver, &F, lower, upper);

  if (status == GSL_EINVAL) {
    throw RootNotBracketedError(
        "Function does not change sign over [" + std::to_string(lower) + ", " +
        std::to_string(upper) + "]."
    );
  }
  check_function_status(status, {});

  int iter = 0, max_iter = 100;
  double result = 0, x_lo = 0, x_hi = 0;
  do {
    iter++;
    status = gsl_root_fsolver_iterate(solver_state.fsolver);
    check_function_status(status, {});
    result = gsl_root_fsolver_root(solver_state.fsolver);
    x_lo = gsl_root_fsolver_x_lower(solver_state.fsolver);
    x_hi = gsl_root_fsolver_x_upper(solver_state.fsolver);
    status = gsl_root_test_interval(x_lo, x_hi, kEpsAbs, kEpsRel);
  } while (status == GSL_CONTINUE && iter < max_iter);

  if (status == GSL_CONTINUE) {
    throw SolverConvergenceError(
        "Root solver did not converge after " + std::to_string(max_iter) +
        " iterations."
    );
  }
  check_function_status(status, {});
  check_finite_result(result, "brentSolver");

  return result;
}
