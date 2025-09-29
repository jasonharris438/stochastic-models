#include "stochastic_models/numeric_utils/solvers.h"

#include <math.h>

#include <iostream>

#include "stochastic_models/exceptions/errors.h"
#include "stochastic_models/exceptions/gsl_errors.h"
#include "stochastic_models/numeric_utils/helpers.h"

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

const double brentSolver(ModelFunc fn, void* model, double& lower,
                         double& upper) {
    if (lower >= upper) {
        throw std::invalid_argument(
            "Invalid interval: lower bound must be less than upper bound.");
    }

    BrentSolverState solver_state;

    // Catch error if solver cannot be allocated.
    if (solver_state.fsolver == nullptr) {
        std::cerr << "Error: failed to allocate memory for solver."
                  << std::endl;
        throw NoMemoryError();
    }

    gsl_function F;
    F.function = fn;
    F.params = model;

    // Set custom error handler.
    gsl_error_handler_t* old_handler =
        gsl_set_error_handler(&custom_gsl_exception_handler);

    int status = gsl_root_fsolver_set(solver_state.fsolver, &F, lower, upper);

    // We are choosing to ignore an invalid interval as we aren't always
    // straddling y = 0.
    const std::vector<int> ignore_codes = {GSL_EINVAL};
    check_function_status(status, ignore_codes);

    int iter = 0, max_iter = 100;
    double result = 0, x_lo = 0, x_hi = 0;
    do {
        iter++;
        status = gsl_root_fsolver_iterate(solver_state.fsolver);
        result = gsl_root_fsolver_root(solver_state.fsolver);
        x_lo = gsl_root_fsolver_x_lower(solver_state.fsolver);
        x_hi = gsl_root_fsolver_x_upper(solver_state.fsolver);
        status = gsl_root_test_interval(x_lo, x_hi, 0, 0.0001);

    } while (status == GSL_CONTINUE && iter < max_iter);

    // Restore original handler.
    gsl_set_error_handler(old_handler);

    const double value = result;

    return value;
}
