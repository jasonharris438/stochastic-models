#ifndef STOCHASTIC_MODELS_NUMERIC_UTILS_HELPERS_H
#define STOCHASTIC_MODELS_NUMERIC_UTILS_HELPERS_H
#include "stochastic_models/sde/stochastic_model.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <vector>

/**
 * @file
 * @brief Numeric helper utilities used across the project.
 *
 * This header contains convenience helpers (rounding, squared-values,
 * and solver bounds) and an adapter to convert GSL error status codes into
 * exceptions used by the test-suite and library.
 */

/**
 * @brief Checks the status returned by a GSL routine and maps failures to
 * exceptions.
 *
 * The function prints a diagnostic message to stderr and throws an exception
 * when the status indicates an error that is not listed in
 * @p ignore_codes.
 *
 * @param status Status code returned by a GSL routine (GSL_SUCCESS == 0).
 * @param ignore_codes Vector of GSL error codes that should be ignored.
 */
void check_function_status(
    const int& status, const std::vector<int> ignore_codes
);

/**
 * @brief Return a vector containing element-wise squares of the input.
 *
 * This is a small convenience wrapper around std::transform.
 *
 * @param vec Input vector of values.
 * @return const std::vector<double> New vector where each element is the
 *         square of the corresponding input element.
 */
const std::vector<double> valuesSquared(const std::vector<double>& vec);

/**
 * @brief Round a floating-point value to a fixed number of decimal places.
 *
 * The rounding uses std::round after scaling by 10^decimals. The function is
 * stable for typical decimal ranges (e.g. decimals between 0 and 15).
 *
 * @param value Value to be rounded.
 * @param decimals Number of decimal places to round to (>=0).
 * @return const double Rounded value.
 */
const double roundToDecimals(const double& value, const int& decimals);

/**
 * @brief Provide an upper bound to give to numerical solvers.
 *
 * This helper uses the model mean plus four times the unconditional standard
 * deviation which tries to empirically provide a safe search interval for root
 * finding in the library's entrypoint solvers.
 *
 * @param model Pointer to a StochasticModel instance used for obtaining
 *              mean/variance information.
 * @return const double Upper bound for solver search intervals.
 */
const double upperSolverBound(const StochasticModel* model);

/**
 * @brief Provide a lower bound to give to numerical solvers.
 *
 * Mirrors `upperSolverBound` but returns mean - 4 * stddev.
 *
 * @param model Pointer to a StochasticModel instance used for obtaining
 *              mean/variance information.
 * @return const double Lower bound for solver search intervals.
 */
const double lowerSolverBound(const StochasticModel* model);
#endif // STOCHASTIC_MODELS_NUMERIC_UTILS_HELPERS_H
