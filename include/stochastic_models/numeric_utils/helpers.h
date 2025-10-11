#ifndef HELPERS_H
#define HELPERS_H
#include "stochastic_models/sde/stochastic_model.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <vector>
/**
 * @brief Checks the status of the integration routine and throws a runtime
 * error exception if the status is not zero.
 *
 * @param status Reference to the status returned by the integration routine.
 * @param ignore_codes Vector of error codes to ignore.
 */
void check_function_status(
    const int& status, const std::vector<int> ignore_codes
);

/**
 * @brief Squares all values in a vector in-place.
 *
 * @param vec Vector of values to be squared.
 */
void valuesSquared(std::vector<double>& vec);
/**
 * @brief Rounds a value to a given number of decimals.
 *
 * @param value Value to be rounded.
 * @param decimals Number of decimals to round to.
 */
const double roundToDecimals(const double& value, const int& decimals);
/**
 * @brief Returns an upper bound value to provide to a solver for optimal entry
 * and exit levels.
 *
 * @param model The model to use in calculating the upper bound value.
 * @return const double The upper bound value for optimal entry and exit
 * levels.
 */
const double upperSolverBound(const StochasticModel* model);
/**
 * @brief Returns a lower bound value to provide to a solver for optimal entry
 * and exit levels.
 *
 * @param model The model to use in calculating the lower bound value.
 * @return const double The lower bound value for optimal entry and exit
 * levels.
 */
const double lowerSolverBound(const StochasticModel* model);
#endif // HELPERS_H
