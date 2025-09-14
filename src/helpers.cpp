#include "stochastic_models/numeric_utils/helpers.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include "stochastic_models/exceptions/errors.h"

void check_function_status(const int& status,
                           const std::vector<int> ignore_codes) {
    if (status && std::find(ignore_codes.begin(), ignore_codes.end(), status) ==
                      ignore_codes.end()) {
        if (status == GSL_EMAXITER) {
            std::cerr << "error: " << gsl_strerror(status) << std::endl;
            throw IntegrationMaxIterationError(
                "Maximum number of subdivisions reached during numerical "
                "integration.");
        } else if (status == GSL_ENOMEM) {
            std::cerr << "error: " << gsl_strerror(status) << std::endl;
            throw NoMemoryError();
        } else if (status == GSL_EBADFUNC) {
            std::cerr << "error: " << gsl_strerror(status) << std::endl;
            throw NoSolutionError(
                "Root finding solver failed due to no solution at a single "
                "point.");
        } else if (status == GSL_EZERODIV) {
            std::cerr << "error: " << gsl_strerror(status) << std::endl;
            throw ZeroDivError("Derivative vanished in root finding solver.");
        } else if (status == GSL_ESING) {
            std::cerr << "error: " << gsl_strerror(status) << std::endl;
            throw SingularityError(
                "A singularity was encountered during numerical integration.");
        } else if (status != GSL_SUCCESS) {
            std::cerr << "error: " << gsl_strerror(status) << std::endl;
            throw GslFailedError("GSL mathematical function failed.");
        }
    }
}
void valuesSquared(std::vector<double>& vec) {
    for (double& item : vec) {
        item *= item;
    }
}
const double roundToDecimals(const double& value, const int& decimals) {
    const double factor = std::pow(10, decimals);
    return std::round(value * factor) / factor;
}
const double upperSolverBound(const StochasticModel* model) {
    return model->getMean() +
           (4 * std::sqrt(model->getUnconditionalVariance()));
}
const double lowerSolverBound(const StochasticModel* model) {
    return model->getMean() -
           (4 * std::sqrt(model->getUnconditionalVariance()));
}
