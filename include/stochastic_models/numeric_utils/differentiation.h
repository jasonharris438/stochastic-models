#ifndef STOCHASTIC_MODELS_NUMERIC_UTILS_DIFFERENTIATION_H
#define STOCHASTIC_MODELS_NUMERIC_UTILS_DIFFERENTIATION_H
#include "stochastic_models/numeric_utils/types.h"

/**
 * @file
 * @brief Numeric differentiation helpers.
 */

/**
 * @brief Compute the derivative of fn at x using an adaptive central
 * difference.
 *
 * @param fn Function pointer to evaluate.
 * @param model Opaque context pointer passed through to fn. Contains the model
 * instance that is being used.
 * @param x Point at which to compute derivative (may be adjusted by the
 *          adaptive routine in some implementations).
 * @return const double Approximated derivative value f'(x).
 */
const double
adaptiveCentralDifferentiation(ModelFunc fn, void* model, double& x);
#endif // STOCHASTIC_MODELS_NUMERIC_UTILS_DIFFERENTIATION_H
