#ifndef STOCHASTIC_MODELS_NUMERIC_UTILS_TYPES_H
#define STOCHASTIC_MODELS_NUMERIC_UTILS_TYPES_H
/**
 * @brief Function pointer type used by the numeric integration helpers.
 *
 * The function takes the integration point x and an opaque pointer to a
 * user-provided model/context structure and returns the function value at x.
 *
 * GSL's convention is heavily dependent on passing context via void pointers.
 * This type alias encapsulates that pattern for use across the numeric
 * utilities.
 */
typedef double (*ModelFunc)(double x, void* model);
#endif // STOCHASTIC_MODELS_NUMERIC_UTILS_TYPES_H
