#ifndef STOCHASTIC_MODELS_EXCEPTIONS_GSL_ERRORS_H
#define STOCHASTIC_MODELS_EXCEPTIONS_GSL_ERRORS_H
#include <gsl/gsl_errno.h>

/**
 * @brief Custom error handler for GSL errors. Does not call abort() after
 * reporting the error.
 *
 * @param reason The reason for the error.
 * @param file The file in which the error occurred.
 * @param line The line number at which the error occurred.
 * @param gsl_errno The GSL error number.
 */
void custom_gsl_exception_handler(
    const char* reason, const char* file, int line, int gsl_errno
);

/**
 * @brief RAII guard that installs a GSL error handler and restores the previous
 * one on destruction (i.e. on every exit path, including exceptions).
 */
struct GslHandlerGuard {
  gsl_error_handler_t* old;
  explicit GslHandlerGuard(gsl_error_handler_t* handler)
      : old(gsl_set_error_handler(handler)) {}
  ~GslHandlerGuard() {
    gsl_set_error_handler(old);
  }
  GslHandlerGuard(const GslHandlerGuard&) = delete;
  GslHandlerGuard& operator=(const GslHandlerGuard&) = delete;
};

#endif // STOCHASTIC_MODELS_EXCEPTIONS_GSL_ERRORS_H
