#ifndef GSL_ERRORS_H
#define GSL_ERRORS_H

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

#endif // GSL_ERRORS_H
