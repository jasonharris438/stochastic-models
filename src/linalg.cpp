#include "stochastic_models/numeric_utils/linalg.h"

#include <iostream>

#include "stochastic_models/exceptions/gsl_errors.h"
#include "stochastic_models/numeric_utils/helpers.h"
void BoostMartixInverter::linalgLuDecomp(gsl_matrix* gsl_mat,
                                         gsl_permutation* perm) const {
    int signum;
    const int status = gsl_linalg_LU_decomp(gsl_mat, perm, &signum);
    if (status != 0) {
        std::cerr << "error: " << gsl_strerror(status) << std::endl;
        std::cerr << "LU decomposition failed.\n";
        gsl_permutation_free(perm);
        gsl_matrix_free(gsl_mat);
        throw std::runtime_error("GSL LU decomposition failed.");
    }
}
void BoostMartixInverter::linalgInvertMatrix(gsl_matrix* gsl_mat,
                                             gsl_permutation* perm,
                                             gsl_matrix* gsl_inv) const {
    const int status = gsl_linalg_LU_invert(gsl_mat, perm, gsl_inv);
    if (status != 0) {
        std::cerr << "error: " << gsl_strerror(status) << std::endl;
        std::cerr << "Matrix inversion failed.\n";
        gsl_permutation_free(perm);
        gsl_matrix_free(gsl_mat);
        gsl_matrix_free(gsl_inv);
        throw std::runtime_error("GSL matrix inversion failed.");
    }
}
void BoostMartixInverter::copyBoostToGslMatrix(
    const boost::numeric::ublas::matrix<double>& boost_matrix,
    gsl_matrix* gsl_mat) const {
    std::size_t rows = boost_matrix.size1();
    std::size_t cols = boost_matrix.size2();
    for (std::size_t i = 0; i < rows; i++) {
        for (std::size_t j = 0; j < cols; j++) {
            gsl_matrix_set(gsl_mat, i, j, boost_matrix(i, j));
        }
    }
}
const boost::numeric::ublas::matrix<double>
BoostMartixInverter::copyGslToBoostMatrix(const gsl_matrix* gsl_mat) const {
    std::size_t rows = gsl_mat->size1;
    std::size_t cols = gsl_mat->size2;
    boost::numeric::ublas::matrix<double> boost_matrix(rows, cols);
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            boost_matrix(i, j) = gsl_matrix_get(gsl_mat, i, j);
        }
    }
    return boost_matrix;
}

const boost::numeric::ublas::matrix<double> BoostMartixInverter::invertMatrix(
    const boost::numeric::ublas::matrix<double>& boost_matrix) const {
    // Allocate GSL matrix memory space for inverted matrix.
    std::size_t rows = boost_matrix.size1();
    std::size_t cols = boost_matrix.size2();
    gsl_matrix* gsl_mat = gsl_matrix_alloc(rows, cols);

    // Set custom GSL error handler.
    gsl_error_handler_t* old_handler =
        gsl_set_error_handler(&custom_gsl_exception_handler);

    // Copy Boost matrix to GSL matrix.
    copyBoostToGslMatrix(boost_matrix, gsl_mat);

    // Allocate permutation matrix for LU decomposition.
    gsl_permutation* perm = gsl_permutation_alloc(rows);

    // Perform LU decomposition on the GSL matrix.
    linalgLuDecomp(gsl_mat, perm);

    // Allocate GSL matrix for the final inverse matrix.
    gsl_matrix* gsl_inv = gsl_matrix_alloc(rows, cols);

    // Compute the final GSL inverse matrix.
    linalgInvertMatrix(gsl_mat, perm, gsl_inv);

    // Copy elements from the GSL inverse matrix to a boost::numeric::ublas
    // matrix.
    const boost::numeric::ublas::matrix<double> boost_inv_matrix =
        copyGslToBoostMatrix(gsl_inv);

    // Free memory that was allocated for the GSL procedures.
    gsl_permutation_free(perm);
    gsl_matrix_free(gsl_mat);
    gsl_matrix_free(gsl_inv);

    // Restore original GSL error handler.
    gsl_set_error_handler(old_handler);

    return boost_inv_matrix;
}
