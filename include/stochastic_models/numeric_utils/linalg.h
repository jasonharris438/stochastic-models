#ifndef LINALG_H
#define LINALG_H
#include <boost/numeric/ublas/matrix.hpp>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>

/**
 * @brief Class that handles inverting a boost::numeric::ublas matrix using GSL.
 *
 * This logic is unfortunately required due to boost not providing a native
 * matrix inversion method.
 */
class BoostMartixInverter {
private:
  /**
   * @brief Method to perform LU decomposition on a matrix with GSL.
   *
   * This method will throw a std::runtime_error if the LU decomposition
   * fails. The permutation matrix and gsl matrix will be freed before an
   * error is thrown.
   *
   * @param gsl_mat The matrix to be decomposed.
   * @param perm The GSL permutation matrix.
   * @throws std::runtime_error if the LU decomposition fails.
   */
  void linalgLuDecomp(gsl_matrix* gsl_mat, gsl_permutation* perm) const;
  /**
   * @brief Method to invert a matrix with GSL.
   *
   * This method will throw a std::runtime_error if the matrix inversion
   * fails. The permutation matrix and gsl matrices will be freed before an
   * error is thrown.
   *
   * @param gsl_mat The matrix to be inverted.
   * @param perm The GSL permutation matrix memory allocation.
   * @param gsl_inv The inverted matrix memory that has been allocated.
   * @throws std::runtime_error if the matrix inversion fails.
   */
  void linalgInvertMatrix(
      gsl_matrix* gsl_mat, gsl_permutation* perm, gsl_matrix* gsl_inv
  ) const;
  /**
   * @brief Method to copy a boost matrix to a GSL matrix.
   *
   * @param boost_matrix The boost matrix to be copied.
   * @param gsl_mat The GSL matrix memory allocation to be copied to.
   */
  void copyBoostToGslMatrix(
      const boost::numeric::ublas::matrix<double>& boost_matrix,
      gsl_matrix* gsl_mat
  ) const;
  /**
   * @brief Method to copy a GSL matrix to a boost matrix.
   *
   * @param gsl_mat The GSL matrix to be copied.
   * @return The boost matrix that has been copied.
   */
  const boost::numeric::ublas::matrix<double>
  copyGslToBoostMatrix(const gsl_matrix* gsl_mat) const;

public:
  /**
   * @brief Method to calculate the inverse of a matrix.
   * @param matrix The matrix to be inverted.
   * @return The inverse of the matrix.
   */
  const boost::numeric::ublas::matrix<double>
  invertMatrix(const boost::numeric::ublas::matrix<double>& boost_matrix) const;
};

#endif // LINALG_H
