#ifndef STOCHASTIC_MODELS_KALMAN_FILTER_TYPE_CONVERSION_H
#define STOCHASTIC_MODELS_KALMAN_FILTER_TYPE_CONVERSION_H
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <nlohmann/json.hpp>
#include <vector>
/**
 * @file
 * @brief Helpers to convert between boost::ublas and std::vector types.
 */

/**
 * @brief Convert a boost matrix into a vector-of-vectors (row-major).
 * @param matrix Input boost matrix.
 * @return std::vector<std::vector<double>> Row-major vector representation.
 */
const std::vector<std::vector<double>> copy_matrix_elements_to_vector(
    const boost::numeric::ublas::matrix<double>& matrix
);
/**
 * @brief Create a boost matrix from row-major vector-of-vectors.
 * @param matrix_as_vectors Row-major vector of vectors.
 * @return boost::numeric::ublas::matrix<double> Constructed boost matrix.
 */
boost::numeric::ublas::matrix<double> create_boost_matrix_from_vectors(
    const std::vector<std::vector<double>>& matrix_as_vectors
);
/**
 * @brief Create a boost vector from an std::vector (copying contents).
 * @param vector_as_vector Input std::vector<double>.
 * @return boost::numeric::ublas::vector<double> Constructed boost vector.
 */
boost::numeric::ublas::vector<double>
create_boost_vector_from_vector(const std::vector<double>& vector_as_vector);
/**
 * @brief Add a scalar in-place to every element of a boost matrix.
 */
void add_scalar_to_matrix(
    boost::numeric::ublas::matrix<double>& boost_matrix, const double& scalar
);
/**
 * @brief Add a scalar in-place to every element of a boost vector.
 */
void add_scalar_to_vector(
    boost::numeric::ublas::vector<double>& boost_vector, const double& scalar
);

#endif // STOCHASTIC_MODELS_KALMAN_FILTER_TYPE_CONVERSION_H
