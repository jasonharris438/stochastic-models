#ifndef TYPE_CONVERSION_H
#define TYPE_CONVERSION_H
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <nlohmann/json.hpp>
#include <vector>
/**
 * @brief Function to copy the elements of a boost matrix to a std::vector.
 *
 * @param matrix The boost matrix to copy.
 * @return std::vector<std::vector<double>>. The std::vector containing the data
 * from the matrix.
 */
const std::vector<std::vector<double>> copy_matrix_elements_to_vector(
    const boost::numeric::ublas::matrix<double>& matrix
);
/**
 * @brief Copies elements of a std::vector<std::vector<double>> to a newly
 * initialized boost matrix<double>.
 *
 * @param matrix_as_vectors The vector to copy elements from into a new boost
 * matrix.
 * @return boost::numeric::ublas::matrix<double>. The boost matrix containing
 * the data from the std::vector.
 */
boost::numeric::ublas::matrix<double> create_boost_matrix_from_vectors(
    const std::vector<std::vector<double>>& matrix_as_vectors
);
/**
 * @brief Copies elements of a std::vector<double> to a newly initialized boost
 * vector<double>.
 *
 * @param vector_as_vector The vector to copy elements from into a new boost
 * vector.
 * @return boost::numeric::ublas::vector<double>. The boost vector containing
 * the data from the std::vector.
 */
boost::numeric::ublas::vector<double>
create_boost_vector_from_vector(const std::vector<double>& vector_as_vector);
/**
 * @brief Adds a scalar to each element of a boost matrix.
 *
 * @param boost_matrix The matrix to add the scalar to.
 * @param scalar The scalar to add to each element of the matrix.
 */
void add_scalar_to_matrix(
    boost::numeric::ublas::matrix<double>& boost_matrix, const double& scalar
);
/**
 * @brief Adds a scalar to each element of a boost vector.
 *
 * @param boost_vector The vector to add the scalar to.
 * @param scalar The scalar to add to each element of the vector.
 */
void add_scalar_to_vector(
    boost::numeric::ublas::vector<double>& boost_vector, const double& scalar
);

#endif // TYPE_CONVERSION_H
