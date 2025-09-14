#include "stochastic_models/kalman_filter/type_conversion.h"

#include <vector>
const std::vector<std::vector<double>> copy_matrix_elements_to_vector(
    const boost::numeric::ublas::matrix<double>& matrix) {
    // Create a std::vector to hold the matrix values.
    std::vector<std::vector<double>> result(
        matrix.size1(), std::vector<double>(matrix.size2()));

    // Copy the values from the matrix to the std::vector by row.
    for (int i{0}; i < matrix.size1(); i++) {
        std::copy(row(matrix, i).begin(), row(matrix, i).end(),
                  result.at(i).begin());
    }
    return result;
}
boost::numeric::ublas::matrix<double> create_boost_matrix_from_vectors(
    const std::vector<std::vector<double>>& matrix_as_vectors) {
    // Allocate a boost matrix with the same dimensions as the input vectors.
    boost::numeric::ublas::matrix<double> boost_matrix(
        matrix_as_vectors.size(), matrix_as_vectors[0].size());

    // Copy the values from the std::vector<std::vector<>>.
    for (int i{0}; i < matrix_as_vectors.size(); i++) {
        std::copy(matrix_as_vectors.at(i).begin(),
                  matrix_as_vectors.at(i).end(), row(boost_matrix, i).begin());
    }
    return boost_matrix;
}
boost::numeric::ublas::vector<double> create_boost_vector_from_vector(
    const std::vector<double>& vector_as_vector) {
    // Allocate a boost vector with the same dimensions as the input vector.
    boost::numeric::ublas::vector<double> boost_vector(vector_as_vector.size());

    // Copy the values from the std::vector<>.
    std::copy(vector_as_vector.begin(), vector_as_vector.end(),
              boost_vector.begin());
    return boost_vector;
}
void add_scalar_to_matrix(boost::numeric::ublas::matrix<double>& boost_matrix,
                          const double& scalar) {
    for (size_t i{0}; i < boost_matrix.size1(); i++)
        for (size_t j{0}; j < boost_matrix.size2(); j++)
            boost_matrix(i, j) += scalar;
}
void add_scalar_to_vector(boost::numeric::ublas::vector<double>& boost_vector,
                          const double& scalar) {
    for (size_t i{0}; i < boost_vector.size(); i++) boost_vector(i) += scalar;
}
