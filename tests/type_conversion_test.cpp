/**
 * @file type_conversion_test.cpp
 * @brief Unit tests for the kalman_filter type_conversion helpers.
 *
 */
#include "stochastic_models/kalman_filter/type_conversion.h"

#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>

/**
 * @test An empty outer vector previously hit `matrix_as_vectors[0]` (UB). It
 * must be rejected.
 */
TEST(TypeConversionTest, createBoostMatrixRejectsEmptyInput) {
  const std::vector<std::vector<double>> empty{};
  EXPECT_THROW(create_boost_matrix_from_vectors(empty), std::invalid_argument)
      << "create_boost_matrix_from_vectors accepted an empty outer vector.";
}

/**
 * @test A ragged inner vector longer than row 0 previously wrote past the row
 * end. Ragged input must be rejected.
 */
TEST(TypeConversionTest, createBoostMatrixRejectsRaggedInput) {
  const std::vector<std::vector<double>> ragged{
      {1.0, 2.0}, {3.0, 4.0, 5.0}
  };
  EXPECT_THROW(create_boost_matrix_from_vectors(ragged), std::invalid_argument)
      << "create_boost_matrix_from_vectors accepted ragged rows.";
}

/**
 * @test Well-formed input still converts correctly.
 */
TEST(TypeConversionTest, createBoostMatrixConvertsWellFormedInput) {
  const std::vector<std::vector<double>> source{{1.0, 2.0}, {3.0, 4.0}};
  const boost::numeric::ublas::matrix<double> result =
      create_boost_matrix_from_vectors(source);
  ASSERT_EQ(result.size1(), 2u);
  ASSERT_EQ(result.size2(), 2u);
  EXPECT_EQ(result(0, 0), 1.0);
  EXPECT_EQ(result(0, 1), 2.0);
  EXPECT_EQ(result(1, 0), 3.0);
  EXPECT_EQ(result(1, 1), 4.0);
}
