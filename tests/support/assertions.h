#ifndef STOCHASTIC_MODELS_TESTS_SUPPORT_ASSERTIONS_H
#define STOCHASTIC_MODELS_TESTS_SUPPORT_ASSERTIONS_H

#include <cstddef>
#include <gtest/gtest.h>
#include <nlohmann/json.hpp>

namespace test_support {

  inline double asDouble(const double value) {
    return value;
  }

  inline double asDouble(const nlohmann::json& value) {
    return value.get<double>();
  }

  template <typename Actual, typename Expected>
  void expectVectorNear(
      const Actual& actual,
      const Expected& expected,
      const char* field,
      const double tolerance
  ) {
    ASSERT_EQ(actual.size(), expected.size())
        << "Field '" << field << "' has the wrong length.";
    for (std::size_t i{0}; i < expected.size(); i++) {
      EXPECT_NEAR(asDouble(actual.at(i)), expected.at(i), tolerance)
          << "Field '" << field << "' differs at index " << i << ".";
    }
  }

  template <typename Actual, typename Expected>
  void expectMatrixNear(
      const Actual& actual,
      const Expected& expected,
      const char* field,
      const double tolerance
  ) {
    ASSERT_EQ(actual.size(), expected.size())
        << "Field '" << field << "' has the wrong row count.";
    for (std::size_t i{0}; i < expected.size(); i++) {
      ASSERT_NO_FATAL_FAILURE(
          expectVectorNear(actual.at(i), expected.at(i), field, tolerance)
      );
    }
  }

} // namespace test_support

#endif
