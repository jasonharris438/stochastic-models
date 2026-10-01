#ifndef STOCHASTIC_MODELS_TESTS_SUPPORT_ASSERTIONS_H
#define STOCHASTIC_MODELS_TESTS_SUPPORT_ASSERTIONS_H

#include <cstddef>
#include <gtest/gtest.h>
#include <limits>
#include <nlohmann/json.hpp>

namespace test_support {

  inline double asDouble(const double value) {
    return value;
  }

  inline double asDouble(const nlohmann::json& value) {
    return value.is_number() ? value.get<double>()
                             : std::numeric_limits<double>::quiet_NaN();
  }

  template <typename Actual, typename Expected>
  void expectVectorNear(
      const Actual& actual,
      const Expected& expected,
      const double tolerance,
      const char* message
  ) {
    ASSERT_EQ(actual.size(), expected.size()) << message << " Wrong length.";
    for (std::size_t i{0}; i < expected.size(); i++) {
      EXPECT_NEAR(asDouble(actual.at(i)), expected.at(i), tolerance)
          << message << " Element [" << i << "] differs.";
    }
  }

  template <typename Actual, typename Expected>
  void expectMatrixNear(
      const Actual& actual,
      const Expected& expected,
      const double tolerance,
      const char* message
  ) {
    ASSERT_EQ(actual.size(), expected.size()) << message << " Wrong row count.";
    for (std::size_t row{0}; row < expected.size(); row++) {
      ASSERT_EQ(actual.at(row).size(), expected.at(row).size())
          << message << " Row " << row << " has the wrong length.";
      for (std::size_t column{0}; column < expected.at(row).size(); column++) {
        EXPECT_NEAR(
            asDouble(actual.at(row).at(column)), expected.at(row).at(column),
            tolerance
        ) << message
          << " Element [" << row << "][" << column << "] differs.";
      }
    }
  }

} // namespace test_support

#endif // STOCHASTIC_MODELS_TESTS_SUPPORT_ASSERTIONS_H
