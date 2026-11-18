#ifndef STOCHASTIC_MODELS_TESTS_SUPPORT_GENERAL_LINEAR_FIXTURES_H
#define STOCHASTIC_MODELS_TESTS_SUPPORT_GENERAL_LINEAR_FIXTURES_H

#include <array>
#include <vector>

namespace test_support {

  /**
   * @brief The 20-point price series shared by the General Linear estimator
   * tests, and the observation that extends it to 21 points.
   */
  inline constexpr std::array<double, 20> general_linear_series{
      1094.1, 1104.1, 1107.7, 1123.6, 1115.6, 1112.7, 1118.4,
      1116.9, 1127.9, 1153.2, 1159.6, 1153.6, 1138.3, 1124.6,
      1122.6, 1134.0, 1132.5, 1139.8, 1133.6, 1124.5
  };
  inline constexpr double general_linear_next_observation{1125.25};

  inline std::vector<double> generalLinearSeries() {
    return std::vector<double>(
        general_linear_series.begin(), general_linear_series.end()
    );
  }

} // namespace test_support

#endif // STOCHASTIC_MODELS_TESTS_SUPPORT_GENERAL_LINEAR_FIXTURES_H
