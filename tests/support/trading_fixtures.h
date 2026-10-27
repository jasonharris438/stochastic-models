#ifndef STOCHASTIC_MODELS_TESTS_SUPPORT_TRADING_FIXTURES_H
#define STOCHASTIC_MODELS_TESTS_SUPPORT_TRADING_FIXTURES_H

#include <gtest/gtest.h>

namespace test_support {

  /**
   * @brief One trading parameter set: the Ornstein-Uhlenbeck model parameters
   * with the trading cost and the discount rate.
   */
  struct TradingParameters {
    double alpha;
    double mu;
    double sigma;
    double cost;
    double rate;
  };

  inline constexpr TradingParameters ornstein_uhlenbeck_set{
      8.0, 0.3, 0.3, 0.02, 0.05
  };
  inline constexpr TradingParameters exponential_set{
      5.0, 1.3499, 0.15, 0.02, 0.05
  };
  inline constexpr TradingParameters stop_loss_value_set{
      16.6677, 0.5388, 0.1599, 0.05, 0.05
  };
  inline constexpr TradingParameters no_solution_set{
      0.000116, 1.818978, 0.006623, 0.02, 0.05
  };

  /**
   * @brief Fixture that exposes 1 parameter set as named members.
   */
  template <const TradingParameters& ParameterSet>
  class TradingFixture : public testing::Test {
  protected:
    const double alpha{ParameterSet.alpha};
    const double mu{ParameterSet.mu};
    const double sigma{ParameterSet.sigma};
    const double cost{ParameterSet.cost};
    const double rate{ParameterSet.rate};
  };

} // namespace test_support

#endif // STOCHASTIC_MODELS_TESTS_SUPPORT_TRADING_FIXTURES_H
