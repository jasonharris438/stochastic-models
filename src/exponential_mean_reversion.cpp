/**
 * @file exponential_mean_reversion.cpp
 * @brief Source file for the ExponentialMeanReversion class.
 *
 */
#include "stochastic_models/trading/exponential_mean_reversion.h"

#include "stochastic_models/numeric_utils/differentiation.h"
#include "stochastic_models/numeric_utils/integration.h"
#include "stochastic_models/sde/ornstein_uhlenbeck.h"
#include "stochastic_models/trading/optimal_mean_reversion.h"
#include "stochastic_models/trading/trading_levels_params.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
const ExponentialMeanReversion* ExponentialMeanReversion::clone() const {
  return new ExponentialMeanReversion(*this);
}
const double ExponentialMeanReversion::F(
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& x,
    const double& r,
    const double& c
) const {
  // First create a deep copy of the model pointer and copy the contents of
  // model into the temporary location. This is because we are going to
  // free that memory in the destructor of the OptimalMeanReversionParams
  // struct.
  const HittingTimeOrnsteinUhlenbeck* temp_ptr = hitting_time_kernel->clone();

  // Create new instance of the OptimalMeanReversionParams struct as a
  // void pointer as this is required by the GSL numerical integration
  // function.
  void* params = new OptimalMeanReversionParams{temp_ptr, x, r};
  ModelFunc fn = funcOptimalMeanReversionF;
  double lower = 0;
  double value{0.0};
  try {
    value = semiInfiniteIntegrationUpper(fn, params, lower);
  } catch (const std::exception& e) {
    std::cout << "Exception " << e.what()
              << " caught in ExponentialMeanReversion::F." << std::endl;
    delete static_cast<OptimalMeanReversionParams*>(params);
    params = nullptr;
    temp_ptr = nullptr;
    throw;
  }

  // Then cast the void pointer back to the original type and free the memory.
  OptimalMeanReversionParams* ptr =
      static_cast<OptimalMeanReversionParams*>(params);
  delete ptr;
  ptr = nullptr;
  params = nullptr;
  temp_ptr = nullptr;

  return value;
}
const double ExponentialMeanReversion::G(
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& x,
    const double& r,
    const double& c
) const {
  // First create a deep copy of the model pointer and copy the contents of
  // model into the temporary location. This is because we are going to free
  // that memory in the destructor of the OptimalMeanReversionParams
  // struct.
  const HittingTimeOrnsteinUhlenbeck* temp_ptr = hitting_time_kernel->clone();

  // Create new instance of the OptimalMeanReversionParams struct as a
  // void pointer as this is required by the GSL numerical integration
  // function.
  void* params = new OptimalMeanReversionParams{temp_ptr, x, r};
  ModelFunc fn = funcOptimalMeanReversionG;
  double lower = 0;
  double value{0.0};
  try {
    value = semiInfiniteIntegrationUpper(fn, params, lower);
  } catch (const std::exception& e) {
    std::cout << "Exception " << e.what()
              << " caught in ExponentialMeanReversion::G." << std::endl;
    delete static_cast<OptimalMeanReversionParams*>(params);
    params = nullptr;
    temp_ptr = nullptr;
    throw;
  }

  // Then cast the void pointer back to the original type and free the memory.
  OptimalMeanReversionParams* ptr =
      static_cast<OptimalMeanReversionParams*>(params);
  delete ptr;
  ptr = nullptr;
  params = nullptr;
  temp_ptr = nullptr;

  return value;
}
const double ExponentialMeanReversion::L_star(
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& r,
    const double& c
) const {
  return hitting_time_kernel->optimalTradingLCore(r, c);
}
const double ExponentialMeanReversion::b(
    const double& value,
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& r,
    const double& c
) const {
  ModelFunc fn = funcIntegrateF;
  const double result =
      std::exp(value) *
          ExponentialMeanReversion::F(hitting_time_kernel, value, r, c) -
      (std::exp(value) - c) *
          instantaneousDifferential(fn, hitting_time_kernel, value, r, c);
  return result;
}
const double ExponentialMeanReversion::b(
    const double& value,
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& stop_loss,
    const double& r,
    const double& c
) const {
  // There is no mathematical definition for the optimal exit level of
  // an exponential model with a stop loss parameter.
  throw std::logic_error(
      "ExponentialMeanReversion does not implement function "
      "b with stop_loss parameter"
  );
}
const double ExponentialMeanReversion::d(
    const double& value,
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& b_star,
    const double& r,
    const double& c
) const {
  ModelFunc funcG = funcIntegrateG;
  ModelFunc funcV = valueFunc;
  const double result =
      (ExponentialMeanReversion::G(hitting_time_kernel, value, r, c) *
       (instantaneousDifferential(
            funcV, hitting_time_kernel, value, b_star, r, c
        ) -
        std::exp(value))) -
      (instantaneousDifferential(
           funcG, hitting_time_kernel, value, b_star, r, c
       ) *
       (ExponentialMeanReversion::V(hitting_time_kernel, value, b_star, r, c) -
        std::exp(value) - c));

  return result;
}
const double ExponentialMeanReversion::d(
    const double& value,
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& b_star,
    const double& stop_loss,
    const double& r,
    const double& c
) const {
  // There is no mathematical definition for the optimal entry level of
  // an exponential model with a stop loss parameter.
  throw std::logic_error(
      "ExponentialMeanReversion does not implement function "
      "d with stop_loss parameter"
  );
}
const double ExponentialMeanReversion::V(
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& x,
    const double& b_star,
    const double& r,
    const double& c
) const {
  if (x < b_star) {
    return (std::exp(b_star) - c) * F(hitting_time_kernel, x, r, c) /
           F(hitting_time_kernel, b_star, r, c);
  } else {
    return std::exp(x) - c;
  }
}
const double ExponentialMeanReversion::V(
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& x,
    const double& b_star,
    const double& stop_loss,
    const double& r,
    const double& c
) const {
  // There is no mathematical definition for the value function V of
  // an exponential model with a stop loss parameter.
  throw std::logic_error(
      "ExponentialMeanReversion does not implement function "
      "V with stop_loss parameter"
  );
}
const double ExponentialMeanReversion::a(
    const double& value,
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& b_star,
    const double& stop_loss,
    const double& r,
    const double& c
) const {
  // There is no mathematical definition for the lower optimal entry level of
  // an exponential model with a stop loss parameter.
  throw std::logic_error(
      "ExponentialMeanReversion does not implement function "
      "a with stop_loss parameter"
  );
}
const double ExponentialMeanReversion::a(
    const double& value,
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& b_star,
    const double& r,
    const double& c
) const {
  ModelFunc funcF = funcIntegrateF;
  ModelFunc funcV = valueFunc;
  const double result =
      (ExponentialMeanReversion::F(hitting_time_kernel, value, r, c) *
       (instantaneousDifferential(
            funcV, hitting_time_kernel, value, b_star, r, c
        ) -
        std::exp(value))) -
      (instantaneousDifferential(
           funcF, hitting_time_kernel, value, b_star, r, c
       ) *
       (ExponentialMeanReversion::V(hitting_time_kernel, value, b_star, r, c) -
        std::exp(value) - c));
  return result;
}
const double ExponentialMeanReversion::instantaneousDifferential(
    ModelFunc fn,
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& x,
    const double& b_star,
    const double& stop_loss,
    const double& r,
    const double& c
) const {
  // There is no mathematical definition for the instantaneous differential of
  // an exponential model with a stop loss parameter.
  throw std::logic_error(
      "ExponentialMeanReversion does not implement function "
      "instantaneousDifferential with stop_loss parameter"
  );
}
const double ExponentialMeanReversion::instantaneousDifferential(
    ModelFunc fn,
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& x,
    const double& b_star,
    const double& r,
    const double& c
) const {
  // Create a deep copy of the model and optimizer pointers to ensure that
  // cleanup does not delete heap memory that is not owned by this scope.
  const HittingTimeOrnsteinUhlenbeck* temp_kernel =
      hitting_time_kernel->clone();
  const ExponentialMeanReversion* temp_optimizer = clone();

  // Create copy to not modify the original value and
  // not provide to function that takes a non-const.
  double x_copied = x;

  // Adaptive differentiation function.
  void* params =
      new EntryLevelParams{temp_optimizer, temp_kernel, b_star, r, c};
  double value = adaptiveCentralDifferentiation(fn, params, x_copied);
  struct EntryLevelParams* p = static_cast<EntryLevelParams*>(params);
  delete p;
  p = nullptr;
  params = nullptr;
  return value;
}
const double ExponentialMeanReversion::instantaneousDifferential(
    ModelFunc fn,
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& x,
    const double& r,
    const double& c
) const {
  // Create a deep copy of the model and optimizer pointers to ensure that
  // cleanup does not delete heap memory that is not owned by this scope.
  const HittingTimeOrnsteinUhlenbeck* temp_model = hitting_time_kernel->clone();
  const ExponentialMeanReversion* temp_optimizer = clone();

  // Create copy to not modify the original value and
  // not provide to function that takes a non-const.
  double x_copied = x;

  // Adaptive differentiation function.
  void* params = new ExitLevelParams{temp_optimizer, temp_model, r, c};
  double value = adaptiveCentralDifferentiation(fn, params, x_copied);
  struct ExitLevelParams* p = static_cast<ExitLevelParams*>(params);
  delete p;
  p = nullptr;
  params = nullptr;
  return value;
}
