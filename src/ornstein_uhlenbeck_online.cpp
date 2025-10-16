#include "stochastic_models/likelihood/ornstein_uhlenbeck_online.h"

#include <cmath>

OrnsteinUhlenbeckUpdater::OrnsteinUhlenbeckUpdater(
    const double lead_series_sum,
    const double lag_series_sum,
    const double lead_sum_of_squares,
    const double lag_sum_of_squares,
    const double last_observation,
    const int n_obs
)
    : lead_series_sum(lead_series_sum), lag_series_sum(lag_series_sum),
      lead_sum_of_squares(lead_sum_of_squares),
      lag_sum_of_squares(lag_sum_of_squares),
      last_observation(last_observation), n_obs(n_obs) {}
