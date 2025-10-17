#ifndef STOCHASTIC_MODELS_LIKELIHOOD_ORNSTEIN_UHLENBECK_ONLINE_H
#define STOCHASTIC_MODELS_LIKELIHOOD_ORNSTEIN_UHLENBECK_ONLINE_H

class OrnsteinUhlenbeckUpdater {
public:
  const double lead_series_sum;
  const double lag_series_sum;
  const double lead_sum_of_squares;
  const double lag_sum_of_squares;
  const double last_observation;
  const int n_obs;

  OrnsteinUhlenbeckUpdater(
      const double lead_series_sum,
      const double lag_series_sum,
      const double lead_sum_of_squares,
      const double lag_sum_of_squares,
      const double last_observation,
      const int n_obs
  );
};

#endif // STOCHASTIC_MODELS_LIKELIHOOD_ORNSTEIN_UHLENBECK_ONLINE_H
