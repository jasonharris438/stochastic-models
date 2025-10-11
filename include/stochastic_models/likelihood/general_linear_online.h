#ifndef GENERAL_LINEAR_ONLINE_H
#define GENERAL_LINEAR_ONLINE_H

class GeneralLinearUpdater {
private:
  double mu_numerator;
  double mu_denominator;
  double sigma_state_value;
  double last_observation;
  double n_obs;

  const double initialSigmaState(const double& sigma, const int& n_obs);
  void applyInternalMuUpdate(
      const double mu_numerator_new,
      const double mu_denominator_new,
      const double new_observation
  );
  void applyInternalSigmaUpdate(
      const double sigma_state_value_new, const int n_obs_new
  );

public:
  GeneralLinearUpdater(
      const double mu_numerator,
      const double mu_denominator,
      const double sigma_state_value,
      const double last_observation,
      const double n_obs
  );
  const double updateMu(const double value);
  const double updateSigma(const double& value, const double& expected_x);
};
#endif // GENERAL_LINEAR_ONLINE_H
