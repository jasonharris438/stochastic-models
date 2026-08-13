#include <iostream>
#include <stochastic_models/entrypoints/ou_model.h>

// Consumes the installed package: simulates an OU path, then fits it back.
int main() {
  const std::vector<double> simulated_path =
      simulateOrnsteinUhlenbeck(0.5, 0.01, 0.0067, 0, 5000, 1);
  std::cout << "Simulated points: " << simulated_path.size() << "\n";

  const std::unordered_map<std::string, const double> estimates =
      ornsteinUhlenbeckMaximumLikelihood(simulated_path);
  for (const auto& estimate : estimates) {
    std::cout << estimate.first << " " << estimate.second << "\n";
  }
  return 0;
}
