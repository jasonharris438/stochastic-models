#include "kca_fuzz_support.h"
#include "stochastic_models/entrypoints/kca_filter.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace {

  const std::vector<double> data_series{10.51255, 10.51985, 10.52405, 10.4656,
                                        10.47,    10.5403,  10.4425,  10.3087,
                                        10.1994,  10.1839,  10.24645, 10.1795,
                                        10.21715, 10.14995, 10.194,   10.22505,
                                        10.27325, 10.25095, 10.30575, 10.27645};
  constexpr double finite_difference_spacing = 1.0;
  constexpr double transition_covariance_scale = 0.001;

} // namespace

extern "C" int
LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size) {
  const std::string dimensions_json(data, data + size);

  std::string output;
  const kca_fuzz::Outcome outcome = kca_fuzz::runClassified([&] {
    output = getInitializedKcaState(
        data_series, finite_difference_spacing, transition_covariance_scale,
        dimensions_json
    );
  });
  if (outcome == kca_fuzz::Outcome::Completed) {
    kca_fuzz::checkRoundTrip(output, dimensions_json);
  }
  return 0;
}
