#include "kca_fuzz_support.h"
#include "stochastic_models/entrypoints/kca_filter.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace {

  // 1 length byte, then the dimensions JSON.
  constexpr std::size_t prefix_size = 1;
  constexpr std::array<double, 20> observations{
      10.51255, 10.51985, 10.52405, 10.4656,  10.47,    10.5403,  10.4425,
      10.3087,  10.1994,  10.1839,  10.24645, 10.1795,  10.21715, 10.14995,
      10.194,   10.22505, 10.27325, 10.25095, 10.30575, 10.27645
  };
  constexpr double finite_difference_spacing = 1.0;
  constexpr double transition_covariance_scale = 0.001;

  /**
   * @brief Builds a series of the given length by cycling the observation
   * table, so that an empty series and 1 longer than the table are both
   * reachable.
   */
  std::vector<double> seriesOfLength(const std::size_t length) {
    std::vector<double> series(length);
    for (std::size_t i{0}; i < length; i++) {
      series[i] = observations[i % observations.size()];
    }
    return series;
  }

} // namespace

extern "C" int
LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size) {
  if (size < prefix_size) {
    return 0;
  }
  const std::vector<double> data_series = seriesOfLength(data[0]);
  const std::string dimensions_json(data + prefix_size, data + size);

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
