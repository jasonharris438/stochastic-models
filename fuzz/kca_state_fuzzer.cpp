#include "kca_fuzz_support.h"
#include "stochastic_models/entrypoints/kca_filter.h"

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <string>

namespace {

  constexpr std::size_t prefix_size = 2 * sizeof(double);

} // namespace

extern "C" int
LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size) {
  if (size < prefix_size) {
    return 0;
  }
  double observation{};
  double innovation_sigma{};
  std::memcpy(&observation, data, sizeof(double));
  std::memcpy(&innovation_sigma, data + sizeof(double), sizeof(double));
  const std::string state(data + prefix_size, data + size);

  std::string output;
  const kca_fuzz::Outcome outcome = kca_fuzz::runClassified([&] {
    output = getUpdatedKcaState(
        state, kca_fuzz::scheme_dimensions, observation, innovation_sigma
    );
  });
  if (outcome == kca_fuzz::Outcome::Completed) {
    kca_fuzz::checkRoundTrip(output, kca_fuzz::scheme_dimensions);
  }
  return 0;
}
