#include "stochastic_models/distributions/gaussian.h"

#include <cstdint>
#include <gtest/gtest.h>

/**
 * @test Two instances built with the same explicit seed produce identical
 * sequences (reproducibility).
 */
TEST(GaussianSeedingTest, SameSeedIsReproducible) {
  GaussianDistribution a(0.0, 1.0, std::uint64_t{12345});
  GaussianDistribution b(0.0, 1.0, std::uint64_t{12345});
  EXPECT_EQ(a.sample(50), b.sample(50));
}

/**
 * @test Different seeds produce different sequences.
 */
TEST(GaussianSeedingTest, DifferentSeedsDiffer) {
  GaussianDistribution a(0.0, 1.0, std::uint64_t{1});
  GaussianDistribution b(0.0, 1.0, std::uint64_t{2});
  EXPECT_NE(a.sample(50), b.sample(50));
}

/**
 * @test Successive sample() calls on one instance advance the engine rather
 * than restarting it (pre-fix code re-seeded each call).
 */
TEST(GaussianSeedingTest, EngineAdvancesAcrossCalls) {
  GaussianDistribution a(0.0, 1.0, std::uint64_t{999});
  const std::vector<double> first = a.sample(50);
  const std::vector<double> second = a.sample(50);
  EXPECT_NE(first, second);
}

/**
 * @test Default-constructed instances are seeded nondeterministically, so two
 * of them produce different sequences (pre-fix: identical within the same
 * wall-clock second).
 */
TEST(GaussianSeedingTest, DefaultConstructedInstancesDiffer) {
  GaussianDistribution a;
  GaussianDistribution b;
  EXPECT_NE(a.sample(50), b.sample(50));
}
