#include "stochastic_models/sde/ornstein_uhlenbeck.h"

#include <gtest/gtest.h>

#include <cstdlib>

#include "stochastic_models/hitting_times/hitting_time_ornstein_uhlenbeck.h"
#include "stochastic_models/numeric_utils/helpers.h"
/**
 * @test Tests the output of the
 * OrnsteinUhlenbeckModel::getUnconditionalVariance method and asserts that it
 * is near the expected value.
 *
 */
TEST(OrnsteinUhlenbeckModelTest, getUnconditionalVarianceOutputTest) {
    const float tolerance = 1e-5;
    // Instantiate model and generate output of hitting time core function.
    OrnsteinUhlenbeckModel *model = new OrnsteinUhlenbeckModel(0.5, 0.02, 0.05);
    const double output = model->getUnconditionalVariance();
    delete model;

    // Assert that the value is near the expected value.
    EXPECT_LE(abs(roundToDecimals(output, 8) - 0.0625), tolerance)
        << "OrnsteinUhlenbeckModel not calculating "
           "correct value for getUnconditionalVariance method.";
}
/**
 * @test Tests the output of the OrnsteinUhlenbeckModel::getMean
 * method and asserts that it is near the expected value.
 *
 */
TEST(OrnsteinUhlenbeckModelTest, getMeanOutputTest) {
    const float tolerance = 1e-5;
    // Instantiate model and generate output of hitting time core function.
    OrnsteinUhlenbeckModel *model = new OrnsteinUhlenbeckModel(0.5, 0.02, 0.05);
    const double output = model->getMean();
    delete model;

    // Assert that the value is near the expected value.
    EXPECT_LE(abs(roundToDecimals(output, 8) - 0.5), tolerance)
        << "OrnsteinUhlenbeckModel not calculating "
           "correct value for getMean method.";
}
/**
 * @test Tests the output of the
 * HittingTimeOrnsteinUhlenbeck::hittingTimeDensityCore method and asserts that
 * it is near the expected value.
 *
 */
TEST(HittingTimeOrnsteinUhlenbeckTest, hittingTimeCoreOutputTest) {
    const float tolerance = 1e-5;
    // Instantiate model and generate output of hitting time core function.
    HittingTimeOrnsteinUhlenbeck *hitting_time_kernel =
        new HittingTimeOrnsteinUhlenbeck(0.5, 0.02, 0.05);
    const double output = hitting_time_kernel->hittingTimeDensityCore(0.3);
    delete hitting_time_kernel;

    // Assert that the value is near the expected value.
    EXPECT_LE(abs(roundToDecimals(output, 8) - 0.18637397), tolerance)
        << "OrnsteinUhlenbeckModel not calculating "
           "correct value for hitting time density.";
}
/**
 * @test Tests the output of the
 * HittingTimeOrnsteinUhlenbeck::optimalTradingFCore method and asserts that it
 * is near the expected value.
 *
 */
TEST(HittingTimeOrnsteinUhlenbeckTest, optimalTradingFCoreOutputTest) {
    const float tolerance = 1e-5;
    // Instantiate model and generate output of first passage time core
    // function F.
    HittingTimeOrnsteinUhlenbeck *hitting_time_kernel =
        new HittingTimeOrnsteinUhlenbeck(0.5, 0.02, 0.05);
    const double output =
        hitting_time_kernel->optimalTradingFCore(0.3, 0.1, 0.02);
    delete hitting_time_kernel;
    // Expect equality for output.
    EXPECT_LE(abs(roundToDecimals(output, 8) - 0.91851228), tolerance)
        << "HittingTimeOrnsteinUhlenbeck not calculating "
           "correct value for first passage time F function.";
}
/**
 * @test Tests the output of the
 * HittingTimeOrnsteinUhlenbeck::optimalTradingGCore method and asserts that it
 * is near the expected value.
 *
 */
TEST(HittingTimeOrnsteinUhlenbeckTest, optimalTradingGCoreOutputTest) {
    const float tolerance = 1e-5;
    // Instantiate model and generate output of first passage time core
    // function G.
    HittingTimeOrnsteinUhlenbeck *hitting_time_kernel =
        new HittingTimeOrnsteinUhlenbeck(0.5, 0.02, 0.05);
    const double output =
        hitting_time_kernel->optimalTradingGCore(0.3, 0.1, 0.02);
    delete hitting_time_kernel;
    // Expect equality for output.
    EXPECT_LE(abs(roundToDecimals(output, 8) - 1.07788415), tolerance)
        << "HittingTimeOrnsteinUhlenbeck not calculating "
           "correct value for first passage time G function.";
}
/**
 * @test Tests the output of the
 * HittingTimeOrnsteinUhlenbeck::optimalTradingLCore method and asserts that it
 * is near the expected value.
 *
 */
TEST(HittingTimeOrnsteinUhlenbeckTest, optimalTradingLCoreOutputTest) {
    const float tolerance = 1e-5;
    // Instantiate model and generate output of core function L.
    HittingTimeOrnsteinUhlenbeck *hitting_time_kernel =
        new HittingTimeOrnsteinUhlenbeck(0.3, 8, 0.3);
    const double output = hitting_time_kernel->optimalTradingLCore(0.05, 0.02);
    delete hitting_time_kernel;
    // Expect equality for output.
    EXPECT_LE(abs(roundToDecimals(output, 8) - 0.298261), tolerance)
        << "HittingTimeOrnsteinUhlenbeck not calculating correct value for L "
           "function.";
}
