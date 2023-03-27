#include "sde/ornstein_uhlenbeck.h"

#include <gtest/gtest.h>
// Demonstrate some basic assertions.
TEST(OrnsteinUhlenbeckModelTest, HittingTimeTest) {
    // Instantiate model and generate output of hitting time core function.
    OrnsteinUhlenbeckModel *model = new OrnsteinUhlenbeckModel(0.5, 0.02, 0.05);
    const double output = model->HittingTimeDensityCore(0.3);
    delete model;
    // Expect equality for output.
    EXPECT_EQ(output, 0.18637397603941006)
        << "OrnsteinUhlenbeckModel not calculating "
           "correct value for hitting time density.";
}
