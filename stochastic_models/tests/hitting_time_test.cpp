#include <gtest/gtest.h>

#include "numeric_utils/integration.h"
// Demonstrate some basic assertions.
TEST(HittingTimeDensityTest, OutputTest) {
    // Run integration function.
    double value = integration();
    EXPECT_EQ(value, 0.003993143831817661)
        << "The value of the integration is not equal to the expected value.";
}
