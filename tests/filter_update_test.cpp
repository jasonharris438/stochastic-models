
#include <gtest/gtest.h>

#include "stochastic_models/kalman_filter/kca.h"
/**
 * @brief Test that the KineticComponents initialiseFilter sets a
 * KinetiComponents instance to the correct initial state. Members are private,
 * however logically the method must set the isInitialised flag to true. We can
 * only test the public interface.
 */
TEST(KalmanFilterUpdateTest, KineticComponentsinitialiseFilterTest) {
    // Mock data series to initialise the KineticComponents object.
    const std::vector<double> data_series{
        10.51255, 10.51985, 10.52405, 10.4656,  10.47,    10.5403,  10.4425,
        10.3087,  10.1994,  10.1839,  10.24645, 10.1795,  10.21715, 10.14995,
        10.194,   10.22505, 10.27325, 10.25095, 10.30575, 10.27645};
    const double h{1.0};
    const double q{0.001};

    // Create the KineticComponents object and check the initial state is set
    // correctly.
    const FilterSystemDimensions dimensions(3, 3, 3, 1, 3, 1, 1, 0.0);
    KineticComponents kinetic_components(dimensions);
    kinetic_components.initialiseFilter(data_series, h, q);
    EXPECT_TRUE(kinetic_components.isInitialised())
        << "The KineticComponents object must correctly indicate whether it is "
           "initialised.";

    std::vector<double> expected_current_state_mean{10.288741828687053, 0.0,
                                                    0.0};
    EXPECT_EQ(kinetic_components.getCurrentState(), expected_current_state_mean)
        << "The KineticComponents object must correctly set the current state "
           "mean at initialisation.";
}
/**
 * @brief Test that the KineticComponents updatePriors sets a KineticComponents
 * instance to the correct state with valid prior values. Members are private,
 * however logically the method must set the isPriorStateValid flag to true. We
 * can only test the public interface.
 */
TEST(KalmanFilterUpdateTest, KineticComponentsupdatePriorsTest) {
    // Mock data series to initialise a KCA states object.
    std::vector<std::vector<double>> transition_matrix{
        {1.0011961162353782, 1.0, 0.5}, {0.0, 1.0, 1.0}, {0.0, 0.0, 1.0}};
    std::vector<std::vector<double>> transition_covariance{
        {0.12695229227341848, 0, 0}, {0, 0.001, 0}, {0, 0, 0.001}};
    std::vector<std::vector<double>> current_state_covariance{
        {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    std::vector<double> current_state_mean{10.288741828687053, 0.0, 0.0};
    std::vector<std::vector<double>> observation_matrix{{1.0, 0.0, 0.0}};
    const double observation_offset{0.0};

    // Create the KCA states object and set the initial state flag.
    const FilterSystemDimensions dimensions(3, 3, 3, 1, 3, 1, 1, 0.0);
    KcaStates kca_states(dimensions);
    kca_states.setTransitionMatrix(transition_matrix);
    kca_states.setTransitionCovariance(transition_covariance);
    kca_states.setCurrentStateCovariance(current_state_covariance);
    kca_states.setCurrentStateMean(current_state_mean);
    kca_states.setObservationMatrix(observation_matrix);
    kca_states.setObservationOffset(observation_offset);
    kca_states.setInitialized();

    // Create the KineticComponents object and set the KcaState object as its
    // filter_state.
    KineticComponents kinetic_components(dimensions);
    kinetic_components.setFilterState(kca_states);

    // Update the priors and check the prior state is valid.
    kinetic_components.updatePriors();
    EXPECT_TRUE(kinetic_components.isPriorStateValid())
        << "The KineticComponents object must correctly indicate whether its "
           "prior state is valid after prior update.";

    std::vector<double> expected_current_state_mean{10.288741828687053, 0.0,
                                                    0.0};
    EXPECT_EQ(kinetic_components.getCurrentState(), expected_current_state_mean)
        << "The KineticComponents object must not change the value of the the "
           "current state when updating the prior state.";
}

/**
 * @brief Test that the KineticComponents updatePosteriors sets a
 * KineticComponents instance to the correct state and returns a the corrent
 * state objects.
 */
TEST(KalmanFilterUpdateTest, KineticComponentsupdatePosteriorsTest) {
    // Mock data series to initialise a KCA states object.
    std::vector<std::vector<double>> transition_matrix{
        {1.0011961162353782, 1.0, 0.5}, {0.0, 1.0, 1.0}, {0.0, 0.0, 1.0}};
    std::vector<std::vector<double>> transition_covariance{
        {0.12695229227341848, 0, 0}, {0, 0.001, 0}, {0, 0, 0.001}};
    std::vector<std::vector<double>> current_state_covariance{
        {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    std::vector<double> current_state_mean{10.288741828687053, 0.0, 0.0};
    std::vector<std::vector<double>> observation_matrix{{1.0, 0.0, 0.0}};
    const double observation_offset{0.0};
    std::vector<double> predicted_state_mean{10.301048359829961, 0.0, 0.0};
    std::vector<std::vector<double>> predicted_state_covariance{
        {0.12695229227341848, 0, 0}, {0, 0.001, 0}, {0, 0, 0.001}};

    // Create the KCA states object and check the initial state is set
    // correctly.
    const FilterSystemDimensions dimensions(3, 3, 3, 1, 3, 1, 1, 0.0);
    KcaStates kca_states(dimensions);
    kca_states.setTransitionMatrix(transition_matrix);
    kca_states.setTransitionCovariance(transition_covariance);
    kca_states.setCurrentStateCovariance(current_state_covariance);
    kca_states.setCurrentStateMean(current_state_mean);
    kca_states.setObservationMatrix(observation_matrix);
    kca_states.setObservationOffset(observation_offset);
    kca_states.setPredictedStateMean(predicted_state_mean);
    kca_states.setPredictedStateCovariance(predicted_state_covariance);
    kca_states.setInitialized();
    kca_states.setPriorsTrue();

    // Create the KineticComponents object and set the KcaState object as its
    // filter_state.
    KineticComponents kinetic_components(dimensions);
    kinetic_components.setFilterState(kca_states);

    // Update the priors and check the prior state is valid.
    const double observation = 10.3;
    const double innovation_sigma = 0.1;
    kinetic_components.updatePosteriors(observation, innovation_sigma);
    EXPECT_FALSE(kinetic_components.isPriorStateValid())
        << "The KineticComponents object must correctly indicate whether its "
           "prior state is invalid after posterior update.";

    std::vector<double> expected_current_state_mean{10.3000765492722, 0.0, 0.0};
    EXPECT_EQ(kinetic_components.getCurrentState(), expected_current_state_mean)
        << "The KineticComponents object must correctly set the current state "
           "mean at initialisation.";
}
/**
 * @brief Test that the KineticComponents correctly completes a full initialize,
 * predict and update round with the correct current state values.
 */
TEST(KalmanFilterUpdateTest, KineticComponentsFullFilterTest) {
    // Mock data series to initialise the KineticComponents object.
    const std::vector<double> data_series{
        10.51255, 10.51985, 10.52405, 10.4656,  10.47,    10.5403,  10.4425,
        10.3087,  10.1994,  10.1839,  10.24645, 10.1795,  10.21715, 10.14995,
        10.194,   10.22505, 10.27325, 10.25095, 10.30575, 10.27645};
    const double h{1.0};
    const double q{0.001};
    const double observation = 10.3;
    const double innovation_sigma = 0.1;

    // Create the KineticComponents object.
    const FilterSystemDimensions dimensions(3, 3, 3, 1, 3, 1, 1, 0.0);
    KineticComponents kinetic_components(dimensions);
    // Initialise the filter with the data series.
    kinetic_components.initialiseFilter(data_series, h, q);
    // Update the priors and make prediction.
    kinetic_components.updatePriors();
    // Update the posteriors with the observation.
    kinetic_components.updatePosteriors(observation, innovation_sigma);

    std::vector<double> expected_current_state_mean{10.3000765492722, 0.0, 0.0};
    EXPECT_EQ(kinetic_components.getCurrentState(), expected_current_state_mean)
        << "The KineticComponents object must correctly finish a kalman filter "
           "predict and update round with the correct current state.";
}
