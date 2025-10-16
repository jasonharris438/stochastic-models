#include "stochastic_models/likelihood/likelihood.h"
LikelihoodCalculator::LikelihoodCalculator() {
  components = nullptr;
};
LikelihoodCalculator::~LikelihoodCalculator() {
  delete components;
  components = nullptr;
};
