#include "stochastic_models/likelihood/likelihood.h"
LikelihoodCalculatuor::LikelihoodCalculatuor() {
  components = nullptr;
};
LikelihoodCalculatuor::~LikelihoodCalculatuor() {
  delete components;
  components = nullptr;
};
