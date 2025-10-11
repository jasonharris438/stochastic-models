#include "stochastic_models/trading/trading_levels_params.h"
ExitLevelStopLossParams::~ExitLevelStopLossParams() {
  delete optimizer;
  optimizer = nullptr;
  delete hitting_time_kernel;
  hitting_time_kernel = nullptr;
}
ExitLevelParams::~ExitLevelParams() {
  delete optimizer;
  optimizer = nullptr;
  delete hitting_time_kernel;
  hitting_time_kernel = nullptr;
}
EntryLevelStopLossParams::~EntryLevelStopLossParams() {
  delete optimizer;
  optimizer = nullptr;
  delete hitting_time_kernel;
  hitting_time_kernel = nullptr;
}
EntryLevelParams::~EntryLevelParams() {
  delete optimizer;
  optimizer = nullptr;
  delete hitting_time_kernel;
  hitting_time_kernel = nullptr;
}
double funcOptimalMeanReversionA(double x, void* params) {
  struct EntryLevelParams* p = static_cast<EntryLevelParams*>(params);
  return p->optimizer->a(x, p->hitting_time_kernel, p->b_star, p->r, p->c);
}
double funcOptimalMeanReversionStopLossA(double x, void* params) {
  struct EntryLevelStopLossParams* p =
      static_cast<EntryLevelStopLossParams*>(params);
  return p->optimizer->a(
      x, p->hitting_time_kernel, p->b_star, p->stop_loss, p->r, p->c
  );
}
double funcOptimalMeanReversionStopLossD(double x, void* params) {
  struct EntryLevelStopLossParams* p =
      static_cast<EntryLevelStopLossParams*>(params);
  return p->optimizer->d(
      x, p->hitting_time_kernel, p->b_star, p->stop_loss, p->r, p->c
  );
}
double funcOptimalMeanReversionStopLossB(double x, void* params) {
  struct ExitLevelStopLossParams* p =
      static_cast<ExitLevelStopLossParams*>(params);
  return p->optimizer->b(x, p->hitting_time_kernel, p->stop_loss, p->r, p->c);
}
double funcOptimalMeanReversionD(double x, void* params) {
  struct EntryLevelParams* p = static_cast<EntryLevelParams*>(params);
  return p->optimizer->d(x, p->hitting_time_kernel, p->b_star, p->r, p->c);
}
double funcOptimalMeanReversionB(double x, void* params) {
  struct ExitLevelParams* p = static_cast<ExitLevelParams*>(params);
  return p->optimizer->b(x, p->hitting_time_kernel, p->r, p->c);
}
