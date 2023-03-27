/// @brief Parameters to be passed to the integrand function.
struct params {
    double alpha;
    double mu;
    double sigma;
};
/**
 * @brief The solution to the backward equation, the scale function @f[ S(x) @f]
 * @param x Point at which to evaluate the scale function.
 * @param p The parameters used in the scale function.
 * @return double. The value of the scale function at @f[ x @f].
 */
double f(double x, void* p);

/**
 * @brief Integrates the scale function over a given interval.
 *
 * @param alpha The parameter @f[ \alpha @f] in the Ornstein-Uhlenbeck model.
 * @param mu The parameter @f[ \mu @f] in the Ornstein-Uhlenbeck model.
 * @param sigma The parameter @f[ \sigma @f] in the Ornstein-Uhlenbeck model.
 * @param lower The lower bound of the integration interval.
 * @param upper The upper bound of the integration interval.
 * @return const double. The value of the integral of F.
 */
const double integration(double& alpha, double& mu, double& sigma,
                         double& lower, double& upper);
/**
 * @brief The hitting time density function @f[ \frac{S(x) - S(b)}{S(a) - S(b)}
 * @f].
 *
 * @param x Point at which to evaluate the hitting time density function.
 * @param alpha The parameter @f[ \alpha @f] in the Ornstein-Uhlenbeck model.
 * @param mu The parameter @f[ \mu @f] in the Ornstein-Uhlenbeck model.
 * @param sigma The parameter @f[ \sigma @f] in the Ornstein-Uhlenbeck model.
 * @param first The first hitting time.
 * @param second The second hitting time.
 * @return const double. The value of the hitting time density function at @f[
 * x@f].
 */
const double hitting_time_density(double& x, double& alpha, double& mu,
                                  double& sigma, double& first, double& second);