/// @brief Parameters to be passed to the integrand function.
struct params {
    double alpha;
    double mu;
    double sigma;
};
/// The solution to the backward equation, the scale function @f[ S(x) @f]
/// @param x Point at which to evaluate the scale function.
/// @param p The parameters used in the scale function.
/// @return The value of the scale function at @f[ x @f].
double f(double x, void *p);
double integration();
