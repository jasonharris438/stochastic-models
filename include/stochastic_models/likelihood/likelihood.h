#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H
#include <string>
#include <unordered_map>
#include <vector>
/**
 * Abstract base likelihood calculator class.
 */
class LikelihoodCalculatuor {
   protected:
    /**
     * @brief Empty pointer to a unordered_map that stores likelihood
     * calculation components used in generating process parameters.
     *
     */
    const std::unordered_map<std::string, const double>* components;

   public:
    LikelihoodCalculatuor();
    /**
     * Calculates and returns the maximised likelihood parameters.
     *
     * @param data The data series to calculate the likelihood of.
     * @returns Const reference to mean of Gaussian Distribution.
     */
    virtual const std::unordered_map<std::string, const double> calculate(
        const std::vector<double>& data) = 0;
    ~LikelihoodCalculatuor();
};
#endif  // LIKELIHOOD_H
