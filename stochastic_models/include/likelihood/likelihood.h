#ifndef _LIKELIHOOD_H_
#define _LIKELIHOOD_H_
#include <map>
#include <string>
#include <vector>

/**
 * Base likelihood calculator class.
 */
class LikelihoodCalculatuor {
   protected:
    const std::map<std::string, double> coefs;
    LikelihoodCalculatuor();

   public:
    /**
     * Calculates and returns the maximised likelihood parameters.
     *
     * @param data The data series to calculate the likelihood of.
     * @returns Const reference to mean of Gaussian Distribution.
     */
    virtual const std::map<std::string, double>& calculate(
        const std::vector<double>& data) const = 0;
};

#endif  // _LIKELIHOOD_H_
