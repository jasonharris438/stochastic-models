#include "ornstein_uhlenbeck.h"

#include <iostream>

int main()
{
    OrnsteinUhlenbeckModel *model = new OrnsteinUhlenbeckModel(0, 0.0108, 0.0067);
    std::vector<double> vec = (*model).Simulate(1000);

    for(size_t i{0} ; i < 5 ; i++){
    std::cout << vec[i] << std::endl;
    }
    return 0;
}