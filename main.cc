#include <vector>
#include <algorithm>
#include <random>
#include <cmath>
#include <functional>
#include <iostream>
#include "Brownian_motion_sphere.h"
#include "AdaptiveMultilevelSplitting.h"
#include "Brownian_motion_killed.h"
#include "Branching_PDMP_Time_dependent.h"

int main(int argc, char const *argv[])
{
    double timestep = 0.05;
    double radius = 5.0;
    double target_radius = 1;
    std::vector<double> target_location = {0, 0, 10};
    Sampler_BM sampler(timestep, radius, target_radius, target_location);
    Sampler_BMK sampler_killed(timestep, radius);
    Sampler_TDPM sampler_TDPM(timestep, radius, {0,0,0}, 5, 0.95);
    int N = 10000;
    AdaptiveMultilevelSplitting<Path_TDPM, Sampler_TDPM> ams(N, sampler_TDPM); 
    double probability_estimate = ams.run();
    std::cout << "Estimated probability: " << probability_estimate << std::endl;
    return 0;
}
