#include <vector>
#include <algorithm>
#include <random>
#include <cmath>
#include <functional>
#include <iostream>
#include "Brownian_motion_sphere.h"
#include "Brownian_motion_killed.h"
#include "Branching_PDMP_Time_dependent.h"

#ifndef ADAPTIVE_MULTILEVEL_SPLITTING_H
#define ADAPTIVE_MULTILEVEL_SPLITTING_H

template <typename Path, typename Simulator>
class AdaptiveMultilevelSplitting {
public:
    AdaptiveMultilevelSplitting(
        int N, Simulator simulator)
        : N_(N),
          min_score_(1),
          scores_(N, 0.0), simulator_(simulator)
    {   
        std::random_device rd;
        rng_ = std::mt19937(rd());
        particles_=std::vector<Path>(N_,Path({},{},-1));
    }
    double run(double threshold = std::numeric_limits<double>::infinity()) {
        for (int i = 0; i < N_; i++) {
            simulator_.Resample(particles_[i],particles_[i],rng_);
            scores_[i] = particles_[i].importance;
        }

        double log_estimator = 0.0;
        int iterations = 0;

        while (true) {
            // Sort indices by score
            std::vector<int> idx(N_);
            std::iota(idx.begin(), idx.end(), 0);

            std::sort(idx.begin(), idx.end(),
                      [&](int a, int b) { return scores_[a] < scores_[b]; });

            if (scores_[idx[0]] >= threshold)
                break;
            double level = scores_[idx[0]];
            if (iterations % 100 == 0){
                std::cout << "Iteration " << iterations << ", level: " << level << std::endl;
            }
            int k = 0;
            while (k < N_ && scores_[idx[k]] <= level) {
                k++;
            }
            if (k==N_) {
                std::cout<< "All particles have minimal importance, terminating simulation." << std::endl;
                break;
            }
            log_estimator += std::log(1.0 - double(k) / N_);
            iterations++;

            std::uniform_int_distribution<int> pick(k, N_ - 1);

            for (int i = 0; i < k; i++) {
                int parent = idx[pick(rng_)];
                int child = idx[i];
                simulator_.Resample(particles_[child], particles_[parent], rng_);
                scores_[child] = particles_[child].importance;

            }
        }

        return std::exp(log_estimator);
    }

private:
    int N_;
    double min_score_;
    std::vector<double> scores_;
    Simulator simulator_;
    std::mt19937 rng_;
    std::vector<Path> particles_;
};

#endif