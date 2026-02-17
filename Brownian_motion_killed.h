#include <vector>
#include <algorithm>
#include <random>
#include <cmath>
#include <functional>
#include <iostream>

#ifndef BROWNIAN_MOTION_KILLED_H
#define BROWNIAN_MOTION_KILLED_H
struct Path_BMK {
Path_BMK(std::vector<std::vector<double>> Path_vals_,std::vector<double> Score_vals_,double importance_): importance(importance_), Path_vals(Path_vals_), Score_vals(Score_vals_){};
double importance;
std::vector<std::vector<double>> Path_vals;
std::vector<double> Score_vals;
};

class Sampler_BMK {
public:
Sampler_BMK(double timestep_, double radius_): 
timestep(timestep_), domain_radius(radius_) {};
void Resample(Path_BMK &path, Path_BMK &path_clone, std::mt19937 &rng){
    std::normal_distribution<double> normal_dist(0.0, std::sqrt(timestep));
    std::exponential_distribution<double> exp_dist(1.0);
    if (path.importance == -1) {
        std::vector<double> temp_path = {0,0,0}; 
        while (euclidean_distance(temp_path) < domain_radius+5 && path.importance < 1) {
                path.Path_vals.push_back(temp_path);
                path.Score_vals.push_back(Scorefunction(temp_path));
                path.importance = fmax(path.Score_vals.back(),path.importance);
                for (double &v : temp_path) {
                    v += normal_dist(rng);
                }
            if (exp_dist(rng) < timestep) {
                temp_path = {100,100,100};
                path.Path_vals.push_back(temp_path);
                path.Score_vals.push_back(-1);
            } 
            
        }
        path.Path_vals.push_back(temp_path);
        path.Score_vals.push_back(-1);
    }
    else {
        int counter = -1;
        path.Path_vals = path_clone.Path_vals;
        path.Score_vals = path_clone.Score_vals;
        do {
            counter++;
        } while(path.Score_vals[counter]<=path.importance);
        path.importance = path.Score_vals[counter];
        std::vector<double> temp_path = path.Path_vals[counter];
        for (double &v : temp_path) {
                v += normal_dist(rng);
            }
        while (++counter < path.Path_vals.size() && euclidean_distance(temp_path) < domain_radius+5 && path.importance < 1) {
            if (exp_dist(rng) < timestep) {
                temp_path = {100,100,100};
                path.Path_vals[counter]=temp_path;
                path.Score_vals[counter]=-1;
            } 
            else {
                path.Path_vals[counter]=temp_path;
                path.Score_vals[counter]=Scorefunction(temp_path);
                path.importance = fmax(path.Score_vals[counter],path.importance);
                for (double &v : temp_path) {
                    v += normal_dist(rng);
                }
            }   
        }
        if (counter < path.Path_vals.size()) {
            path.Path_vals[counter]=temp_path;
            path.Score_vals[counter]=-1;
        } 
        else {
            while(euclidean_distance(temp_path) < domain_radius+5 && path.importance < 1) {
                if (exp_dist(rng) < timestep) {
                temp_path = {100,100,100};
                path.Path_vals.push_back(temp_path);
                path.Score_vals.push_back(-1);
                }
                else { 
                path.Path_vals.push_back(temp_path);
                path.Score_vals.push_back(Scorefunction(temp_path));
                path.importance = fmax(path.Score_vals.back(),path.importance);
                for (double &v : temp_path) {
                    v += normal_dist(rng);
                }
            }
            }
            path.Path_vals.push_back(temp_path);
            path.Score_vals.push_back(-1);
        }
    }
}
double euclidean_distance(std::vector<double> vec_1, std::vector<double> vec_2={0,0,0}) {
    double sum = 0.0;
    for (size_t i=0; i<vec_1.size(); i++) {
        double diff = vec_1[i] - vec_2[i];
        sum += diff * diff;
    }
    return std::sqrt(sum);
}
double Scorefunction(std::vector<double> point) {
    return fmin(1.0, euclidean_distance(point) / domain_radius);

}

double timestep, domain_radius;
};
#endif