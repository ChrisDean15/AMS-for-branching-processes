#include <vector>
#include <algorithm>
#include <random>
#include <cmath>
#include <functional>
#include <iostream>

#ifndef BROWNIAN_MOTION_SPHERE_H
#define BROWNIAN_MOTION_SPHERE_H
struct Path_BM {
Path_BM(std::vector<std::vector<double>> Path_vals_,std::vector<double> Score_vals_,double importance_): importance(importance_), Path_vals(Path_vals_), Score_vals(Score_vals_){};
double importance;
std::vector<std::vector<double>> Path_vals;
std::vector<double> Score_vals;
};

class Sampler_BM {
public:
Sampler_BM(double timestep_, double radius_, double target_radius_, std::vector<double> target_location_): 
timestep(timestep_), domain_radius(radius_), target_radius(target_radius_), target_location(target_location_){};
void Resample(Path_BM &path, Path_BM &path_clone, std::mt19937 &rng){
    std::normal_distribution<double> normal_dist(0.0, std::sqrt(timestep));
    if (path.importance == -1) {
        std::vector<double> temp_path = {0,0,0}; 
        while (euclidean_distance(temp_path) < domain_radius && path.importance < 1) {
            path.Path_vals.push_back(temp_path);
            path.Score_vals.push_back(Scorefunction(temp_path));
            path.importance = fmax(path.Score_vals.back(),path.importance);
            for (double &v : temp_path) {
                v += normal_dist(rng);
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
        while (++counter < path.Path_vals.size() && euclidean_distance(temp_path) < domain_radius && path.importance < 1) {
            path.Path_vals[counter]=temp_path;
            path.Score_vals[counter]=Scorefunction(temp_path);
            path.importance = fmax(path.Score_vals[counter],path.importance);
            for (double &v : temp_path) {
                v += normal_dist(rng);
            }
        }
        if (counter < path.Path_vals.size()) {
            path.Path_vals[counter]=temp_path;
            path.Score_vals[counter]=-1;
        } 
        else {
            while(euclidean_distance(temp_path) < domain_radius && path.importance < 1) {
                path.Path_vals.push_back(temp_path);
                path.Score_vals.push_back(Scorefunction(temp_path));
                path.importance = fmax(path.Score_vals.back(),path.importance);
                for (double &v : temp_path) {
                    v += normal_dist(rng);
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
    double temp = fmax(euclidean_distance(point, target_location)-target_radius, 0.0);
    return fmax(0.0, 1.0 - (temp / (2*domain_radius)));

}

double timestep, domain_radius, target_radius;
std::vector<double> target_location;
};
#endif