#include <vector>
#include <algorithm>
#include <random>
#include <cmath>
#include <functional>
#include <iostream>

#ifndef BRANCHING_TDPDMP_H
#define BRANCHING_TDPDMP_H
struct Path_TDPM {
Path_TDPM(std::vector<std::vector<std::vector<double>>> Path_vals_,std::vector<double> Score_vals_,double importance_): importance(importance_), Path_vals(Path_vals_), Score_vals(Score_vals_){};
double importance;
std::vector<std::vector<std::vector<double>>> Path_vals;
std::vector<double> Score_vals;
};

class Sampler_TDPM {
public:
Sampler_TDPM(double timestep_, double radius_, std::vector<double> initial_location_, double branching_rate_, double kill_prob_): INF_temp(std::numeric_limits<double>::infinity()), kill_prob(kill_prob_), initial_location(initial_location_), branching_rate(branching_rate_),
timestep(timestep_), domain_radius(radius_) {};

void Resample(Path_TDPM &path, Path_TDPM &path_clone, std::mt19937 &rng){
    std::normal_distribution<double> normal_dist(0.0, std::sqrt(timestep));
    if (path.importance == -1) {
        double x = normal_dist(rng);
        double y = normal_dist(rng);
        double z = normal_dist(rng);
        double dist = sqrt(x*x+y*y+z*z);
        std::vector<std::vector<double>> temp_path = {{initial_location[0],initial_location[1],initial_location[3],x/dist,y/dist,z/dist}}; 
        while (!temp_path.empty() && path.importance != INF_temp) {
                path.Path_vals.push_back(temp_path);
                path.Score_vals.push_back(Scorefunction(temp_path));
                path.importance = fmax(path.Score_vals.back(),path.importance);
                Path_update(temp_path,rng);   
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
        std::vector<std::vector<double>> temp_path = path.Path_vals[counter];
        Path_update(temp_path,rng);
        while (++counter < path.Path_vals.size() && !temp_path.empty() && path.importance != INF_temp) {
                path.Path_vals[counter]=temp_path;
                path.Score_vals[counter]=Scorefunction(temp_path);
                path.importance = fmax(path.Score_vals[counter],path.importance);
                Path_update(temp_path,rng);
            }   
        if (counter < path.Path_vals.size()) {
            path.Path_vals[counter]=temp_path;
            path.Score_vals[counter]=-1;
        } 
        else {
            while(!temp_path.empty() && path.importance != INF_temp) {
                path.Path_vals.push_back(temp_path);
                path.Score_vals.push_back(Scorefunction(temp_path));
                path.importance = fmax(path.Score_vals.back(),path.importance);
                Path_update(temp_path,rng);
            }
            path.Path_vals.push_back(temp_path);
            path.Score_vals.push_back(-1);
        }
    }
}
void Path_update(std::vector<std::vector<double>> &temp_path, std::mt19937 &rng) {
    std::vector<std::vector<double>> new_paths;
    std::exponential_distribution<double> exp_dist(branching_rate);
    for (auto v: temp_path) {
        if (exp_dist(rng) < timestep) {
            Branch_event(v,new_paths,rng);
        }
        else {
            new_paths.push_back({v[0]+v[3]*timestep,v[1]+v[4]*timestep,v[2]+v[5]*timestep,v[3],v[4],v[5]});
        }

    }
    temp_path = new_paths;
}
void Branch_event(std::vector<double> v, std::vector<std::vector<double>> &new_paths, std::mt19937 &rng) {
    std::normal_distribution<double> normal_dist(0.0, 1);
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    if (dist(rng) < kill_prob) {
    } else {
            for (int i=0; i<2; i++) {
                double x = normal_dist(rng);
                double y = normal_dist(rng);
                double z = normal_dist(rng);
                double dist = sqrt(x*x+y*y+z*z);
                new_paths.push_back({v[0]+v[3]*timestep,v[1]+v[4]*timestep,v[2]+v[5]*timestep,x/dist,y/dist,z/dist});
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
double Scorefunction(std::vector<std::vector<double>> point) {
    double score = 0.0;
    for (auto v: point) {
        if (euclidean_distance({v[0],v[1],v[2]}) >= domain_radius) {
            score = INF_temp;
            break;
        }
        score += euclidean_distance({v[0],v[1],v[2]});
    }
    return score;

}
double timestep, domain_radius, branching_rate, kill_prob, INF_temp;
std::vector<double> initial_location;
};
#endif