#ifndef GUARD_gran_h
#define GUARD_gran_h

#include "vec.hpp"
#include <unordered_map>
#include <tuple>
#include <vector>
#include <mutex>
#include <iostream>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

using Eigen::Ref;
typedef std::tuple<int, int> key_tt;

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

struct key_hash : public std::unary_function<key_tt, std::size_t> {
    std::size_t operator()(const key_tt& k) const {
        //return std::get<0>(k)*100 + std::get<1>(k);
        size_t h = (size_t(std::get<0>(k))<<32) + size_t(std::get<1>(k));
        h *= 1231231557ull;
        h ^= (h>>32); 
        return h; 
    }
};

class Circle {
    public:
        Circle(vec2 position, double radius, double mass, double young_mod, double friction, double damp_normal, double damp_tangent);


    public:
        vec2 position, velocity, rd2, rd3, rd4, force;
        double radius, mass, young_mod, friction, damp_normal, damp_tangent;
};

class Wall2d {
    public:
        Wall2d(vec2 point, vec2 normal);

    public:
        vec2 point, normal, tangent;
};

void interact(Circle& c1, Circle& c2);
void interact(Circle& c, const Wall2d& w);

class granular_media_2d {
    public:
        granular_media_2d(double dt);

        void step();
        py::array_t<double> get_position() {
            auto result = py::array_t<double>({Nparticles,2});
            auto r = result.mutable_unchecked<2>();
            for (int i=0; i<Nparticles; i++) {
                const auto& grain = (i < Rparticles) ? d_grains[i] : s_grains[i-Rparticles];
                r(i,0) = grain.position(0);
                r(i,1) = grain.position(1);
            }
            return result;
        };

        void add_wall(vec2 point, vec2 normal);
        void add_grains(py_arr position, py_arr radii, py_arr mass, py_arr young_mod, py_arr friction, py_arr damp_normal, py_arr damp_tangent);
        void add_static_grains(py_arr position, py_arr radii, py_arr mass, py_arr young_mod, py_arr friction, py_arr damp_normal, py_arr damp_tangent);
        //circle_collection add_grains(py_arr position, py_arr radii, py_arr mass, py_arr young_mod, py_arr friction, py_arr damp_normal, py_arr damp_tangent);
        //circle_collection add_static_circles(py_arr position, py_arr radii, py_arr mass, py_arr young_mod, py_arr friction, py_arr damp_normal, py_arr damp_tangent);
        //void update_position(const Matrix& new_position);

    private:
        void predict();
        void correct();
        void compute_force();
        void initialize_voxels();
        void assign_voxels();

    public:
        int Nparticles;
        double time;
        double dt;
        vec2 gravity;

    private:
        std::vector<Circle> d_grains;
        std::vector<Circle> s_grains;
        std::vector<Wall2d> walls;
        int Rparticles, Vparticles;

        std::unordered_map<key_tt, std::vector<int>, key_hash> voxels;
        double voxel_size;
        std::vector<key_tt> voxel_idx;
        std::mutex mtx;
};

#endif
