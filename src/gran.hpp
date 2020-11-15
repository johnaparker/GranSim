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

class GranSim {
    public:
        GranSim(py_arr position, py_arr rradii, 
                py_arr _mass, double young_mod, double friction, 
                double damp_normal, double damp_tangent, double dt,
                py_arr vposition, py_arr vradii);

        void step();
        void update_position(const Matrix& new_position);
        py::array_t<double> get_position() {
            auto result = py::array_t<double>({Nparticles,2});
            auto r = result.mutable_unchecked<2>();
            for (int i=0; i<Nparticles; i++) {
                r(i,0) = position[i](0);
                r(i,1) = position[i](1);
            }
            return result;
        };

    private:
        void predict();
        void correct();
        void compute_force();
        void assign_voxels();

    public:
        std::vector<vec2> position, velocity;
        int Nparticles;
        std::vector<double> radii;
        double time;
        double dt;

    private:
        std::vector<double> mass;
        double young_mod;
        double friction;
        double damp_normal;
        double damp_tangent;

        std::vector<vec2> rd2, rd3, rd4;
        std::vector<vec2> force;

        int Rparticles, Vparticles;

        std::unordered_map<key_tt, std::vector<int>, key_hash> voxels;
        double voxel_size;
        std::vector<key_tt> voxel_idx;
        std::mutex mtx;
};

#endif
