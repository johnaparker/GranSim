#ifndef GUARD_gran_3d_h
#define GUARD_gran_3d_h

#include "vec.hpp"
#include <unordered_map>
#include <tuple>
#include <vector>
#include <mutex>
#include <iostream>

#include "geometry_3d.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

typedef std::tuple<int, int, int> key_tt;

struct key_hash : public std::unary_function<key_tt, std::size_t> {
    const int N1 = 100;
    const int N2 = N1*N1;
    std::size_t operator()(const key_tt& k) const {
        return std::get<0>(k) + std::get<1>(k)*N1 + std::get<2>(k)*N2;
        size_t h1 = (size_t(std::get<0>(k))<<32) + size_t(std::get<2>(k));
        size_t h2 = (size_t(std::get<0>(k))<<32) + size_t(std::get<1>(k));
        h1 *= 1231231557ull;
        h1 ^= (h1>>32); 
        h2 *= 1231231557ull;
        h2 ^= (h2>>32); 
        return h1+h2; 
    }
};

class granular_media {
    public:
        granular_media(double dt);

        void step();

        void add_wall(vec3 point, vec3 normal);
        void add_grains(py_arr position, py_arr radii, py_arr mass, py_arr young_mod, py_arr friction, py_arr damp_normal, py_arr damp_tangent);
        void add_static_grains(py_arr position, py_arr radii, py_arr mass, py_arr young_mod, py_arr friction, py_arr damp_normal, py_arr damp_tangent);

        py::array_t<double> get_position();

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
        vec3 gravity;

    private:
        std::vector<Sphere> d_grains;
        std::vector<Sphere> s_grains;
        std::vector<Wall> walls;
        int Rparticles, Vparticles;

        std::unordered_map<key_tt, std::vector<int>, key_hash> voxels;
        double voxel_size;
        std::vector<key_tt> voxel_idx;
        std::mutex mtx;
};

#endif
