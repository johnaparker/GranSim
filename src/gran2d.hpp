#ifndef GUARD_gran_2d_h
#define GUARD_gran_2d_h

#include "vec.hpp"
#include <unordered_map>
#include <tuple>
#include <vector>
#include <mutex>
#include <iostream>

#include "geometry_2d.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

typedef std::tuple<int, int> key_tt;

struct key_hash : public std::unary_function<key_tt, std::size_t> {
    std::size_t operator()(const key_tt& k) const {
        //return std::get<0>(k)*100 + std::get<1>(k);
        size_t h = (size_t(std::get<0>(k))<<32) + size_t(std::get<1>(k));
        h *= 1231231557ull;
        h ^= (h>>32); 
        return h; 
    }
};

class granular_media_2d {
    public:
        granular_media_2d(double dt);

        void step();

        void add_wall(vec2 point, vec2 normal);
        void add_grains(py_arr position, py_arr radii, py_arr mass, py_arr young_mod, py_arr friction, py_arr damp_normal, py_arr damp_tangent);
        void add_static_grains(py_arr position, py_arr radii, py_arr mass, py_arr young_mod, py_arr friction, py_arr damp_normal, py_arr damp_tangent);
        //circle_collection add_grains(py_arr position, py_arr radii, py_arr mass, py_arr young_mod, py_arr friction, py_arr damp_normal, py_arr damp_tangent);
        //circle_collection add_static_circles(py_arr position, py_arr radii, py_arr mass, py_arr young_mod, py_arr friction, py_arr damp_normal, py_arr damp_tangent);
        //void update_position(const Matrix& new_position);

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
