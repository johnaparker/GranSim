#ifndef GUARD_gran_h
#define GUARD_gran_h

#include "vec.hpp"
#include <unordered_map>
#include <tuple>
#include <vector>

using Eigen::Ref;
typedef std::tuple<int, int> key_t;

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

struct key_hash : public std::unary_function<key_t, std::size_t> {
    std::size_t operator()(const key_t& k) const {
        //return std::get<0>(k)*100 + std::get<1>(k);
        size_t h = (size_t(std::get<0>(k))<<32) + size_t(std::get<1>(k));
        h *= 1231231557ull;
        h ^= (h>>32); 
        return h; 
    }
};

class GranSim {
    public:
        GranSim(const Matrix& position, const Array& radii, 
                const Array& mass, double young_mod, double friction, 
                double damp_normal, double damp_tangent, double dt);

        void step();

    private:
        void predict();
        void correct();
        void compute_force();
        void assign_voxels();

    public:
        Matrix position, velocity;

    private:
        Array radii;
        Array mass;
        double young_mod;
        double friction;
        double damp_normal;
        double damp_tangent;

        Matrix rd2, rd3, rd4;
        Matrix force;

        double dt;
        double time;
        int Nparticles;

        std::unordered_map<key_t, std::vector<int>, key_hash> voxels;
        double voxel_size;
        std::vector<key_t> voxel_idx;
};

#endif
