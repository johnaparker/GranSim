#ifndef GUARD_gran_h
#define GUARD_gran_h

#include "vec.hpp"
using Eigen::Ref;

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

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

    public:
        Matrix position, velocity;

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
};

#endif
