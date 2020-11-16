#ifndef GUARD_geometry_3d_h
#define GUARD_geometry_3d_h

#include "vec.hpp"

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}


class Sphere {
    public:
        Sphere(vec3 position, double radius, double mass, double young_mod, double friction, double damp_normal, double damp_tangent);


    public:
        vec3 position, velocity, rd2, rd3, rd4, force;
        double radius, mass, young_mod, friction, damp_normal, damp_tangent;
};

class Wall {
    public:
        Wall(vec3 point, vec3 normal);

    public:
        vec3 point, normal;
};

void interact(Sphere& c1, Sphere& c2);
void interact(Sphere& c, const Wall& w);

#endif
