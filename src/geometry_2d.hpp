#ifndef GUARD_geometry_2d_h
#define GUARD_geometry_2d_h

#include "vec.hpp"

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}


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

#endif
