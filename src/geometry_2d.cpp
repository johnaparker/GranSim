#include "geometry_2d.hpp"

Circle::Circle(vec2 position, double radius, double mass, double young_mod, double friction, double damp_normal, double damp_tangent): position(position), radius(radius), mass(mass), young_mod(young_mod), friction(friction), damp_normal(damp_normal), damp_tangent(damp_tangent) {

    velocity = vec2::Zero();
    rd2 = vec2::Zero();
    rd3 = vec2::Zero();
    rd4 = vec2::Zero();
    force = vec2::Zero();
}

Wall2d::Wall2d(vec2 point, vec2 normal): point(point), normal(normal) {
    normal.normalize();
    tangent = vec2(normal(1), normal(0));
}

void interact(Circle& c1, Circle& c2) {
    vec2 dr = c1.position - c2.position;
    double dr_norm = dr.norm();
    double overlap = c1.radius + c2.radius - dr_norm;
    dr /= dr_norm;
    vec2 dt(-dr(1), dr(0));
    vec2 dv = c1.velocity - c2.velocity;

    double dv_n = -dr.dot(dv);
    double dv_t = dt.dot(dv);

    double reff = c1.radius*c2.radius/(c1.radius + c2.radius);
    double young_mod = c1.young_mod*c2.young_mod/(c1.young_mod + c2.young_mod);
    double damp_normal = .5*(c1.damp_normal + c2.damp_normal);
    double damp_tangent = std::min(c1.damp_tangent, c2.damp_tangent);
    double friction = std::min(c1.friction, c2.friction);

    double Fn_mag = std::max(0.0, sqrt(reff)*young_mod*sqrt(overlap)*(overlap + damp_normal*dv_n));
    double Ft_mag = std::min(friction*Fn_mag, damp_tangent*std::abs(dv_t));

    auto F = dr*Fn_mag - dt*sgn(dv_t)*Ft_mag;

    c1.force += F;
    c2.force -= F;
}

void interact(Circle& c, const Wall2d& w) {
    double overlap = c.radius - (c.position - w.point).dot(w.normal);
    if (overlap > 0) {
        double dv_n = -c.velocity.dot(w.normal);
        double dv_t = c.velocity.dot(w.tangent);
        double reff = c.radius;
        double Fn_mag = std::max(0.0, sqrt(reff)*c.young_mod*sqrt(overlap)*(overlap + c.damp_normal*dv_n));
        double Ft_mag = std::min(c.friction*Fn_mag, c.damp_tangent*std::abs(dv_t));
        c.force += w.normal*Fn_mag - w.tangent*Ft_mag*sgn(dv_t);
    }
}
