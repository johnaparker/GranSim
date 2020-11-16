#include "geometry_3d.hpp"

Sphere::Sphere(vec3 position, double radius, double mass, double young_mod, double friction, double damp_normal, double damp_tangent): position(position), radius(radius), mass(mass), young_mod(young_mod), friction(friction), damp_normal(damp_normal), damp_tangent(damp_tangent) {

    velocity = vec3::Zero();
    rd2 = vec3::Zero();
    rd3 = vec3::Zero();
    rd4 = vec3::Zero();
    force = vec3::Zero();
}

Wall::Wall(vec3 point, vec3 normal): point(point), normal(normal) {
    normal.normalize();
}

void interact(Sphere& s1, Sphere& s2) {
    vec3 dr = s1.position - s2.position;
    double dr_norm = dr.norm();
    double overlap = s1.radius + s2.radius - dr_norm;
    dr /= dr_norm;
    vec3 dv = s1.velocity - s2.velocity;

    double dv_n = -dr.dot(dv);
    vec3 dv_t = dv - dv_n*dr;
    double dv_t_norm = dv_t.norm();

    double reff = s1.radius*s2.radius/(s1.radius + s2.radius);
    double young_mod = s1.young_mod*s2.young_mod/(s1.young_mod + s2.young_mod);
    double damp_normal = .5*(s1.damp_normal + s2.damp_normal);
    double damp_tangent = std::min(s1.damp_tangent, s2.damp_tangent);
    double friction = std::min(s1.friction, s2.friction);

    double Fn_mag = std::max(0.0, sqrt(reff)*young_mod*sqrt(overlap)*(overlap + damp_normal*dv_n));
    double Ft_mag = std::min(friction*Fn_mag, damp_tangent*dv_t_norm);
    auto F = dr*Fn_mag - Ft_mag*dv_t/dv_t_norm;

    s1.force += F;
    s2.force -= F;
}

void interact(Sphere& s, const Wall& w) {
    double overlap = s.radius - (s.position - w.point).dot(w.normal);
    if (overlap > 0) {
        double dv_n = -s.velocity.dot(w.normal);
        double reff = s.radius;
        double Fn_mag = std::max(0.0, sqrt(reff)*s.young_mod*sqrt(overlap)*(overlap + s.damp_normal*dv_n));

        vec3 dv_t = s.velocity - dv_n*w.normal;
        double dv_t_norm = dv_t.norm();
        double Ft_mag = std::min(s.friction*Fn_mag, s.damp_tangent*dv_t_norm);
        s.force += w.normal*Fn_mag - Ft_mag*dv_t/dv_t_norm;
    }
}
