#include "gran.hpp"
#include <cmath>
#include <iostream>
#include <algorithm>

Circle::Circle(vec2 position, double radius, double mass, double young_mod, double friction, double damp_normal, double damp_tangent): position(position), radius(radius), mass(mass), young_mod(young_mod), friction(friction), damp_normal(damp_normal), damp_tangent(damp_tangent) {

    velocity = vec2(0,0);
    rd2 = vec2(0,0);
    rd3 = vec2(0,0);
    rd4 = vec2(0,0);
    force = vec2(0,0);
}

Wall2d::Wall2d(vec2 point, vec2 normal): point(point), normal(normal) {
    normal.normalize();
    tangent = vec2(normal(1), normal(0));
}

granular_media_2d::granular_media_2d(double dt): dt(dt) {
    time = 0;
    Nparticles = 0;
    Rparticles = 0;
    Vparticles = 0;

    gravity = vec2(0,-9.8);
}

void granular_media_2d::add_wall(vec2 point, vec2 normal) {
    walls.push_back(Wall2d(point, normal));
}

void granular_media_2d::add_grains(py_arr position, py_arr radii, py_arr mass, py_arr young_mod, py_arr friction, py_arr damp_normal, py_arr damp_tangent) {
    auto position_ = position.unchecked<2>();
    auto radii_ = radii.unchecked<1>();
    auto mass_ = mass.unchecked<1>();
    auto young_mod_ = young_mod.unchecked<1>();
    auto friction_ = friction.unchecked<1>();
    auto damp_normal_ = damp_normal.unchecked<1>();
    auto damp_tangent_ = damp_tangent.unchecked<1>();

    const int Nnew = position.shape(0);
    Rparticles += Nnew;
    Nparticles += Nnew;

    for (int i=0; i < Nnew; i++) {
        Circle circle(vec2(position_(i,0), position_(i,1)),
                       radii_(i), mass_(i), young_mod_(i), friction_(i),
                       damp_normal_(i), damp_tangent_(i));

        d_grains.push_back(circle);
    }

    initialize_voxels();
}

void granular_media_2d::predict() {
	const double a1 = dt;
    const double a2 = a1*dt/2.0;
    const double a3 = a2*dt/3.0;
    const double a4 = a3*dt/4.0;

    #pragma omp parallel for
    for (int i=0; i<Rparticles; i++) {
        auto& grain = d_grains[i];

        grain.position += a1*grain.velocity 
                          + a2*grain.rd2
                          + a3*grain.rd3
                          + a4*grain.rd4;

        grain.velocity += a1*grain.rd2
                     + a2*grain.rd3
                     + a4*grain.rd4;

        grain.rd2 += a1*grain.rd3
                + a2*grain.rd4;

        grain.rd3 += a1*grain.rd4;
    }
}

void granular_media_2d::correct() {
    const double c0 = 19.0/180.0*pow(dt,2);
    const double c1 = 3.0/8.0*dt;
    const double c3 = 3.0/2.0/dt;
    const double c4 = 1.0/pow(dt,2);

    #pragma omp parallel for
    for (int i=0; i<Rparticles; i++) {
        auto& grain = d_grains[i];

        auto accel = grain.force/grain.mass;
        auto corr = accel - grain.rd2;
        grain.position += c0*corr;
        grain.velocity += c1*corr;
        grain.rd2 = accel;
        grain.rd2 += c3*corr;
        grain.rd4 += c4*corr;
    }

    time += dt;
}

void granular_media_2d::compute_force() {
    #pragma omp parallel for
    for (int i=0; i<Rparticles; i++) {
        d_grains[i].force = vec2(0,0);
    }

    #pragma omp parallel for
    for (int i=0; i<Rparticles; i++) {
        auto& g1 = d_grains[i];
        int ix = int(g1.position(0)/voxel_size);
        int jx = int(g1.position(1)/voxel_size);

        for (int ix2=ix-1; ix2<ix+2; ix2++) {
            for (int jx2=jx-1; jx2<jx+2; jx2++) {
                key_tt key(ix2,jx2);
                auto loc = voxels.find(key);
                if (loc == voxels.end()) continue;

                for (int j: loc->second) {
                    if (i >= j) continue;

                    auto& g2 = d_grains[j];
                    vec2 dr = g1.position - g2.position;
                    bool condition = (dr.squaredNorm() < (g1.radius + g2.radius)*(g1.radius + g2.radius));

                    if (condition) {
                        interact(g1, g2);
                    }
                }
            }
        }
    }

    #pragma omp parallel for
    for (int i=0; i<Rparticles; i++) {
        auto& grain = d_grains[i];

        // gravity
        grain.force += grain.mass*gravity;

        // wall collisions
        for (const auto& wall: walls)
            interact(grain, wall);
    }
}

void granular_media_2d::initialize_voxels() {
    voxel_size = 2*d_grains[0].radius;

    for (int i=0; i<Rparticles; i++) {
        const auto& grain = d_grains[i];

        int ix = int(grain.position(0)/voxel_size);
        int jx = int(grain.position(1)/voxel_size);
        key_tt key(ix,jx);
        voxel_idx.push_back(key);

        auto loc = voxels.find(key);
        if (loc == voxels.end()) {
            voxels.insert({{key, {i}}});
        }
        else {
            loc->second.push_back(i);
        }
    }
}

void granular_media_2d::assign_voxels() {
    for (int i=0; i<Rparticles; i++) {
        auto& grain = d_grains[i];

        int ix = int(grain.position(0)/voxel_size);
        int jx = int(grain.position(1)/voxel_size);

        key_tt key(ix,jx);
        key_tt key_prev(voxel_idx[i]);
        if (key == key_prev) {
            continue;
        }
        else {
            auto loc = voxels.find(key);
            if (loc == voxels.end()) {
                voxels.insert({{key, {i}}});
            }
            else {
                loc->second.push_back(i);
            }
            
            auto loc_prev = voxels.find(key_prev);
            if (loc_prev != voxels.end()) {
                auto& vec = loc_prev->second;
                vec.erase(std::remove(vec.begin(), vec.end(), i), vec.end());
            }

            voxel_idx[i] = key;
        }
    }
}

void granular_media_2d::step() {
    predict();
    assign_voxels();
    compute_force();
    correct();
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
    double Fn_mag = std::max(0.0, sqrt(reff)*c1.young_mod*sqrt(overlap)*(overlap + c1.damp_normal*dv_n));
    double Ft_mag = std::min(c1.friction*Fn_mag, c1.damp_tangent*std::abs(dv_t));
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
