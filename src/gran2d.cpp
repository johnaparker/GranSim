#include "gran2d.hpp"
#include <cmath>
#include <iostream>
#include <algorithm>

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

void granular_media_2d::add_static_grains(py_arr position, py_arr radii, py_arr mass, py_arr young_mod, py_arr friction, py_arr damp_normal, py_arr damp_tangent) {
    auto position_ = position.unchecked<2>();
    auto radii_ = radii.unchecked<1>();
    auto mass_ = mass.unchecked<1>();
    auto young_mod_ = young_mod.unchecked<1>();
    auto friction_ = friction.unchecked<1>();
    auto damp_normal_ = damp_normal.unchecked<1>();
    auto damp_tangent_ = damp_tangent.unchecked<1>();

    const int Nnew = position.shape(0);
    Vparticles += Nnew;
    Nparticles += Nnew;

    for (int i=0; i < Nnew; i++) {
        Circle circle(vec2(position_(i,0), position_(i,1)),
                       radii_(i), mass_(i), young_mod_(i), friction_(i),
                       damp_normal_(i), damp_tangent_(i));

        s_grains.push_back(circle);
    }

    initialize_voxels();
}

py::array_t<double> granular_media_2d::get_position() {
    auto result = py::array_t<double>({Nparticles,2});
    auto r = result.mutable_unchecked<2>();

    #pragma omp parallel for
    for (int i=0; i<Nparticles; i++) {
        const auto& grain = (i < Rparticles) ? d_grains[i] : s_grains[i-Rparticles];
        r(i,0) = grain.position(0);
        r(i,1) = grain.position(1);
    }
    return result;
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

                    auto& g2 = (j < Rparticles) ? d_grains[j] : s_grains[j-Rparticles];
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
    voxels.clear();
    voxel_idx.clear();
    voxel_size = 2*d_grains[0].radius;

    for (int i=0; i<Nparticles; i++) {
        const auto& grain = (i < Rparticles) ? d_grains[i] : s_grains[i-Rparticles];

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
    for (int i=0; i<Nparticles; i++) {
        const auto& grain = (i < Rparticles) ? d_grains[i] : s_grains[i-Rparticles];

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
