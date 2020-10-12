#include "gran.hpp"
#include <cmath>
#include <iostream>
#include <algorithm>

GranSim::GranSim(const Matrix& position, const Array& radii, 
        const Array& mass, double young_mod, double friction, double damp_normal,
        double damp_tangent, double dt): position(position),
        radii(radii), mass(mass), young_mod(young_mod), friction(friction),
        damp_normal(damp_normal), damp_tangent(damp_tangent), dt(dt) {

    time = 0;
    Nparticles = position.rows();

    velocity = Matrix::Zero(Nparticles,2);
    rd2 = Matrix::Zero(Nparticles,2);
    rd3 = Matrix::Zero(Nparticles,2);
    rd4 = Matrix::Zero(Nparticles,2);
    force = Matrix::Zero(Nparticles,2);

    voxel_size = 2*radii(0);

    for (int i=0; i<Nparticles; i++) {
        int ix = int(position(i,0)/voxel_size);
        int jx = int(position(i,1)/voxel_size);
        key_tt key(ix,jx);
        voxel_idx.push_back(key);
    }
}

void GranSim::predict() {
	const double a1 = dt;
    const double a2 = a1*dt/2.0;
    const double a3 = a2*dt/3.0;
    const double a4 = a3*dt/4.0;

    #pragma omp parallel for
    for (int i=0; i<Nparticles; i++) {
        position.row(i) += a1*velocity.row(i) 
                         + a2*rd2.row(i)
                         + a3*rd3.row(i)
                         + a4*rd4.row(i);

        velocity.row(i) += a1*rd2.row(i)
                         + a2*rd3.row(i)
                         + a4*rd4.row(i);

        rd2.row(i) += a1*rd3.row(i)
                    + a2*rd4.row(i);

        rd3.row(i) += a1*rd4.row(i);
    }
}

void GranSim::correct() {
    const double c0 = 19.0/180.0*pow(dt,2);
    const double c1 = 3.0/8.0*dt;
    const double c3 = 3.0/2.0/dt;
    const double c4 = 1.0/pow(dt,2);

    #pragma omp parallel for
    for (int i=0; i<Nparticles; i++) {
        auto accel = force.row(i)/mass(i);
        auto corr = accel - rd2.row(i);
        position.row(i) += c0*corr;
        velocity.row(i) += c1*corr;
        rd2.row(i) = accel;
        rd2.row(i) += c3*corr;
        rd4.row(i) += c4*corr;
    }

    time += dt;
}

void GranSim::compute_force() {
    force = 0;

    #pragma omp parallel for
    for (int i=0; i<Nparticles; i++) {
        int ix = int(position(i,0)/voxel_size);
        int jx = int(position(i,1)/voxel_size);

        for (int ix2=ix-1; ix2<ix+2; ix2++) {
            for (int jx2=jx-1; jx2<jx+2; jx2++) {
                key_tt key(ix2,jx2);
                auto loc = voxels.find(key);
                if (loc == voxels.end()) continue;

                for (int j: loc->second) {
                    if (i <= j) continue;

                    vec2 dr = position.row(i) - position.row(j);
                    bool condition = (dr.squaredNorm() < (radii(i) + radii(j))*(radii(i) + radii(j)));


                    if (condition) {
                        double dr_norm = dr.norm();
                        double overlap = radii(i) + radii(j) - dr_norm;
                        dr /= dr_norm;
                        vec2 dt(-dr(1), dr(0));
                        vec2 dv = velocity.row(i) - velocity.row(j);

                        double dv_n = -dr.dot(dv);
                        double dv_t = dt.dot(dv);
                        double reff = radii(i)*radii(j)/(radii(i) + radii(j));
                        double Fn_mag = std::max(0.0, sqrt(reff)*young_mod*sqrt(overlap)*(overlap + damp_normal*dv_n));
                        double Ft_mag = std::min(friction*Fn_mag, damp_tangent*std::abs(dv_t));
                        auto F = dr.array()*Fn_mag - dt.array()*sgn(dv_t)*Ft_mag;
                        force.row(i) += F;
                        force.row(j) -= F;
                    }
                }
            }
        }
    }

    #pragma omp parallel for
    for (int i=0; i<Nparticles; i++) {
        force.row(i) += mass(i)*vec2(0, -9.8).array();

        double overlap = radii(i) - position(i,1);
        if (overlap > 0) {
            double dv_n = -velocity(i,1);
            double dv_t = velocity(i,0);
            double reff = radii(i);
            double Fn_mag = std::max(0.0, sqrt(reff)*young_mod*sqrt(overlap)*(overlap + damp_normal*dv_n));
            double Ft_mag = std::min(friction*Fn_mag, damp_tangent*std::abs(dv_t));
            force.row(i) += vec2(-sgn(dv_t)*Ft_mag, Fn_mag).array();
        }
    }
}

void GranSim::assign_voxels() {
    for (int i=0; i<Nparticles; i++) {
        int ix = int(position(i,0)/voxel_size);
        int jx = int(position(i,1)/voxel_size);

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

void GranSim::step() {
    predict();
    assign_voxels();
    compute_force();
    correct();
}
