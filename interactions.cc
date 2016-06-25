#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "system.h"
#include "vec.h"
#include "molecule.h"
#include "interactions.h"
#include <eigen3/Eigen/Geometry>

using Eigen::Matrix3d;	

using namespace std;
using std::vector;
using std::pow;

void interact(sphere& sphere1, sphere& sphere2, bool tanForce, bool staticFriction, msystem* sys) {
	static double k_n = sys->params.k_n;
	static double k_t = sys->params.k_t;
	static double k_static = sys->params.k_static;
	static double cutoff = sys->params.cutoff;

	molecule* mol1 = sphere1.mol;
	molecule* mol2 = sphere2.mol;
	if (mol1 == mol2) return;

	Vector3d dr = sphere1.pos - sphere2.pos;
	double overlap = pow(dr.dot(dr),.5) - (sphere1.r + sphere2.r);
	if (overlap < 0) {
		dr.normalize();
		Vector3d interactPoint = sphere2.pos + sphere2.r*dr;
		Vector3d dv = mol1->pointVelocity(interactPoint) - mol2->pointVelocity(interactPoint);
		double dv_n = -dv.dot(dr);
		double Fn_mag = max(0.0,pow(-overlap,.5)*(-overlap*k_n + dv_n*k_t));
		Vector3d Fn = Fn_mag*dr;

		Vector3d F = Fn;

		if (tanForce) {
			Vector3d dv_t = dv + dv_n*dr;
			double dv_t_len = dv_t.norm();
			if (dv_t_len) dv_t /= dv_t_len;
			else dv_t = Vector3d();
			double Ft_mag;
			Vector3d elong;
			double elong_len;

			if (staticFriction) {
				contact c = mol1->num < mol2->num ? contact(sphere1, sphere2):contact(sphere2, sphere1);
				sys->contacts_mutex.lock();
				elong = (sys->contacts[c] += dv_t*dv_t_len*sys->timeStep);
				sys->contacts_mutex.unlock();
				elong_len = elong.norm();
				Ft_mag = min(k_static*elong_len,cutoff*Fn_mag);
			}
			else Ft_mag = cutoff*Fn_mag;

			Vector3d Ft = -Ft_mag*dv_t;

			F += Ft;
		}

		mol1->push(F,interactPoint);
		mol2->push(-F,interactPoint);
	}
}

void virtualInteract(sphere& sphere1, Vector3d loc, msystem* sys) {
	static double k_n = sys->params.k_n;
	static double k_t = sys->params.k_t;
	static double k_static = sys->params.k_static;

	molecule* mol1 = sphere1.mol;
	Vector3d dr = sphere1.pos - loc;
	double overlap = pow(dr.dot(dr),.5) - (2*sphere1.r);
	if (overlap < 0) {
		dr.normalize();
		Vector3d interactPoint = loc + sphere1.r*dr;
		Vector3d dv = mol1->pointVelocity(interactPoint);
		double dv_n = -dv.dot(dr);
		// Vector3d dv_t = dv + dv_n*dr;
		// double dv_t_len = dv_t.norm();
		// if (dv_t_len) dv_t /= dv_t_len;
		// else dv_t = Vector3d::Zero();
		double Fn_mag = max(0.0,pow(-overlap,.5)*(-overlap*k_n + dv_n*k_t));
		// double Ft_mag = 0.1*Fn_mag;
		Vector3d Fn = Fn_mag*dr;
		// Vector3d Ft = -Ft_mag*dv_t;
		// Vector3d F = Fn + Ft;
		Vector3d F = Fn;
		mol1->push(F,interactPoint);
	}
}

void wallInteract(sphere& sphere1, const vector<double> dim, msystem* mySys) {
	double x = sphere1[0], y = sphere1[1], z = sphere1[2], r = sphere1.r;
	static double omega = mySys->params.freq*M_PI/180.0;
	static double L = 5;
	double t = mySys->time;
	double theta = omega*t;
	double cost = cos(theta); double sint = sin(theta);

	if (y < dim[2] + r) {
		Vector3d pos(x,dim[2]-r,z);
		virtualInteract(sphere1,pos,mySys);
	}
	else if (y > dim[3] - r) {
		Vector3d pos(x,dim[3]+r,z);
		virtualInteract(sphere1,pos,mySys);
	}
	double D = L/2.0 + (x-L/2.0)*cost + z*sint;
	if (D < r) {
		Vector3d pos(x-(D+r)*cost,y,z - (D+r)*sint);
		virtualInteract(sphere1,pos,mySys);
	}
	D = L/2.0 - (x-L/2.0)*cost - z*sint;
	if (D < r) {
		Vector3d pos(x+(D+r)*cost,y,z + (D+r)*sint);
		virtualInteract(sphere1,pos,mySys);
	}
	D = (L/2.0-x)*sint + z*cost;
	if (D < r) {
		Vector3d pos(x+(D+r)*sint,y,z - (D+r)*cost);
		virtualInteract(sphere1,pos,mySys);
	}

	// if (x < dim[0] + r) {
	// 	Vector3d pos(dim[0]-r,y,z);
	// 	virtualInteract(sphere1,pos,mySys);
	// }
	// else if (x > dim[1] - r) {
	// 	Vector3d pos(dim[1]+r,y,z);
	// 	virtualInteract(sphere1,pos,mySys);
	// }
	// double D = L/2.0 + (y-L/2.0)*cost + z*sint;
	// if (D < r) {
	// 	Vector3d pos(x,y-(D+r)*cost,z - (D+r)*sint);
	// 	virtualInteract(sphere1,pos,mySys);
	// }
	// D = L/2.0 - (y-L/2.0)*cost - z*sint;
	// if (D < r) {
	// 	Vector3d pos(x,y+(D+r)*cost,z + (D+r)*sint);
	// 	virtualInteract(sphere1,pos,mySys);
	// }
	// D = (L/2.0-y)*sint + z*cost;
	// if (D < r) {
	// 	Vector3d pos(x,y+(D+r)*sint,z - (D+r)*cost);
	// 	virtualInteract(sphere1,pos,mySys);
	// }

	// if (x < dim[0] + r) {
	// 	Vector3d pos(dim[0]-r,y,z);
	// 	virtualInteract(sphere1,pos,mySys);
	// }
	// else if (x > dim[1] - r) {
	// 	Vector3d pos(dim[1]+r,y,z);
	// 	virtualInteract(sphere1,pos,mySys);
	// }
	// if (y < dim[2] + r) {
	// 	Vector3d pos(x,dim[2]-r,z);
	// 	virtualInteract(sphere1,pos,mySys);
	// }
	// else if (y > dim[3] - r) {
	// 	Vector3d pos(x,dim[3]+r,z);
	// 	virtualInteract(sphere1,pos,mySys);
	// }
	// if (z < dim[4] + r) {
	// 	Vector3d pos(x,y,dim[4]-r);
	// 	virtualInteract(sphere1,pos,mySys);
	// }
	// else if (z > dim[5] - r) {
	// 	Vector3d pos(x,y,dim[5]+r);
	// 	virtualInteract(sphere1,pos,mySys);
	// }
}