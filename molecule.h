#ifndef GUARD_molecule_h
#define GUARD_molecule_h

#include <eigen3/Eigen/Dense>
#include <vector>
#include "vec.h"


Matrix3d RotMatrix(Vector3d &, double);
class molecule;

class sphere {
public:
	Vector3d pos;
	double r;
	double m;
	int cell_num;
	molecule* mol;
public:
	sphere(Vector3d pos,double r,double m):pos(pos),r(r),m(m) {};
	double operator[] (int i) {return pos[i];};
};


class moleculeType {
public:
	std::vector<sphere> spheres;
	std::vector<sphere>::size_type size;
	double M, vol;
	Vector3d orientation, evals, COM;
	Matrix3d I, evecs;
public:
	moleculeType(std::vector<sphere>);
	sphere operator[] (int i) {return spheres[i];};
};

class molecule{
public:
	moleculeType *type;          
	std::vector<sphere> spheres;
	Matrix3d evecs;
	Vector3d pos, orientation, vel, angVel, rcm;
	Vector3d rd0, rd1, rd2, rd3, rd4;
	Vector3d od0, od1, od2, od3, od4;
	Vector3d F, tau, *evals;
	double* M;
	int num;
public:
	molecule(moleculeType &, Vector3d pos = Vector3d(0,0,0), Vector3d vel = Vector3d(0,0,0),
				Vector3d orientation = Vector3d(0,0,1), Vector3d angVel = Vector3d(0,0,0));
	void displace(Vector3d);
	void setPosition(Vector3d);
	void rotate(Vector3d axis);
	void rotate(Vector3d axis, double theta);
	void rotate(Vector3d axis, double theta, Vector3d rotPoint);
	void clearForces();
	void push(Vector3d F, Vector3d r);
	void predict(double dt);
	void correct(double dt);
	void updateCOM();
	Vector3d pointVelocity(Vector3d);
	sphere operator[] (int i) {return spheres[i];};
};

#endif