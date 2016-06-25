#include <iostream>
#include <stdexcept>
#include <cmath>
#include <vector>
#include <numeric>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/Eigenvalues>
#include "vec.h"
#include "molecule.h"

using namespace std;						using std::vector;
using std::pow;				
using Eigen::SelfAdjointEigenSolver;

Matrix3d RotMatrix(Vector3d & axis, double theta) {
	if (!axis.dot(axis)) throw invalid_argument("Cannot rotate about the null vector");

	axis.normalize();
	double& t = theta;
	double ux = axis[0], uy = axis[1], uz = axis[2];

	Matrix3d R;
	Vector3d v1(cos(t)+ux*ux*(1-cos(t)),  ux*uy*(1-cos(t))-uz*sin(t),  ux*uz*(1-cos(t))+uy*sin(t));
	Vector3d v2(uy*ux*(1-cos(t))+uz*sin(t),  cos(t)+uy*uy*(1-cos(t)),  uy*uz*(1-cos(t))-ux*sin(t));
	Vector3d v3(uz*ux*(1-cos(t))-uy*sin(t),  uz*uy*(1-cos(t))+ux*sin(t),  cos(t)+uz*uz*(1-cos(t)));
	Matrix3d(v1,v2,v3);

	return Matrix3d(v1,v2,v3);
}

moleculeType::moleculeType(vector<sphere> spheres):spheres(spheres) {
	size = spheres.size();
	M = 0.0;
	for (vector<sphere>::const_iterator sp = spheres.begin(); sp != spheres.end(); ++sp) {
		M += sp->m;
		COM += sp->pos*sp->m;
		vol += 4.0/3.0*M_PI*pow(sp->r,3);
	}
	COM /= M;
	orientation.set(0,0,1);

	I = Matrix3d();
	for (vector<sphere>::const_iterator sp = spheres.begin(); sp != spheres.end(); ++sp) {
		Vector3d R = sp->pos - COM;
		I += (Matrix3d(2.0/5.0*sp->m*pow(sp->r,2)) + (Matrix3d(R.dot(R)) - R.outer(R))*sp->m);
	}

	I.compute(evals,evecs);
}

molecule::molecule(moleculeType & molType, Vector3d pos, Vector3d vel, Vector3d orientation,
						Vector3d angVel): pos(pos), vel(vel), orientation(orientation), angVel(angVel){
	spheres = molType.spheres;
	type = &molType;
	evecs = molType.evecs;
	evals = &molType.evals;
	M = &(molType.M);

	Matrix3d A = evecs.inverse();
	rd0 = pos; rd1 = vel; rd2 = Vector3d(); rd3 = Vector3d(); rd4 = Vector3d();
	od0 = orientation; od1 = A*angVel; od2 = Vector3d(); od3 = Vector3d(); od4 = Vector3d();
	F = Vector3d(); tau = Vector3d(); updateCOM(); 

	Vector3d o_pos = pos; Vector3d o_orient = orientation;
	(*this).pos = Vector3d(); (*this).orientation = Vector3d(0,0,1);
	displace(o_pos);
	rotate(o_orient);
	(*this).pos = o_pos;
	(*this).orientation = o_orient;

	for (vector<sphere>::iterator sp = (*this).spheres.begin(); sp != (*this).spheres.end(); ++sp) {
		sp-> mol = this;
	}
}

void molecule::displace(Vector3d dr) {
	for (vector<sphere>::iterator sp = spheres.begin(); sp != spheres.end(); ++sp) {
			sp->pos += dr;
	}
	pos += dr;
	updateCOM();
}
void molecule::setPosition(Vector3d r) {
	Vector3d dr = r - pos;
	displace(dr);
}

void molecule::rotate(Vector3d axis) {
	Vector3d rotAxis = orientation.cross(axis);
	double phi = acos(orientation.dot(axis)/pow((axis.dot(axis)*orientation.dot(orientation)),0.5));
	rotate(rotAxis,phi);
}
void molecule::rotate(Vector3d axis, double theta) {
	rotate(axis,theta,rcm);
}
void molecule::rotate(Vector3d axis, double theta, Vector3d r) {
	Matrix3d R;
	try{R = RotMatrix(axis,theta);}
	catch(invalid_argument){return;}
	for (vector<sphere>::iterator sp = spheres.begin(); sp != spheres.end(); ++sp) {
			sp->pos = R*(sp->pos - r) + r;		
	}
	evecs = R*evecs;
	orientation = R*orientation;
	updateCOM();
}

void molecule::clearForces() {
	F = Vector3d(); tau = Vector3d();
}

void molecule::push(Vector3d Force, Vector3d r) {
	F += Force;
	tau += (r-rcm).cross(Force);
}

void molecule::updateCOM() {
	rcm = Vector3d();
	for (vector<sphere>::iterator sp = spheres.begin(); sp != spheres.end(); ++sp) {
			rcm += sp->m*sp->pos;	
	}
	rcm /= *M;
}

Vector3d molecule::pointVelocity(Vector3d r) {
	Vector3d dr = r - rcm;
	Vector3d om = evecs*od1;
	return om.cross(dr) + rd1;
}

void molecule::predict(double dt) {
	static double a1 = dt, a2 = a1*dt/2.0, a3 = a2*dt/3.0, a4 = a3*dt/4.0;

	rd0 += (a1*rd1 + a2*rd2 + a3*rd3 + a4*rd4);
	rd1 += (a1*rd2 + a2*rd3 + a3*rd4);
	rd2 += (a1*rd3 + a2*rd4);
	rd3 += (a1*rd4);

	od1 += (a1*od2 + a2*od3 + a3*od4);
	od2 += (a1*od3 + a2*od4);
	od3 += (a1*od4);

	angVel = evecs*od1;
	setPosition(rd0);
	rotate(angVel,dt*pow(angVel.dot(angVel),.5));
}

void molecule::correct(double dt) {

	static double c0 = 19.0/180.0*dt*dt, c1 = 3.0/8.0*dt, c2 = 1.0,
					c3 = 3.0/(2.0*dt), c4 = 1.0/(dt*dt);
	Vector3d accel = F/ *M;
	Vector3d corr = accel - rd2;
	rd0 += c0*corr;   rd1 += c1*corr;
	rd2 = accel;	  rd3 += c3*corr;	rd4 += c4*corr;

	Matrix3d A = evecs.inverse();
	Vector3d temp = od1.cross((*evals)*od1);
	Vector3d angAccel = (A*tau - temp)/(*evals);
	corr = angAccel - od2;

	od1 += c1*corr;	  od2 = angAccel;
	od3 += c3*corr;	  od4 += c4*corr;

	angVel = evecs*od1;
	setPosition(rd0);
	rotate(angVel,dt*pow(angVel.dot(angVel),.5));

	// Vector3d accel = F/ *M;
	// Matrix3d A = evecs.inverse();
	// Vector3d temp = od1.cross((*evals)*od1);
	// Vector3d angAccel = (A*tau - temp)/(*evals);

	// rd1 += accel*dt;
	// rd0 += rd1*dt;
	// od1 += angAccel*dt;

	// angVel = evecs*od1;
	// setPosition(rd0);
	// rotate(angVel,dt*pow(angVel.dot(angVel),.5));
}


