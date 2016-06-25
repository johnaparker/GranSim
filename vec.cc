#include <iostream>
#include "vec.h"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>

using Eigen::SelfAdjointEigenSolver;
using namespace std;

Vector3d::Vector3d() {
	r[0] = 0; r[1] = 0; r[2] = 0;
}

Vector3d::Vector3d(Eigen::Vector3d v2) {
	r[0] = v2(0); r[1] = v2(1); r[2] = v2(2);
}

Vector3d::Vector3d(double a, double b, double c) {
	r[0] = a; r[1] = b; r[2] = c;
}

Matrix3d Vector3d::outer(Vector3d & v2) const {
	Vector3d a1(r[0]*v2[0],r[0]*v2[1],r[0]*v2[2]);
	Vector3d a2(r[1]*v2[0],r[1]*v2[1],r[1]*v2[2]);
	Vector3d a3(r[2]*v2[0],r[2]*v2[1],r[2]*v2[2]);
	return Matrix3d(a1,a2,a3);
}

Matrix3d::Matrix3d() {
	rows[0] = Vector3d();
	rows[1] = Vector3d();
	rows[2] = Vector3d();
}

Matrix3d::Matrix3d(double c) {
	rows[0] = Vector3d(c,0,0);
	rows[1] = Vector3d(0,c,0);
	rows[2] = Vector3d(0,0,c);
}

Matrix3d::Matrix3d(Eigen::Matrix3d M) {
	rows[0] = Vector3d(M(0,0),M(0,1),M(0,2));
	rows[1] = Vector3d(M(1,0),M(1,1),M(1,2));
	rows[2] = Vector3d(M(2,0),M(2,1),M(2,2));
}

Matrix3d::Matrix3d(Vector3d & a, Vector3d & b, Vector3d & c) {
	rows[0] = a;
	rows[1] = b;
	rows[2] = c;
}

void Matrix3d::compute(Vector3d & evals, Matrix3d & evecs) const {

	Eigen::Matrix3d I;
	I << rows[0][0], rows[0][1], rows[0][2],
		 rows[1][0], rows[1][1], rows[1][2],
		 rows[2][0], rows[2][1], rows[2][2];

	SelfAdjointEigenSolver<Eigen::Matrix3d> es;
	es.compute(I);
	Eigen::Vector3d v1 = es.eigenvalues();
	evals = Vector3d(v1);
	Eigen::Matrix3d m1 = es.eigenvectors();
	evecs = Matrix3d(m1);
}

Matrix3d Matrix3d::inverse() const {
	Eigen::Matrix3d A;
	A << rows[0][0], rows[0][1], rows[0][2],
	     rows[1][0], rows[1][1], rows[1][2],
	     rows[2][0], rows[2][1], rows[2][2];
	Matrix3d Ainv(A.inverse());
	return Ainv;
}
