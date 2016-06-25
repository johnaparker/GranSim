#ifndef GUARD_vec_h
#define GUARD_vec_h

#include <iostream>
#include <cmath>
#include <eigen3/Eigen/Dense>
class Matrix3d;


class Vector3d{
public:
	double r[3];
public:
	Vector3d();
	Vector3d(Eigen::Vector3d);
	Vector3d(double,double,double);
	double operator[] (int i) const {return r[i];}
	Vector3d operator+ (const Vector3d & v) const {
		Vector3d v1(r[0]+v.r[0], r[1]+v.r[1], r[2]+v.r[2]);
		return v1;
	}
	Vector3d operator- (const Vector3d & v) const {
		Vector3d v1(r[0]-v.r[0], r[1]-v.r[1], r[2]-v.r[2]);
		return v1;
	}
	Vector3d operator- () const {
		Vector3d v1(-r[0], -r[1], -r[2]);
		return v1;
	}
	Vector3d operator* (const Vector3d & v) {
		Vector3d v1(r[0]*v[0],r[1]*v[1],r[2]*v[2]);
		return v1;
	}
	friend Vector3d operator* (const Vector3d & v, const double c) {
		Vector3d v1(c*v[0], c*v[1], c*v[2]);
		return v1;
	}
	friend Vector3d operator* (const double c, const Vector3d & v) {
		Vector3d v1(c*v[0], c*v[1], c*v[2]);
		return v1;
	}
	Vector3d operator/ (const Vector3d & v) const{
		Vector3d v1(r[0]/v[0],r[1]/v[1],r[2]/v[2]);
		return v1;
	}
	Vector3d operator/ (const double c) const {
		Vector3d v1(r[0]/c, r[1]/c, r[2]/c);
		return v1;
	}
	Vector3d & operator=(double a) {
		r[0] = a;
		r[1] = a;
		r[2] = a;
		return *this;
	}

	Vector3d & operator=(const Vector3d v2) {
		r[0] = v2[0];
		r[1] = v2[1];
		r[2] = v2[2];
		return *this;
	}

	Vector3d & operator=(const Eigen::Vector3d v2) {
		r[0] = v2(0);
		r[1] = v2(1);
		r[2] = v2(2);
		return *this;
	}

	Vector3d & operator+= (const Vector3d & v2) {
		r[0] += v2[0]; r[1] += v2[1]; r[2] += v2[2];
		return *this;
	}

	Vector3d & operator+= (const double c) {
		r[0] += c; r[1] += c; r[2] += c;
		return *this;
	}

	Vector3d & operator-= (const Vector3d & v2) {
		r[0] -= v2[0]; r[1] -= v2[1]; r[2] -= v2[2];
		return *this;
	}

	Vector3d & operator-= (const double c) {
		r[0] -= c; r[1] -= c; r[2] -= c;
		return *this;
	}

	Vector3d & operator*= (const Vector3d & v2) {
		r[0] *= v2[0]; r[1] *= v2[1]; r[2] *= v2[2];
		return *this;
	}

	Vector3d & operator*= (const double c) {
		r[0] *= c; r[1] *= c; r[2] *= c;
		return *this;
	}

	Vector3d & operator/= (const Vector3d & v2) {
		r[0] /= v2[0]; r[1] /= v2[1]; r[2] /= v2[2];
		return *this;
	}

	Vector3d & operator/= (const double c) {
		r[0] /= c; r[1] /= c; r[2] /= c;
		return *this;
	}

	double dot(Vector3d v2) const{
		double x = r[0]*v2[0] + r[1]*v2[1] + r[2]*v2[2];
		return x;
	}
	Vector3d cross(Vector3d v2) const {
		double x = r[1]*v2[2] - r[2]*v2[1];
		double y = r[2]*v2[0] - r[0]*v2[2];
		double z = r[0]*v2[1] - r[1]*v2[0];
		return Vector3d(x,y,z);

	}
	void print() {std::cout << "(" << r[0] << ", " << r[1] << ", " << r[2] << ")" << std::endl;}
	double norm() {
		return std::pow(r[0]*r[0] + r[1]*r[1] + r[2]*r[2],0.5);
	}
	void normalize() {
		double l = (*this).norm();
		(*this) /= l;
	}
	void set(double a, double b, double c){
		r[0] = a; r[1] = b; r[2] = c;
	}
	Matrix3d outer(Vector3d & v2) const;

};

class Matrix3d {
public:
	Vector3d rows[3];
public:
	Matrix3d();
	Matrix3d(double);
	Matrix3d(Eigen::Matrix3d);
	Matrix3d(Vector3d &, Vector3d &, Vector3d &);

	Matrix3d operator* (const Matrix3d & v) const {
		Vector3d v1(rows[0][0]*v[0][0] + rows[0][1]*v[1][0] + rows[0][2]*v[2][0],
						rows[0][0]*v[0][1] + rows[0][1]*v[1][1] + rows[0][2]*v[2][1],
						rows[0][0]*v[0][2] + rows[0][1]*v[1][2] + rows[0][2]*v[2][2]);
		Vector3d v2(rows[1][0]*v[0][0] + rows[1][1]*v[1][0] + rows[1][2]*v[2][0],
						rows[1][0]*v[0][1] + rows[1][1]*v[1][1] + rows[1][2]*v[2][1],
						rows[1][0]*v[0][2] + rows[1][1]*v[1][2] + rows[1][2]*v[2][2]);
		Vector3d v3(rows[2][0]*v[0][0] + rows[2][1]*v[1][0] + rows[2][2]*v[2][0],
						rows[2][0]*v[0][1] + rows[2][1]*v[1][1] + rows[2][2]*v[2][1],
						rows[2][0]*v[0][2] + rows[2][1]*v[1][2] + rows[2][2]*v[2][2]);
		return Matrix3d(v1,v2,v3);
	}
	Vector3d operator* (const Vector3d & v) const {
		Vector3d v1(rows[0].dot(v), rows[1].dot(v), rows[2].dot(v));
		return v1;
	}
	Matrix3d operator* (const double c) const {
		Vector3d v1 = c*rows[0];
		Vector3d v2 = c*rows[1];
		Vector3d v3 = c*rows[2];
		return Matrix3d(v1,v2,v3);
	}
	Matrix3d & operator+= (const Matrix3d & m2) {
		rows[0] += m2[0]; rows[1] += m2[1]; rows[2] += m2[2];
		return *this;
	}

	Matrix3d & operator-= (const Matrix3d & m2) {
		rows[0] -= m2[0]; rows[1] -= m2[1]; rows[2] -= m2[2];
		return *this;
	}

	Matrix3d operator- (const Matrix3d & m2) const {
		Vector3d v1(rows[0] - m2[0]);
		Vector3d v2(rows[1] - m2[1]);
		Vector3d v3(rows[2] - m2[2]);
		return Matrix3d(v1,v2,v3);
	}

	Matrix3d operator+ (const Matrix3d & m2) const {
		Vector3d v1(rows[0] + m2[0]);
		Vector3d v2(rows[1] + m2[1]);
		Vector3d v3(rows[2] + m2[2]);
		return Matrix3d(v1,v2,v3);
	}

	Matrix3d & operator=(const Matrix3d v2) {
		rows[0] = v2[0];
		rows[1] = v2[1];
		rows[2] = v2[2];
		return *this;
	}

	Vector3d operator[] (int i) const {return rows[i];}

	void compute(Vector3d &,Matrix3d &) const;
	Matrix3d inverse() const;

	void print() {
		std::cout << "[";
		rows[0].print();
		std::cout << " ";
		rows[1].print();
		std::cout << " ";
		rows[2].print();
		std::cout << "]" << std::endl;
	}


};


#endif