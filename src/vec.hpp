#ifndef GUARD_vec_h
#define GUARD_vec_h

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>

using vec2 = Eigen::Vector2d;
using ivec2 = Eigen::Vector2i;
using quat = Eigen::Quaterniond;
using Array = Eigen::ArrayXd;
using iArray = Eigen::ArrayXi;

using Matrix = Eigen::Array<double,Eigen::Dynamic,2,Eigen::RowMajor>;

#endif
