#ifndef GUARD_vec_h
#define GUARD_vec_h

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

using vec2 = Eigen::Vector2d;
using ivec2 = Eigen::Vector2i;

using vec3 = Eigen::Vector3d;
using ivec3 = Eigen::Vector3i;

using quat = Eigen::Quaterniond;
using Array = Eigen::ArrayXd;
using iArray = Eigen::ArrayXi;

using Matrix = Eigen::Array<double,Eigen::Dynamic,2,Eigen::RowMajor>;

using py_arr = py::array_t<double>;

#endif
