#include "gran2d.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>

namespace py = pybind11;
using namespace pybind11::literals;

void bind_gran2d(py::module &m) {
    py::class_<granular_media_2d>(m, "granular_media_2d")
        .def(py::init<double>(), "dt"_a)
        .def("add_wall", &granular_media_2d::add_wall, "point"_a, "normal"_a)
        .def("add_grains", &granular_media_2d::add_grains, "position"_a, "radii"_a, "mass"_a, "young_mod"_a, "friction"_a, "damp_normal"_a, "damp_tangent"_a)
        .def("add_static_grains", &granular_media_2d::add_static_grains, "position"_a, "radii"_a, "mass"_a, "young_mod"_a, "friction"_a, "damp_normal"_a, "damp_tangent"_a)
        .def("step", &granular_media_2d::step)
        .def_property_readonly("position", &granular_media_2d::get_position)
        .def_readwrite("gravity", &granular_media_2d::gravity)
        .def_readonly("Nparticles", &granular_media_2d::Nparticles)
        .def_readonly("time", &granular_media_2d::time)
        .def_readonly("dt", &granular_media_2d::dt);
}
