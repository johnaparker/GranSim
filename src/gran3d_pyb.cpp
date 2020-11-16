#include "gran3d.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>

namespace py = pybind11;
using namespace pybind11::literals;

void bind_gran3d(py::module &m) {
    py::class_<granular_media>(m, "granular_media")
        .def(py::init<double>(), "dt"_a)
        .def("add_wall", &granular_media::add_wall, "point"_a, "normal"_a)
        .def("add_grains", &granular_media::add_grains, "position"_a, "radii"_a, "mass"_a, "young_mod"_a, "friction"_a, "damp_normal"_a, "damp_tangent"_a)
        .def("add_static_grains", &granular_media::add_static_grains, "position"_a, "radii"_a, "mass"_a, "young_mod"_a, "friction"_a, "damp_normal"_a, "damp_tangent"_a)
        .def("step", &granular_media::step)
        .def_property_readonly("position", &granular_media::get_position)
        .def_readwrite("gravity", &granular_media::gravity)
        .def_readonly("Nparticles", &granular_media::Nparticles)
        .def_readonly("time", &granular_media::time)
        .def_readonly("dt", &granular_media::dt);
}
