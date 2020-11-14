#include "gran.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>

namespace py = pybind11;
using namespace pybind11::literals;

void bind_gran_sim(py::module &m) {
    py::class_<GranSim>(m, "GranSim")
        .def(py::init<const Matrix&, const Array&, 
            const Array&, double, double, double, double, double, const Matrix&, const Array&>(), 
            py::arg("position"), py::arg("radii"), py::arg("mass"),
            py::arg("young_mod"), py::arg("friction"), py::arg("damp_normal"),
            py::arg("damp_tangent"), py::arg("dt"), py::arg("vposition"), py::arg("vradii"))
        .def("step", &GranSim::step)
        .def("update_position", &GranSim::update_position)
        .def_readwrite("position", &GranSim::position)
        .def_readonly("radii", &GranSim::radii)
        .def_readonly("Nparticles", &GranSim::Nparticles)
        .def_readonly("time", &GranSim::time)
        .def_readonly("dt", &GranSim::dt);
}
