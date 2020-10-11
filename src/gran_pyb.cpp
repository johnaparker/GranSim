#include "gran.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>

namespace py = pybind11;
using namespace pybind11::literals;

void bind_gran_sim(py::module &m) {
    py::class_<GranSim>(m, "GranSim")
        .def(py::init<const Matrix&, const Array&, 
            const Array&, double, double, double, double, double>(), 
            py::arg("position"), py::arg("radii"), py::arg("mass"),
            py::arg("young_mod"), py::arg("friction"), py::arg("damp_normal"),
            py::arg("damp_tangent"), py::arg("dt"))
        .def("step", &GranSim::step)
        .def_readwrite("position", &GranSim::position);
}
