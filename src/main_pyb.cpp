#define NOMINMAX
#include <pybind11/pybind11.h>

namespace py = pybind11;

void bind_gran2d(py::module &);
void bind_gran3d(py::module &);

PYBIND11_MODULE(cpp, m) {
    m.doc() = R"pbdoc(
        C++ submodule of foo
        -----------------------

        .. currentmodule:: cpp

        .. autosummary::
           :toctree: _generate
    )pbdoc";

    bind_gran2d(m);
    bind_gran3d(m);
}
