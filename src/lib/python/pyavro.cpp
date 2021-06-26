#include <pybind11/pybind11.h>

#include <pybind11/stl.h>

#include <avro.h>
#include <avro_params.h>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

int add(int i, int j) {
    return i + j;
}

void hello_pydg() {
  printf("hello pydg!\n");
}

namespace py = pybind11;

PYBIND11_MODULE(pyavro, m) {
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: pydg

        .. autosummary::
           :toctree: _generate

           add
           subtract
    )pbdoc";

    m.def("add", &add, R"pbdoc(
        Add two numbers

        Some other explanation about the add function.
    )pbdoc");

    m.def("subtract", [](int i, int j) { return i - j; }, R"pbdoc(
        Subtract two numbers

        Some other explanation about the subtract function.
    )pbdoc");

    m.def("hello", &hello_pydg, R"pbdoc(
        Says hello pydg from c++!
    )pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif

  /*
  py::class_<ElementResidual>(m,"ElementResidual")
    .def(py::init<int>())
    .def_readonly("data",&ElementResidual::data_);*/
}
