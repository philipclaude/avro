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

using namespace avro;

std::pair< std::vector<real_t> , std::vector<index_t> >
retrieve_mesh( const Context& context ) {
  std::vector<real_t> coordinates;
  std::vector<index_t> connectivity;
  context.retrieve_mesh(coordinates,connectivity);
  return {coordinates,connectivity};
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

  m.def("retrieve_mesh",&retrieve_mesh, "");

  py::class_<Context>(m,"Context")
    .def(py::init<coord_t,coord_t,coord_t>())
    .def("define_mesh",&Context::define_mesh)
    .def("define_geometry",static_cast<void (Context::*)(const std::string&)>(&Context::define_geometry), "define geometry from string")
    .def("attach_geometry",&Context::attach_geometry)
    .def("adapt",&Context::adapt);
}
