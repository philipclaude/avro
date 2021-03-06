//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include <pybind11/pybind11.h>

#include <pybind11/stl.h>

#include <avro.h>
#include <avro_params.h>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

using namespace avro;

std::pair< std::vector<real_t> , std::vector<index_t> >
retrieve_mesh( const Context& context ) {
  std::vector<real_t> coordinates;
  std::vector<index_t> connectivity;
  context.retrieve_mesh(coordinates,connectivity);
  return {coordinates,connectivity};
}

std::pair< std::vector<int> , std::vector<real_t> >
retrieve_geometry( const Context& context ) {
  std::vector<real_t> param;
  std::vector<int> geometry;
  context.retrieve_geometry(geometry,param);
  return {geometry,param};
}

std::pair< std::vector<std::vector<index_t>> , std::vector<int> >
retrieve_boundary( const Context& context ) {
  std::vector<std::vector<index_t>> faces;
  std::vector<int> geometry;
  context.retrieve_boundary_parallel(faces,geometry); // not necessarily parallel, but is most general
  return {faces,geometry};
}

std::tuple< std::vector<real_t> , std::vector<index_t> , std::vector<index_t> >
retrieve_polytopes( const Context& context ) {
  std::vector<real_t> vertices;
  std::vector<index_t> elements;
  std::vector<index_t> nv_per_elem;
  context.retrieve_polytopes(vertices,elements,nv_per_elem);
  return std::make_tuple(vertices,elements,nv_per_elem);
}

namespace py = pybind11;

PYBIND11_MODULE(pyavro, m) {
    m.doc() = R"pbdoc(
        avro python interface
        ---------------------

        .. currentmodule:: pyavro

        .. autosummary::
           :toctree: _generate

    )pbdoc";

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif

  m.def("retrieve_mesh",&retrieve_mesh, "");
  m.def("retrieve_boundary",&retrieve_boundary,"");
  m.def("retrieve_geometry",&retrieve_geometry,"");
  m.def("retrieve_polytopes",&retrieve_polytopes,"");

  py::class_<Context>(m,"Context")
    .def(py::init<coord_t,coord_t,coord_t>())
    .def("define_mesh",&Context::define_mesh)
    .def("define_geometry",static_cast<void (Context::*)(const std::string&)>(&Context::define_geometry), "define geometry from string")
    .def("attach_geometry",&Context::attach_geometry)
    .def("adapt",&Context::adapt)
    .def("compute_laguerre",&Context::compute_laguerre)
    .def("compute_optimal_transport",&Context::compute_optimal_transport)
    .def("plot",&Context::plot);

}
