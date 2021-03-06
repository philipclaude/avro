//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "adaptation/metric.h"

#include "common/directory.h"

#include "geometry/egads/context.h"
#include "geometry/egads/model.h"

#include "library/ckf.h"
#include "library/egads.h"
#include "library/factory.h"
#include "library/field.h"
#include "library/meshb.h"
#include "library/metric.h"
#include "library/tesseract.h"

#include "mesh/mesh.h"
#include "mesh/points.h"

#include <fstream>

namespace avro
{

namespace library
{

std::shared_ptr<MetricAttachment>
get_metric( const std::string& name , Points& points , bool& is_analytic ,
           const std::vector<real_t>& params ) {

  // default to discrete metric field
  is_analytic = false;

  // get the file extension
  std::string ext = get_file_ext(name);

  // if there is a file extension, read the mesh
  if (ext == "sol" || ext == "solb") {
    // read the file
    std::shared_ptr<MetricAttachment> pfld = std::make_shared<MetricAttachment>(points);
    pfld->from_solb(name);
    return pfld;
  }

  is_analytic = true;
  if (name == "Uniform") {

    avro_assert(params.size() == 2);
    coord_t dim = 0;
    real_t h = 1;
    if (params.size() > 0) {
      dim = coord_t(params[0]);
      h   = params[1];
    }
    MetricField_Uniform analytic(dim,h);
    return std::make_shared<MetricAttachment>(analytic,points);
  }
  if (name == "Linear-2d") {
    MetricField_UGAWG_Linear2d analytic;
    return std::make_shared<MetricAttachment>(analytic,points);
  }
  if (name == "Linear-3d") {
    MetricField_UGAWG_Linear analytic;
    return std::make_shared<MetricAttachment>(analytic,points);
  }
  if (name == "Polar1") {
    MetricField_UGAWG_Polar1 analytic;
    return std::make_shared<MetricAttachment>(analytic,points);
  }
  if (name == "Polar2") {
    MetricField_UGAWG_Polar2 analytic;
    return std::make_shared<MetricAttachment>(analytic,points);
  }
  if (name == "Linear-4d") {
    real_t hmin = MetricField_Tesseract_Linear::hmin_default;
    if (params.size() > 0) hmin = params[0];
    MetricField_Tesseract_Linear analytic(hmin);
    return std::make_shared<MetricAttachment>(analytic,points);
  }
  if (name == "RotatingBL-3d") {
    MetricField_Cube_RotatingBoundaryLayer analytic;
    return std::make_shared<MetricAttachment>(analytic,points);
  }
  if (name == "RotatingBL-4d") {
    MetricField_Tesseract_RotatingBoundaryLayer analytic;
    return std::make_shared<MetricAttachment>(analytic,points);
  }
  if (name == "MovingCylinder-4d") {
    MetricField_Tesseract_MovingCylinder analytic;
    return std::make_shared<MetricAttachment>(analytic,points);
  }
  if (name == "Wave-3d") {
    MetricField_Cube_Wave analytic;
    return std::make_shared<MetricAttachment>(analytic,points);
  }
  if (name == "Wave-4d") {
    MetricField_Tesseract_Wave analytic;
    return std::make_shared<MetricAttachment>(analytic,points);
  }
  avro_assert_not_reached;
  return nullptr;
}

std::shared_ptr<Model>
get_geometry( const std::string& name , bool& curved ) {

  // get the file extension
  std::string ext = get_file_ext(name);

  // check if this is an egads file
  if (ext == "egads" || ext == "legads") {
    curved = true;
    std::shared_ptr<EGADS::Model> pmodel = std::make_shared<EGADS::Model>(name);
    return pmodel;
  }

  // lookup the geometry in the library
  std::shared_ptr<Model> pmodel = std::make_shared<Model>(0);
  std::shared_ptr<Body> pbody = nullptr;

  if (name == "tesseract" || name == "tesseract-closingwall" || name == "tesseract-expansion") {
    curved = false;
    std::vector<real_t> c(4,0.5);
    std::vector<real_t> lengths(4,1.0);
    std::shared_ptr<Model> pmodel = std::make_shared<Model>(3);
    pbody = std::make_shared<library::Tesseract>(c,lengths);

    library::Tesseract& body = static_cast<library::Tesseract&>(*pbody);
    if (name == "tesseract-closingwall") {
      body.map_to( &library::Tesseract::closingwall);
    }
    else if (name == "tesseract-expansion") {
      body.map_to( &library::Tesseract::expansion );
    }

  }
  if (name == "box") {
    curved = false;
    real_t x0[3] = {0,0,0};
    std::shared_ptr<EGADS::Model> emodel = std::make_shared<EGADS::Model>(2);
    const EGADS::Context* context = &emodel->context();
    std::vector<real_t> lengths(3,1);
    pbody  = std::make_shared<EGADS::Cube>(context,lengths,x0);
    emodel->add_body(pbody);
    emodel->build();
    pmodel = emodel;
  }
  if (name == "square") {
    curved = false;
    //real_t x0[3] = {0.5,0.5,0.};
    std::shared_ptr<EGADS::Model> emodel = std::make_shared<EGADS::Model>(1);
    std::vector<real_t> lengths(2,1);
    const EGADS::Context* context = &emodel->context();
    pbody  = std::make_shared<EGADS::Cube>(context,lengths);
    emodel->add_body(pbody);
    emodel->build();
    pmodel = emodel;
  }
  avro_assert( pbody != nullptr );
  pmodel->add_body(pbody);

  return pmodel;
}

std::shared_ptr<Mesh>
get_mesh( const std::string& name , std::shared_ptr<TopologyBase>& ptopology , coord_t number ) {

  // get the file extension
  std::string ext = get_file_ext(name);

  // if there is a file extension, read the mesh
  if (ext == "mesh" || ext == "meshb") {
    std::shared_ptr<Mesh> pmesh = std::make_shared<meshb>(name);
    ptopology = pmesh->retrieve_ptr<Simplex>(0);
    printf("mesh has %lu topologies\n",pmesh->nb_topologies());
    if (pmesh->nb_topologies() > 1) {
      // if there is more than one topology in the meshb file, add them as children to the first one
      for (index_t k=1;k<pmesh->nb_topologies();k++)
        static_cast<Topology<Simplex>*>(ptopology.get())->add_child( pmesh->retrieve_ptr<Simplex>(k) );
    }
    return pmesh;
  }
  if (ext == "avro" || ext == "json") {

    std::ifstream file(name);
    nlohmann::json jm;
    file >> jm;

    coord_t dim = jm["dim"];
    number = jm["number"];
    std::shared_ptr<Mesh> pmesh  = std::make_shared<Mesh>(number,dim);

    index_t nb_ghost = 1;
    try {
      nb_ghost = jm["nb_ghost"];
    }
    catch (...) {
      printf("[warning] could not read number of ghost points, setting to 1\n");
    }

    std::vector<real_t> vertices = jm["vertices"];
    index_t nb_vertices = vertices.size() / dim;
    for (index_t k = 0; k < nb_vertices; k++) {
      pmesh->points().create( vertices.data() + k*dim );
    }

    if (jm["type"] == "simplex") {
      ptopology = std::make_shared<Topology<Simplex>>(pmesh->points(),number);

      std::vector<index_t> simplices = jm["elements"];
      index_t nb_simplices = simplices.size()/(number+1);
      printf("nb_simplices = %lu\n",nb_simplices);
      printf("nb_vertices = %lu\n",pmesh->points().nb());
      for (index_t k = 0; k < nb_simplices; k++) {
        bool ghost = false;
        std::vector<index_t> simplex( simplices.data() + k*(number+1) , simplices.data() + (k+1)*(number+1) );
        for (index_t j = 0; j < simplex.size(); j++)
          if (simplex[j] < nb_ghost) ghost = true;
        if (ghost) continue;
        ptopology->add( simplex.data(),simplex.size());
      }
    }
    else {
      std::string s = jm["type"];
      printf("unsupported element type %s\n",s.c_str());
      avro_implement;
    }

    pmesh->add(ptopology);
    return pmesh;

  }

  if (ext.empty()) {
    // check if this is a memory address
    if (name.substr(0,2) == "0x") {
      unsigned long address = std::stoul(name,0,16);
      Mesh* mesh0 = (Mesh*) address;
      std::vector<const TopologyBase*> topologies;
      mesh0->retrieve(topologies);
      avro_assert( topologies.size() > 0 );
      index_t maxnumber = 0;
      const TopologyBase* topology = nullptr;
      for (index_t k = 0; k < topologies.size(); k++) {
        //printf("topology(%lu) number = %u with %lu elements\n",k,topologies[k]->number(),topologies[k]->nb());
        if (topologies[k]->number()<=maxnumber) continue;
        maxnumber = topologies[k]->number();
        topology  = topologies[k];
      }
      printf("determined mesh number = %lu, using topology with %lu elements\n",maxnumber,topology->nb());

      std::shared_ptr<Mesh> pmesh = std::make_shared<Mesh>(maxnumber,mesh0->points().dim());
      mesh0->points().copy( pmesh->points() );

      printf("topology has %lu fields to copy\n",topology->fields().nb());
      printf("adding mesh with %lu points\n",pmesh->points().nb());
      ptopology = nullptr;
      for (index_t k = 0;k < mesh0->nb_topologies(); k++) {
        if (mesh0->topology_ptr(k).get() == topology) {
          ptopology = mesh0->topology_ptr(k);
          break;
        }
      }
      avro_assert( ptopology!=nullptr );

      printf("topology type = %s\n",ptopology->type_name().c_str());
      pmesh->add( ptopology );
      return pmesh;
    }
  }

  // if no file extension, it must be a mesh in the library
  std::vector<std::string> s = split(name,"-");
  if (s[0] == "CKF" or s[0] == "CKF-simplex") {
    // count how many hyphens there are to read the dimensions
    index_t offset = 1;
    if (s[1] == "simplex") offset = 2;
    number = s.size() - offset;

    std::vector<real_t> lens(number,1.0);
    std::vector<index_t> dims(number,0);
    for (coord_t i=0;i<number;i++)
      dims[i] = unstringify<index_t>(s[i+offset]);

    std::shared_ptr<Mesh> pmesh = std::make_shared<Mesh>(number,number);
    std::shared_ptr<Topology<Simplex>> ptopology_ckf = std::make_shared<CKF_Triangulation>(dims);
    ptopology_ckf->orient(); // make positive volumes
    ptopology_ckf->points().copy( pmesh->points() );
    ptopology = ptopology_ckf;
    pmesh->add(ptopology_ckf);
    return pmesh;
  }
  if (s[0] == "CKF-cube") {
    avro_implement;
  }
  printf("cannot find mesh %s\n",name.c_str());
  avro_assert_not_reached;
  return nullptr;
}

} // programs

} // avro
