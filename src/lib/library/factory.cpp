// avro: Adaptive Voronoi Remesher
// Copyright 2017-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "adaptation/metric.h"

#include "common/directory.h"

#include "geometry/egads/context.h"
#include "geometry/egads/model.h"

#include "library/ckf.h"
#include "library/factory.h"
#include "library/field.h"
#include "library/meshb.h"
#include "library/metric.h"
#include "library/tesseract.h"

#include "mesh/mesh.h"
#include "mesh/points.h"



namespace avro
{

namespace library
{

std::shared_ptr<MetricAttachment>
get_metric( const std::string& name , Points& points , bool& is_analytic ,
           const std::vector<real_t>& params )
{
  // default to discrete metric field
  is_analytic = false;

  // get the file extension
  std::string ext = get_file_ext(name);

  // if there is a file extension, read the mesh
  if (ext=="sol" || ext=="solb")
  {
    // read the file
    std::shared_ptr<MetricAttachment> pfld = std::make_shared<MetricAttachment>(points);
    pfld->from_solb(name);
    return pfld;
  }

  is_analytic = true;
  if (name=="Uniform")
  {

    avro_assert(params.size()==2);
    coord_t dim = 0;
    real_t h = 1;
    if (params.size()>0)
    {
      dim = coord_t(params[0]);
      h   = params[1];
    }
    MetricField_Uniform analytic(dim,h);
    return std::make_shared<MetricAttachment>(analytic,points);
  }
  if (name=="Linear-3d")
  {
    MetricField_UGAWG_Linear analytic;
    return std::make_shared<MetricAttachment>(analytic,points);
  }
  if (name=="Polar1")
  {
    MetricField_UGAWG_Polar1 analytic;
    return std::make_shared<MetricAttachment>(analytic,points);
  }
  if (name=="Polar2")
  {
    MetricField_UGAWG_Polar2 analytic;
    return std::make_shared<MetricAttachment>(analytic,points);
  }
  if (name=="Linear-4d")
  {
    real_t hmin = MetricField_Tesseract_Linear::hmin_default;
    if (params.size()>0) hmin = params[0];
    MetricField_Tesseract_Linear analytic(hmin);
    return std::make_shared<MetricAttachment>(analytic,points);
  }
  if (name=="Wave-4d")
  {
    MetricField_Tesseract_Wave analytic;
    return std::make_shared<MetricAttachment>(analytic,points);
  }
  avro_assert_not_reached;
  return nullptr;
}

std::shared_ptr<Model>
getGeometry( const std::string& name , bool& curved )
{
  // get the file extension
  std::string ext = get_file_ext(name);

  // check if this is an egads file
  if (ext=="egads")
  {
    curved = true;
    std::shared_ptr<EGADS::Model> pmodel = std::make_shared<EGADS::Model>(name);
    return pmodel;
  }

  // lookup the geometry in the library
  std::shared_ptr<Model> pmodel = std::make_shared<Model>(0);
  std::shared_ptr<Body> pbody = nullptr;
/*
  if (name=="tesseract")
  {
    curved = false;
    std::vector<real_t> c(4,0.5);
    std::vector<real_t> L(4,1.0);
    pbody = std::make_shared<library::Tesseract>(c.data(),L.data());
  }
  if (name=="box")
  {
    curved = false;
    std::vector<real_t> x0(3,0.0);
    std::vector<real_t> L(3,1.0);
    pbody = std::make_shared<library::EGADSBox>(&context,x0.data(),L.data());
  }
  if (name=="square")
  {
    curved = false;
    std::vector<real_t> x0(3,0.5);
    x0[2] = 0.0;
    pbody = std::make_shared<library::EGADSSquare>(&context,x0.data(),1,1);
  }
*/

  avro_assert( pbody!=nullptr );
  pmodel->add_body(pbody);

  return pmodel;
}

template<typename type>
std::shared_ptr<Mesh>
get_mesh( const std::string& name , std::shared_ptr<Topology<type>>& ptopology , coord_t number )
{
  // get the file extension
  std::string ext = get_file_ext(name);

  // if there is a file extension, read the mesh
  if (ext=="mesh" || ext=="meshb")
  {
    std::shared_ptr<Mesh> pmesh = std::make_shared<meshb>(name);
    ptopology = pmesh->retrieve_ptr<type>(0);
    return pmesh;
  }
  if (ext=="json")
  {
    /*
    std::shared_ptr<Mesh> pmesh = std::make_shared<Mesh<type>>(number,number);
    io::readMesh(name,*pmesh);
    ptopology = pmesh->topology_ptr(0);
    return pmesh;
    */
    avro_implement;
  }

  // if no file extension, it must be a mesh in the library
  std::vector<std::string> s = split(name,"-");
  if (s[0]=="CKF")
  {
    // count how many hyphens there are to read the dimensions
    index_t nb_hyphen = s.size() -1;
    number = nb_hyphen;

    std::vector<real_t> lens(number,1.0);
    std::vector<index_t> dims(number,0);
    for (coord_t i=0;i<number;i++)
      dims[i] = unstringify<index_t>(s[i+1]);

    std::shared_ptr<Mesh> pmesh = std::make_shared<Mesh>(number,number);
    ptopology = std::make_shared<CKF_Triangulation>(dims);
    ptopology->orient(); // make positive volumes
    return pmesh;
  }

  avro_assert_not_reached;
  return nullptr;
}

template std::shared_ptr<Mesh> get_mesh( const std::string& , std::shared_ptr<Topology<Simplex>>& , coord_t );

} // programs

} // avro
