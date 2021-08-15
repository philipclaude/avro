//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "programs.h"

#include "common/directory.h"
#include "common/process.h"
#include "common/tools.h"

#include "graphics/application.h"

#include "library/factory.h"
#include "library/library.h"

#include "mesh/mesh.h"

#include "numerics/predicates.h"

#include "avro.h"

typedef avro::real_t REAL;
#include <tetgen1.5.0/predicates.h>

#include <triangle/predicates.h>

#include <stdio.h>

namespace avro
{

namespace programs
{

#if 0

void
random_point_in_simplex( const Points& points , const index_t* v , index_t nv , std::vector<real_t>& p )
{
  coord_t dim = points.dim();

  // calculate random barycentric coordinates
  std::vector<real_t> alpha(nv);
  alpha[nv-1] = 1.0;
  for (coord_t j=0;j<nv-1;j++)
  {
    alpha[j]    = real_t(rand())/real_t(RAND_MAX);
    alpha[nv-1] -= alpha[j];
  }

  // use barycentric coordinates to compute p
  std::fill(p.begin(),p.end(),0.);
  for (coord_t i=0;i<nv;i++)
  {
    for (coord_t d=0;d<dim;d++)
      p[d] += alpha[i] * points[v[i]][d];
  }
}

void
CentroidalVoronoiTessellation::sample_geometry( Entity* entity , index_t nb_samples )
{
  int idx = -1;
  for (index_t k=0;k<entities_.size();k++)
  {
    if (entities_[k] == entity)
    {
      idx = k;
      break;
    }
  }
  avro_assert( idx>=0 );
  const Topology<Simplex>& topology = *topologies_[idx].get();

  coord_t np = topology.number() +1;
  index_t simplex_begin = 0;
  index_t simplex_end = topology.nb();
  coord_t dim = sites_.dim();

  avro_assert( dim == topology.points().dim() );
  avro_assert( topology.nb()>0 );

  if (np==1)
  {
    avro_assert(topology.nb()==1);
    index_t id = sites_.nb();
    sites_.create(topology.points()[ topology(0)[0] ]);
    sites_.set_entity(id,entity);
    return;
  }

  // generate random numbers for every sample
  std::vector<real_t> s(nb_samples);
  for (index_t i=0;i<nb_samples;i++)
    s[i] = real_t( double(rand())/double(RAND_MAX) );
  std::sort(s.begin(),s.end());

  // measure the volume of each element and the total volume
  std::vector<real_t> vt(topology.nb());
  real_t vtot = 0.;
  topology.get_volumes(vt);
  for (index_t k=0;k<vt.size();k++)
    vtot += vt[k];

  int first_s = -1;
  int last_s  = 0;

  index_t cur_t = simplex_begin;
  real_t cur_s = vt[0]/vtot;
  for (index_t i=0;i<nb_samples;i++)
  {
    avro_assert( i<s.size() );
    while (s[i] > cur_s && cur_t < simplex_end-1 )
    {
      cur_t++;
      avro_assert( cur_t < simplex_end );
      cur_s += vt[cur_t]/vtot;
    }
    if (first_s==-1)
      first_s = int(cur_t);
    last_s = ( last_s > int(cur_t) ) ? last_s : int(cur_t);

    // generate random point in simplex
    std::vector<real_t> cur_p(dim);
    random_point_in_simplex( topology.points() , topology(cur_t) , topology.nv(cur_t) , cur_p );

    // do something with the point
    index_t id = sites_.nb();
    sites_.create(cur_p.data());
    sites_.set_entity(id,entity);
  }
}

#endif

int
voronoi( int nb_input , const char** inputs )
{

#if 0

  // so far only simplex adaptation is supported
  typedef Simplex type;

  if (nb_input<4 || nb_input==-1)
  {
    printf("\t\tvoronoi [input mesh] [sites] [geometry] [optional]\n");
    printf("\t\t--> optional can be:\n");
    printf("\t\t\thierarchical=(bool, whether the CVT is computed in a hierarchical manner)\n");
    printf("\t\t\tnb_iter=(int, # CVT iterations)\n");
    printf("\t\t\tnb_sites=(int, # sites for random only)\n");
    return 1;
  }

  ProcessCPU::initialize();
  ProcessGPU::initialize();

  // initialize the predicates
  exactinit(1,0,0,10,10,10);
  exactinit();
  GEO::PCK::initialize();

  // options
  bool found; UNUSED(found);
  const char **options = inputs +3;
  int  nb_options = nb_input -3;

  bool hierarchical = false;

  // retrieve the number of smoothing iterations
  index_t nb_iter = 0;
  if (nb_input>3)
    found = parse(lookfor(options,nb_options,"nb_iter"),nb_iter);

  // get the original input mesh
  std::string meshname( inputs[0] );
  std::shared_ptr<TopologyBase> ptopology = nullptr;
  std::shared_ptr<Mesh> pmesh = library::get_mesh(meshname,ptopology);

  // get the topology and add it to the input mesh
  Topology<type>& topology = *static_cast<Topology<type>*>(ptopology.get());
  topology.orient();
  pmesh->add(ptopology);

  // get the input geometry
  bool curved = false;
  std::string geometryname( inputs[2] );
  std::shared_ptr<Model> pmodel;
  if (geometryname!="none")
  {
    pmodel = library::get_geometry( geometryname , curved );

    // check the points are on the geometry...
    // option to project them
    pmesh->points().attach(*pmodel);
  }

  // get the input metric
  std::string sitesname( inputs[1] );
  Points sites( pmesh->points().dim() );

  if (sitesname=="points")
    pmesh->points().copy( sites );
  else if (sitesname=="random")
  {
    // retrieve the number of sites (if random)
    index_t nb_sites = 100;
    if (nb_input>3)
      found = parse(lookfor(options,nb_options,"nb_sites"),nb_sites);

    std::vector<real_t> xmin( sites.dim() ,  1e20 );
    std::vector<real_t> xmax( sites.dim() , -1e20 );
    for (index_t k=0;k<pmesh->points().nb();k++)
    {
      for (coord_t d=0;d<pmesh->points().dim();d++)
      {
        if (pmesh->points()[k][d]<xmin[d]) xmin[d] = pmesh->points()[k][d];
        if (pmesh->points()[k][d]>xmax[d]) xmax[d] = pmesh->points()[k][d];
      }
    }

    std::vector<real_t> x(sites.dim());
    for (index_t k=0;k<nb_sites;k++)
    {
      for (coord_t d=0;d<sites.dim();d++)
        x[d] = random_within( xmin[d] , xmax[d] );
      sites.create(x.data());
    }
  }
  else if (sitesname=="exact")
  {
    avro_implement;
  }

  std::shared_ptr<delaunay::CentroidalVoronoiTessellation> cvt;
  cvt = std::make_shared<delaunay::CentroidalVoronoiTessellation>( topology , sites , hierarchical );

  if (sitesname=="sample")
  {
    cvt->generate_sites(100);
  }

  cvt->compute(nb_iter);

  std::shared_ptr<Mesh> pmesh_out = std::make_shared<Mesh>(cvt->number(),cvt->points().dim());
  pmesh_out->add( cvt );
  cvt->points().copy( pmesh_out->points() );

  Library* lib = Library::get();
  lib->add_mesh_ptr(pmesh_out);
  char mesh_label[128];
  sprintf(mesh_label,"%p",(void*)pmesh_out.get());
  lib->add_mesh(mesh_label);

  printf("done computing CVT!\n");

#else
  printf("this functionality needs to be updated\n");
#endif

  return 0;
}

} // program

} // avro
