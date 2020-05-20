#include "programs.h"

#include "common/directory.h"
#include "common/process.h"
#include "common/tools.h"

#include "graphics/application.h"

#include "library/factory.h"
#include "library/library.h"

#include "mesh/mesh.h"

#include "numerics/predicates.h"

#include "voronoi/cvt.h"

#include "avro.h"

typedef avro::real_t REAL;
#include <tetgen1.5.0/predicates.h>

#include <triangle/predicates.h>

#include <stdio.h>

namespace avro
{

namespace programs
{

int
voronoi( int nb_input , const char** inputs )
{

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

  return 0;
}

} // program

} // avro
