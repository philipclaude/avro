#include "programs.h"

#include "common/directory.h"
#include "common/process.h"
#include "common/tools.h"

#include "geometry/entity.h"

#include "graphics/application.h"

#include "library/factory.h"
#include "library/library.h"

#include "mesh/boundary.h"
#include "mesh/mesh.h"

#include <stdio.h>

namespace avro
{

namespace programs
{

int
plot( int nb_input , const char** inputs )
{
  // so far only simplex adaptation is supported
  typedef Simplex type;

  if (nb_input<1 || nb_input==-1)
  {
    printf("\t\tplot [mesh] [optional]\n");
    printf("\t\t--> optional can be:\n");
    printf("\t\t\tgeometry=(string, EGADS filename or library geometry)\n");
    // TODO printf("\t\t\tfield=(string, .sol/solb file)\n");
    printf("\t\t\tall_bnd=(bool,whether all dimensional geometry children are plotted)\n\t\t\t\tdefault=false...only dim-1\n");
    return 1;
  }

  ProcessCPU::initialize();
  ProcessGPU::initialize();

  // the global library that we will add to
  Library* lib = Library::get();

  // options
  bool found;
  const char **options = inputs +1;
  int  nb_options = nb_input -1;

  bool curved = true;
  bool all_bnd = false; // only number-1 (dim-1) children may be plotted
  found = parse<bool>(lookfor(options,nb_options,"all_bnd"),all_bnd);
  UNUSED(found);

  // get the input mesh
  std::string meshname( inputs[0] );
  std::shared_ptr<TopologyBase> ptopology = nullptr;
  std::shared_ptr<Mesh> pmesh = library::get_mesh(meshname,ptopology);
  Mesh& mesh = *pmesh;
  Topology<type>& topology = *static_cast<Topology<type>*>(ptopology.get());
  coord_t number = mesh.number();

  lib->add_mesh_ptr(pmesh);
  char mesh_label[128];
  sprintf(mesh_label,"%p",(void*)pmesh.get());
  lib->add_mesh(mesh_label);

  // get the input geometry if provided
  std::shared_ptr<Model> pmodel = nullptr;
  std::string geometryname = lookfor(options,nb_options,"geometry");
  if (!geometryname.empty())
    pmodel = library::get_geometry( geometryname , curved );

  // check the points are on the geometry...
  if (pmodel!=nullptr)
  {
    // check if the points are already on the geometry

    // if not...project them
    mesh.points().attach(*pmodel);
  }

  graphics::Visualizer vis;
  std::shared_ptr<graphics::Widget> toolbar = std::make_shared<graphics::Toolbar>(vis.main_window(),vis);
  vis.main_window().interface().add_widget( toolbar );

  // option to plot the boundary
  Boundary<type> boundary(topology);
  if (pmodel!=nullptr && number<4)
  {
    if (!all_bnd)
      boundary.extract(); // only num-1 children
    else
      boundary.extractall(); // all children
    //bplot = std::make_shared<library::BoundaryPlot<type>>(boundary);
  }

  vis.add_topology( topology );
  vis.run();

  return 0;
}

} // programs

} // avro
