#if 0
#include "programs.h"

#include "common/directory.h"
#include "common/process.h"
#include "common/stringify.h"

#include "geometry/context.h"
#include "geometry/entity.h"

#include "graphics/plotter.h"

#include "library/metric.h"
#include "library/plots.h"

#include "mesh/boundary.h"
#include "mesh/gamma.h"
#include "mesh/mesh.h"

#include "numerics/metric.h"

#include <stdio.h>

namespace avro
{

namespace programs
{

int
plot( int nb_input , char** inputs )
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

  // options
  bool found;
  char **options = inputs +1;
  int  nb_options = nb_input -1;

  bool curved = true;
  bool all_bnd = false; // only number-1 (dim-1) children may be plotted
  found = parse<bool>(lookfor(options,nb_options,"all_bnd"),all_bnd);
  UNUSED(found);

  // get the input mesh
  std::string meshname( inputs[0] );
  std::shared_ptr<Topology<type>> ptopology = nullptr;
  std::shared_ptr<Mesh<type>> pmesh = getMesh<type>(meshname,ptopology);
  Mesh<type>& mesh = *pmesh;
  Topology<type>& topology = *ptopology.get();
  coord_t number = mesh.number();

  // get the input geometry if provided
  Context context;
  std::shared_ptr<Model> pmodel = nullptr;
  std::string geometryname = lookfor(options,nb_options,"geometry");
  if (!geometryname.empty())
    pmodel = getGeometry( geometryname , context , curved );

  // check the vertices are on the geometry...
  if (pmodel!=nullptr)
  {
    // check if the vertices are already on the geometry

    // if not...project them
    mesh.vertices().findGeometry(*pmodel);
  }

  // create a graphics plotter instance
  Server server;
  std::shared_ptr<Plotter> plotter = std::make_shared<Plotter>(&server);

  // save it globally (some internal functions use the plotter)
  __plotter__ = (void*) plotter.get();

  // option to plot the boundary
  Boundary<Simplex> boundary(topology);
  std::shared_ptr<library::BoundaryPlot<type>> bplot = nullptr;
  if (pmodel!=nullptr && number<4)
  {
    if (!all_bnd)
      boundary.extract(); // only num-1 children
    else
      boundary.extractall(); // all children
    bplot = std::make_shared<library::BoundaryPlot<type>>(boundary);
  }

  if (number==4)
    plotter->addSlice(topology);
  else
    plotter->addPlot( mesh );
  if (bplot!=nullptr)
    plotter->addPlot(bplot->mesh());

  plotter->run();

  return 0;
}



} // programs

} // avro

#endif
