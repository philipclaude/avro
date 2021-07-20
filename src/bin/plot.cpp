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

#include <fstream>
#include <stdio.h>

namespace avro
{

namespace programs
{

int
plot( int nb_input , const char** inputs , bool webplot )
{
  // so far only simplex adaptation is supported
  typedef Simplex type;

  if (nb_input<1 || nb_input==-1)
  {
    printf("\t\tplot/webplot [mesh] [optional]\n");
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

  lib->add_mesh_ptr(pmesh);
  char mesh_label[128];
  sprintf(mesh_label,"%p",(void*)pmesh.get());
  lib->add_mesh(mesh_label);

  // get the input geometry if provided
  std::shared_ptr<Model> pmodel = nullptr;
  std::string geometryname = lookfor(options,nb_options,"geometry");
  if (!geometryname.empty())
    pmodel = library::get_geometry( geometryname , curved );

  if (geometryname.empty()) {

    // check if this is an avro file, and if a geomeetyr was specified
    std::string ext = get_file_ext(meshname);
    if (ext == "avro") {

      std::ifstream file(meshname);
      nlohmann::json jm;
      file >> jm;

      try {
        geometryname = jm["geometry"];
        pmodel = library::get_geometry( geometryname,curved );
        printf("--> loaded geometry: %s\n",geometryname.c_str());
      }
      catch (...) {
        printf("could not load geometry from avro file");
      }
    }
  }


  // check the points are on the geometry...
  if (pmodel!=nullptr)
  {
    // check if the points are already on the geometry
    // TODO

    // if not...project them
    mesh.points().attach(*pmodel);
  }


  graphics::Viewer app(webplot);
  app.add(topology);
  app.run();

  return 0;
}

} // programs

} // avro
