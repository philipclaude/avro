#include "programs.h"

#include "common/directory.h"
#include "common/process.h"
#include "common/tools.h"

#include "geometry/entity.h"
#include "geometry/model.h"

#include "graphics/application.h"
#include "graphics/bsp.h"
#include "graphics/postscript.h"

#include "library/factory.h"
#include "library/library.h"
#include "library/tesseract.h"

#include "mesh/boundary.h"
#include "mesh/mesh.h"

#include <fstream>
#include <stdio.h>

namespace avro
{

namespace programs
{

int
image( int nb_input , const char** inputs )
{
  // so far only simplex adaptation is supported
  typedef Simplex type;

  if (nb_input < 3 || nb_input == -1)
  {
    printf("\t\timage [mesh] [view] [outfile] width=600 height=600\n");
    return 1;
  }

  ProcessCPU::initialize();
  ProcessGPU::initialize();

  // options (these are actually required)
  const char **options = inputs +2;
  int  nb_options = nb_input -2;

  std::string meshname( inputs[0] );
  std::string viewname( inputs[1] );
  std::string outfile( inputs[2] );

  index_t width = 600, height = 600;
  parse<index_t>(lookfor(options,nb_options,"width"),width);
  parse<index_t>(lookfor(options,nb_options,"height"),height);

  // get the input mesh
  std::shared_ptr<TopologyBase> ptopology = nullptr;
  std::shared_ptr<Mesh> pmesh = library::get_mesh(meshname,ptopology);
  Mesh& mesh = *pmesh;
  Topology<type>& topology = *static_cast<Topology<type>*>(ptopology.get());

  // get the input geometry if provided
  bool curved;
  std::shared_ptr<Model> pmodel = nullptr;
  std::string geometryname = lookfor(options,nb_options,"geometry");
  if (!geometryname.empty())
    pmodel = library::get_geometry( geometryname , curved );

  if (geometryname.empty()) {

    // check if this is an avro file, and if a geometry was specified
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

  graphics::OpenGL_Application app;
  app.add(topology);

  graphics::Window& window = app.window();

  window.compute_view();
  window.load_view(viewname);
  window.resize(width,height);

  graphics::BSPTriangles bsp_triangles;
  graphics::mat4 ms;
  int i = 0;
  bsp_triangles.build( window.plot(i) , window.camera().view_matrix() , window.camera().projection_matrix() , ms );

  graphics::BSPTree tree;
  tree.build(bsp_triangles);
  std::vector<graphics::BSPTriangle*> triangles;
  const graphics::vec3& eye = window.camera().eye();
  tree.get_triangles( eye , triangles );

  std::string eps_filename(outfile);
  graphics::PostScriptWriter writer(eps_filename);
  writer.begin( window.width() , window.height() );
  writer.write( triangles , window.camera().view_matrix() , window.camera().projection_matrix() );
  writer.end();

  return 0;
}

} // programs

} // avro
