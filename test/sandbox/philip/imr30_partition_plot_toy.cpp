//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "unit_tester.hpp"

#include "graphics/application.h"
#include "graphics/bsp.h"
#include "graphics/colormap.h"
#include "graphics/math.h"
#include "graphics/postscript.h"

#include "library/ckf.h"
#include "library/factory.h"
#include "library/plots.h"

#include "mesh/mesh.h"

#include "numerics/geometry.h"

#include "json/json.hpp"

#include <fstream>
#include <iomanip>

using namespace avro;

UT_TEST_SUITE(  sandbox_partition_plot_toy )

UT_TEST_CASE( test1 )
{
  typedef Simplex type;

  index_t nb_processors = 4;
  std::string base = "/Users/pcaplan/Codes/geocl/avro/build/";

  Colormap colormap;
  colormap.change_style("bwr");
  float lims[] = {0,float(nb_processors-1)};
  colormap.set_limits(lims);

  graphics::OpenGL_Application vis;
  graphics::Window& window = vis.window();

  window.load_view( base + "partition-view.json" );

  for (index_t adapt_iter = 0; adapt_iter <= 5; adapt_iter++) {

    for (index_t pass = 0; pass < 3; pass++) {

      std::vector< std::shared_ptr<graphics::BSPTriangles> > triangles( nb_processors );


      for (index_t k = 0; k < nb_processors; k++) {

        std::string mesh_name = base + "release_mpi/tmp/processor" + std::to_string(k) + "_pass" + std::to_string(pass) + "_" + std::to_string(adapt_iter) + ".mesh";

        // read the topology for this processor
        std::shared_ptr<TopologyBase> ptopology = nullptr;
        std::shared_ptr<Mesh> pmesh = library::get_mesh(mesh_name,ptopology);
        Topology<type>& topology = *static_cast<Topology<type>*>(ptopology.get());
        topology.orient();

        triangles[k] = std::make_shared<graphics::BSPTriangles>();

        graphics::mat4 model_matrix = graphics::glm::identity();
        triangles[k]->build( topology , model_matrix );
      }

      std::string eps_filename("adapt" + std::to_string(adapt_iter) + "-pass" + std::to_string(pass) + ".eps" );
      graphics::PostScriptWriter writer(eps_filename);
      writer.begin(window.width(),window.height());


      for (index_t k = 0; k < nb_processors; k++) {

        float c[3];
        colormap.map( float(k) , c );
        graphics::vec3 color = { c[0] , c[1] , c[2] };

        color.print();

        std::vector<graphics::BSPTriangle*> T( triangles[k]->nb() );
        for (index_t i = 0; i < triangles[k]->nb(); i++)
          T[i] = (*triangles[k])[i].get();

        writer.write( T , window.camera().view_matrix() , window.camera().projection_matrix() , color );
      }
      writer.end();
    }
  }
}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( sandbox_partition_plot_toy )
