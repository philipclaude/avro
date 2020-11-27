#include "unit_tester.hpp"

#include "graphics/application.h"

#include "library/ckf.h"
#include "library/plots.h"

#include "mesh/facets.h"

#include "numerics/geometry.h"

#include "voronoi/delaunay.h"
#include "voronoi/optimal_transport.h"

#include "visualize.h"

#include <fstream>
#include <iomanip>

UT_TEST_SUITE( sandbox_semidiscrete_ot_toy )

UT_TEST_CASE( test1 )
{
  std::string prefix = "/Users/pcaplan/Dropbox/Codes/mach-II/avro/test/tmp/sdot-dim4-10000";

  library::meshb mesh( prefix+"_tet.mesh");
  std::ifstream field( prefix+"_sites.json");
  json J = json::parse(field);
  std::vector<index_t> sites = J["field"];

  Topology<Simplex>& tet = mesh.retrieve<Simplex>(0);
  std::shared_ptr<SliceSites> ts = std::make_shared<SliceSites>(tet,sites);
  tet.fields().make("sites",ts);

  graphics::Visualizer vis;
  vis.add_topology(tet);

  vis.run();
}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( sandbox_semidiscrete_ot_toy )
