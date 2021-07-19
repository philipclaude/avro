#include "unit_tester.hpp"

#include "graphics/application.h"
#include "graphics/managers.h"
#include "graphics/vao.h"
#include "graphics/window.h"

#include "graphics/colormap.h"
#include "graphics/shader.h"

#include "library/ckf.h"
#include "library/egads.h"

#include "mesh/points.h"
#include "mesh/topology.h"

#include "voronoi/delaunay.h"
#include "voronoi/optimal_transport.h"

#if AVRO_HEADLESS_GRAPHICS == 0

using namespace avro;
using namespace avro::graphics;

UT_TEST_SUITE( graphics_decomposition_suite )

class TestField : public Field<Polytope,real_t> {

public:
  TestField( Topology<Polytope>& topology ) :
    Field( topology , 0, DISCONTINUOUS )
  {
    build();

    for (index_t k = 0; k < topology.nb(); k++) {
      this->value(k) = k;
    }
  }

  index_t nb_rank() const { return 1; }
};

UT_TEST_CASE( vao_polytopes_test )
{
  typedef Polytope type;
  coord_t number = 3;
  coord_t dim = number;
  index_t nb_points = 1e2;

  // create random delaunay vertices
  CubeDomain<type> domain(dim,dim,10);
  Delaunay delaunay( dim );
  std::vector<index_t> elems;
  for (index_t k = 0; k < nb_points; k++) {

    index_t elem = 0;
    std::vector<real_t> p(dim,0.);
    for (coord_t d = 0; d < dim; d++)
      p[d] = random_within(0.0,1.0);

    // create the delaunay point and retain which element in the domain it is in
    delaunay.create(p.data());
    elems.push_back( elem );
  }

  // initialize and compute the laguerre diagram
  delaunay::LaguerreDiagram<type> diagram( delaunay , domain );
  diagram.set_elements( elems );

  diagram.compute(false,nullptr);

  std::shared_ptr<TestField> fld = std::make_shared<TestField>(diagram);
  diagram.fields().make("test",fld);

  Viewer app;
  app.add(diagram);
  app.run(AVRO_FULL_UNIT_TEST);

}
UT_TEST_CASE_END( vao_polytopes_test )


UT_TEST_SUITE_END( graphics_decomposition_suite )

#else
UT_TEST_SUITE( graphics_decomposition_suite )
UT_TEST_SUITE_END( graphics_decomposition_suite )

#endif
