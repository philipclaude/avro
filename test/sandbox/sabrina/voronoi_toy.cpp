
#include "unit_tester.hpp"

#include "graphics/application.h"

#include "library/ckf.h"

using namespace avro;

UT_TEST_SUITE( voronoi_toy )

UT_TEST_CASE( test1 )
{
  std::vector<index_t> dims(2,10);
  dims[1] = 5;

  CKF_Triangulation topology(dims);

  graphics::Viewer vis;

  vis.add( topology );

  vis.run();

}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( voronoi_toy )
