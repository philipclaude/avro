#include "unit_tester.hpp"

#include "geometry/entity.h"
#include "geometry/tessellation.h"

#include "geometry/egads/context.h"
#include "geometry/egads/model.h"

#include "graphics/application.h"
#include "graphics/user_interface.h"
#include "graphics/window.h"

using namespace avro;

UT_TEST_SUITE( geometry_tessellation_suite )

UT_TEST_CASE(test1)
{
  EGADS::Context context;
  EGADS::Model model(&context,"data/cube-cylinder.egads" );
  //EGADS::Model model(&context,"/Users/pcaplan/Codes/EngSketchPad/data/basic/import_2.egads" );

  TessellationParameters params;
  params.standard();

  params.min_size() = 0.5;
  params.min_angle() = 20;

  ModelTessellation tess(model,params);

  tess.points().print(true);

  for (index_t k=0;k<tess.points().nb();k++)
  {
    if (tess.points().entity(k)->number()<=1)
      UT_ASSERT( tess.points().u(k,1) > 1e10 );
    else
      UT_ASSERT( tess.points().u(k,1) < 1e10 );
  }

  graphics::Visualizer vis;

  std::shared_ptr<graphics::Widget> toolbar = std::make_shared<graphics::Toolbar>(vis.main_window(),vis);
  vis.main_window().interface().add_widget( toolbar );

  vis.add_topology( tess.topology(0) );

  //vis.run();
}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END( geometry_tessellation_suite )
