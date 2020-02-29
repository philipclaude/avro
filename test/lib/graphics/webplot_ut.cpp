#include "unit_tester.hpp"

#include "graphics/application.h"

#include "library/ckf.h"
#include "library/obj.h"
#include "library/samples.h"

#include "mesh/topology.h"

#include "voronoi/delaunay.h"
#include "voronoi/voronoi.h"

using namespace avro;
using namespace avro::graphics;

UT_TEST_SUITE( web_application_suite )

template<typename type>
class CoordinateField : public Field<type,std::vector<real_t>>
{
public:
  CoordinateField( const Topology<type>& topology ) :
    Field<type,std::vector<real_t>>(topology,topology.order(),CONTINUOUS)
  {
    this->build();
    for (index_t k=0;k<topology.points().nb();k++)
    {
      std::vector<real_t> X(topology.points()[k],topology.points()[k]+topology.points().dim());
      this->value(k) = X;
    }
  }


  index_t nb_rank() const override { return this->topology().points().dim(); }

  std::string get_name( index_t j ) const override
  {
    if (j==0) return "x";
    if (j==1) return "y";
    if (j==2) return "z";
    if (j==3) return "t";
    return std::to_string(j);
  }

};

UT_TEST_CASE( test1 )
{

  library::objFile topology0( "/Users/pcaplan/Google Drive/library/models/obj/spot.obj" );
  CKF_Triangulation topology1( {4,4} );
  CKF_Triangulation topology2( {4,4,4} );

  std::shared_ptr<CoordinateField<Simplex>> fld = std::make_shared<CoordinateField<Simplex>>(topology0);
  topology0.fields().make( "coordinates" , fld );

  WebVisualizer vis;

  vis.add_topology(topology0);
  vis.add_topology(topology1);
  vis.add_topology(topology2);

  vis.run();
}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( web_application_suite )
