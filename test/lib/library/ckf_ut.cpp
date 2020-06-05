//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "unit_tester.hpp"

#include "library/ckf.h"

using namespace avro;

UT_TEST_SUITE( ckf_suite )

class CKF_Triangulation_test : public CKF_Triangulation
{
public:
  CKF_Triangulation_test( const std::vector<index_t>& dims ) :
    CKF_Triangulation(dims)
  {
    nb_simplices_ = 1;
    nb_points_    = 1;
    for (coord_t d=0;d<dims.size();d++)
    {
      nb_points_    *= dims[d];
      nb_simplices_ *= dims[d]-1;
    }
    nb_simplices_ *= numerics::factorial( dims.size() );
  }

  index_t nb_points() const { return nb_points_; }
  index_t nb_simplices() const { return nb_simplices_; }

private:
  index_t nb_points_;
  index_t nb_simplices_;
};

UT_TEST_CASE( test_constant_dims )
{

  for (coord_t dim=2;dim<=4;dim++)
  {
    for (index_t n=2;n<=8;n+=2)
    {
      std::vector<index_t> dims(dim,n);

      CKF_Triangulation_test topology( dims );

      UT_ASSERT_EQUALS( topology.nb() , topology.nb_simplices() );
      UT_ASSERT_EQUALS( topology.points().nb() , topology.nb_points() );
      UT_ASSERT_NEAR( topology.volume() , 1.0 , 1e-10 );
    }
  }

}
UT_TEST_CASE_END( test_constant_dims )

UT_TEST_CASE( test_different_dims )
{

  for (coord_t dim=2;dim<=4;dim++)
  {
    for (index_t n=2;n<=8;n+=2)
    {
      std::vector<index_t> dims(dim,n);

      dims[0] += 1;

      CKF_Triangulation_test topology( dims );

      UT_ASSERT_EQUALS( topology.nb() , topology.nb_simplices() );
      UT_ASSERT_EQUALS( topology.points().nb() , topology.nb_points() );
      UT_ASSERT_NEAR( topology.volume() , 1.0 , 1e-10 );
    }
  }

}
UT_TEST_CASE_END( test_different_dims )

UT_TEST_SUITE_END( ckf_suite )
