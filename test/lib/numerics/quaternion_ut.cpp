#include "unit_tester.hpp"

#include "avro_types.h"

#include "numerics/quaternion.h"
#include "graphics/math.h"

using namespace avro;

UT_TEST_SUITE( quaternion_test_suite )

UT_TEST_CASE( test1 )
{
  real_t tol = 1e-7; // glm stores floats so use single-precision
  real_t axis[3];
  real_t theta = 42.*M_PI/180;

  // rotation about x-axis
  axis[0] = 1; axis[1] = 0; axis[2] = 0;
  Quaternion qx(theta,axis);

  graphics::mat3 rx = qx.rotation_matrix();

  UT_ASSERT_NEAR( rx(0,0) , 1 , tol );
  UT_ASSERT_NEAR( rx(0,1) , 0 , tol );
  UT_ASSERT_NEAR( rx(0,2) , 0 , tol );
  UT_ASSERT_NEAR( rx(1,0) , 0 , tol );
  UT_ASSERT_NEAR( rx(1,1) , cos(theta) , tol );
  UT_ASSERT_NEAR( rx(1,2) , -sin(theta) , tol );
  UT_ASSERT_NEAR( rx(2,0) , 0 , tol );
  UT_ASSERT_NEAR( rx(2,1) , sin(theta) , tol );
  UT_ASSERT_NEAR( rx(2,2) , cos(theta) , tol );

  // rotation about y-axis
  axis[0] = 0; axis[1] = 1; axis[2] = 0;
  Quaternion qy(theta,axis);

  graphics::mat3 ry = qy.rotation_matrix();

  UT_ASSERT_NEAR( ry(0,0) , cos(theta) , tol );
  UT_ASSERT_NEAR( ry(0,1) , 0 , tol );
  UT_ASSERT_NEAR( ry(0,2) , sin(theta) , tol );
  UT_ASSERT_NEAR( ry(1,0) , 0 , tol );
  UT_ASSERT_NEAR( ry(1,1) , 1 , tol );
  UT_ASSERT_NEAR( ry(1,2) , 0 , tol );
  UT_ASSERT_NEAR( ry(2,0) , -sin(theta) , tol );
  UT_ASSERT_NEAR( ry(2,1) , 0 , tol );
  UT_ASSERT_NEAR( ry(2,2) , cos(theta) , tol );

  // rotation about z-axis
  axis[0] = 0; axis[1] = 0; axis[2] = 1;
  Quaternion qz(theta,axis);

  graphics::mat3 rz = qz.rotation_matrix();

  UT_ASSERT_NEAR( rz(0,0) , cos(theta) , tol );
  UT_ASSERT_NEAR( rz(0,1) , -sin(theta) , tol );
  UT_ASSERT_NEAR( rz(0,2) , 0 , tol );
  UT_ASSERT_NEAR( rz(1,0) , sin(theta) , tol );
  UT_ASSERT_NEAR( rz(1,1) , cos(theta) , tol );
  UT_ASSERT_NEAR( rz(1,2) , 0 , tol );
  UT_ASSERT_NEAR( rz(2,0) , 0 , tol );
  UT_ASSERT_NEAR( rz(2,1) , 0 , tol );
  UT_ASSERT_NEAR( rz(2,2) , 1 , tol );

}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( quaternion_test_suite )
