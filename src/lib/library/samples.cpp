#include "library/samples.h"

namespace ursa
{

namespace library
{

TwoTriangles::TwoTriangles() :
  Topology<Simplex<Lagrange>>(vertices_,2),
  edges_(vertices_,1),
  vertices_(3)
{
  real_t x0[3] = {0,0,1};
  real_t x1[3] = {1,0,2};
  real_t x2[3] = {0,1,3};
  real_t x3[3] = {1,1,4};
  vertices_.create( x0 );
  vertices_.create( x1 );
  vertices_.create( x2 );
  vertices_.create( x3 );

  index_t t0[3] = {0,1,2};
  index_t t1[3] = {1,3,2};
  index_t t2[3] = {0,3,2};

  add( t0 , 3 );
  add( t1 , 3 );
  add( t2 , 3 );

  index_t e0[2] = {0,1};
  index_t e1[2] = {1,2};
  index_t e2[2] = {2,0};
  index_t e3[2] = {1,3};
  index_t e4[2] = {3,2};

  edges_.add( e0 , 2 );
  edges_.add( e1 , 2 );
  edges_.add( e2 , 2 );
  edges_.add( e3 , 2 );
  edges_.add( e4 , 2 );

}

} // library

} // ursa
