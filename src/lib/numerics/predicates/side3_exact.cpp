// avro: Adaptive Voronoi Remesher
// Copyright 2017-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "numerics/expansion.h"
#include "numerics/predicates.h"

namespace GEO {

namespace PCK {

Sign avro_side3_nd_exact_pck(const double* p0,const double* p1,const double* p2,const double* p3,
 const double* q0,const double* q1,const double* q2 ,unsigned short dim
){

 expansion& l1 = expansion_sq_dist(p1,p0,dim);
 expansion& l2 = expansion_sq_dist(p2,p0,dim);
 expansion& l3 = expansion_sq_dist(p3,p0,dim);

 expansion& a10 = expansion_dot_at(p1,q0,p0,dim).scale_fast(2.0);
 expansion& a11 = expansion_dot_at(p1,q1,p0,dim).scale_fast(2.0);
 expansion& a12 = expansion_dot_at(p1,q2,p0,dim).scale_fast(2.0);

 expansion& a20 = expansion_dot_at(p2,q0,p0,dim).scale_fast(2.0);
 expansion& a21 = expansion_dot_at(p2,q1,p0,dim).scale_fast(2.0);
 expansion& a22 = expansion_dot_at(p2,q2,p0,dim).scale_fast(2.0);

 expansion& a30 = expansion_dot_at(p3,q0,p0,dim).scale_fast(2.0);
 expansion& a31 = expansion_dot_at(p3,q1,p0,dim).scale_fast(2.0);
 expansion& a32 = expansion_dot_at(p3,q2,p0,dim).scale_fast(2.0);

 expansion& b00 = expansion_diff( expansion_product(a11,a22) , expansion_product(a12,a21) );
 expansion& b01 = expansion_diff(a21,a22);
 expansion& b02 = expansion_diff(a12,a11);

 expansion& b10 = expansion_diff( expansion_product(a12,a20) , expansion_product(a10,a22) );
 expansion& b11 = expansion_diff( a22 , a20 );
 expansion& b12 = expansion_diff(a10,a12);

 expansion& b20 = expansion_diff( expansion_product(a10,a21) , expansion_product(a11,a20) );
 expansion& b21 = expansion_diff( a20 , a21 );
 expansion& b22 = expansion_diff( a11 , a10 );

 expansion& Delta = expansion_sum3( b00,b10,b20 );
 Sign Delta_sign = Delta.sign();
 geo_assert(Delta_sign != ZERO);

 expansion& DeltaLambda0 = expansion_sum3(b00,expansion_product(b01,l1),expansion_product(b02,l2));
 expansion& DeltaLambda1 = expansion_sum3(b10,expansion_product(b11,l1),expansion_product(b12,l2));
 expansion& DeltaLambda2 = expansion_sum3(b20,expansion_product(b21,l1),expansion_product(b22,l2));

 expansion& r0 = expansion_product(Delta,l3);
 expansion& r1 = expansion_product(a30,DeltaLambda0);
 expansion& r2 = expansion_product(a31,DeltaLambda1);
 expansion& r3 = expansion_product(a32,DeltaLambda2);
 expansion& r = expansion_diff( r0 , expansion_sum3( r1 , r2 , r3 ) );

 Sign r_sign = r.sign();

 if (r_sign==ZERO)
 {
   const double* p_sort[4];
   p_sort[0] = p0;
   p_sort[1] = p1;
   p_sort[2] = p2;
   p_sort[3] = p3;
   std::sort(p_sort,p_sort+4);

   for (index_t i=0;i<4;++i)
   {

    if (p_sort[i] == p0)
    {
      const expansion& z0 = expansion_product( a30 , expansion_sum(b01,b02) );
      const expansion& z1 = expansion_product( a31 , expansion_sum(b11,b12) );
      const expansion& z2 = expansion_product( a32 , expansion_sum(b21,b22) );
      const expansion& z = expansion_diff( Delta , expansion_sum3(z0,z1,z2) );
      Sign result = z.sign();
      if(result != ZERO) return Sign(Delta_sign*result);
    }
    else if (p_sort[i] == p1)
    {
      const expansion& z = expansion_sum3( expansion_product(a30,b01) , expansion_product(a31,b11) , expansion_product(a32,b21) );
      Sign result = z.sign();
      if (result != ZERO) return Sign(Delta_sign*result);
    }
    else if (p_sort[i] == p2)
    {
      const expansion& z = expansion_sum3( expansion_product(a30,b02) , expansion_product(a31,b12) , expansion_product(a32,b22) );
      Sign result = z.sign();
      if (result!=ZERO) return Sign(Delta_sign*result);
    }
    else if (p_sort[i] == p3)
    {
      return NEGATIVE;
    }
  }
  geo_assert_not_reached;
 }
 return Sign(Delta_sign*r_sign);
}

} // PCK

} // GEO
