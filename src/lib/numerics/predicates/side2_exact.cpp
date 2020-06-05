//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//

#include "numerics/expansion.h"
#include "numerics/predicates.h"

namespace GEO {

namespace PCK {

Sign avro_side2_nd_exact_pck(const double* p0, const double* p1, const double* p2,
 const double* q0, const double* q1 ,unsigned short dim
){

 expansion& l1 = expansion_sq_dist(p1,p0,dim);
 expansion& l2 = expansion_sq_dist(p2,p0,dim);
 expansion& a10 = expansion_dot_at(p1,q0,p0,dim).scale_fast(2.0);
 expansion& a11 = expansion_dot_at(p1,q1,p0,dim).scale_fast(2.0);

 expansion& a20 = expansion_dot_at(p2,q0,p0,dim).scale_fast(2.0);
 expansion& a21 = expansion_dot_at(p2,q1,p0,dim).scale_fast(2.0);

 expansion& Delta = expansion_diff(a11, a10);
 Sign Delta_sign = Delta.sign();
 geo_assert(Delta_sign != ZERO);

 expansion& DeltaLambda0 = expansion_diff(a11,l1);
 expansion& DeltaLambda1 = expansion_diff(l1,a10);

 expansion& r0 = expansion_product( Delta,l2);
 expansion& r1 = expansion_product( a20,DeltaLambda0).negate();
 expansion& r2 = expansion_product( a21,DeltaLambda1).negate();
 expansion& r = expansion_sum3(r0,r1,r2);

 Sign r_sign = r.sign();

 if (r_sign==ZERO)
 {
   const double* p_sort[3];
   p_sort[0] = p0;
   p_sort[1] = p1;
   p_sort[2] = p2;
   std::sort(p_sort,p_sort+3);
   for (index_t i=0;i<3;++i)
   {
     if (p_sort[i] == p0)
     {
       const expansion& z1 = expansion_diff(Delta,a21);
       const expansion& z = expansion_sum(z1,a20);
       Sign z_sign = z.sign();
       if(z_sign!=ZERO) return Sign(Delta_sign*z_sign);
     }
     if (p_sort[i] == p1)
     {
       const expansion& z = expansion_diff(a21,a20);
       Sign z_sign = z.sign();
       if(z_sign != ZERO) return Sign(Delta_sign*z_sign);
     }
     if (p_sort[i] == p2)
     {
       return NEGATIVE;
     }
   }
 }
 return Sign(Delta_sign*r_sign);
}

} // PCK

} // GEO
