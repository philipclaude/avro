// avro: Adaptive Voronoi Remesher
// Copyright 2017-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "numerics/expansion.h"
#include "numerics/predicates.h"

namespace GEO
{
namespace PCK
{

#if 0
// this takes a really long time to compile :/

double
avro_orient4d(const double* p0,const double* p1,const double* p2,const double* p3,const double* p4)
{
/* const expansion& l01 = expansion_sq_dist(p1,p0,4);
 const expansion& l02 = expansion_sq_dist(p2,p0,4);
 const expansion& l03 = expansion_sq_dist(p3,p0,4);
 const expansion& l04 = expansion_sq_dist(p4,p0,4);
 const expansion& l12 = expansion_sq_dist(p2,p1,4);
 const expansion& l13 = expansion_sq_dist(p3,p1,4);
 const expansion& l14 = expansion_sq_dist(p4,p1,4);
 const expansion& l23 = expansion_sq_dist(p3,p2,4);
 const expansion& l24 = expansion_sq_dist(p4,p2,4);
 const expansion& l34 = expansion_sq_dist(p4,p3,4);

 const expansion& one = expansion_sum(1.,0.);
 const expansion& zero = expansion_sum(0.,0.);

 const expansion& d1 = expansion_det5x5(one,l01,l02,l03,l04,one,zero,l12,l13,l14,one,l12,zero,l23,l24,l03,l13,l23,zero,l34,l04,l14,l24,l34,zero);
 const expansion& d2 = expansion_det5x5(one,zero,l02,l03,l04,one,l01,l12,l13,l14,one,l02,zero,l23,l24,one,l03,l13,zero,l34,one,l04,l24,l34,zero);
 expansion& Delta = expansion_sum(d1,d2);
*/
 
 const expansion& a11 = expansion_diff( p1[0] , p0[0] );
 const expansion& a12 = expansion_diff( p2[0] , p0[0] );
 const expansion& a13 = expansion_diff( p3[0] , p0[0] );
 const expansion& a14 = expansion_diff( p4[0] , p0[0] );

 const expansion& a21 = expansion_diff( p1[1] , p0[1] );
 const expansion& a22 = expansion_diff( p2[1] , p0[1] );
 const expansion& a23 = expansion_diff( p3[1] , p0[1] );
 const expansion& a24 = expansion_diff( p4[1] , p0[1] );

 const expansion& a31 = expansion_diff( p1[2] , p0[2] );
 const expansion& a32 = expansion_diff( p2[2] , p0[2] );
 const expansion& a33 = expansion_diff( p3[2] , p0[2] );
 const expansion& a34 = expansion_diff( p4[2] , p0[2] );

 const expansion& a41 = expansion_diff( p1[3] , p0[3] );
 const expansion& a42 = expansion_diff( p2[3] , p0[3] );
 const expansion& a43 = expansion_diff( p3[3] , p0[3] );
 const expansion& a44 = expansion_diff( p4[3] , p0[3] );

 const expansion& Delta = expansion_det4x4( a11 , a12 , a13 , a14 , a21 , a22 , a23 , a24 , a31 ,a32,a33,a34 , a41,a42,a43,a44 ); 

 return Delta.value();

 //if (Delta_sign==ZERO) return 0.;
 return 1.;

}

#endif

} // PCK

} // GEO
