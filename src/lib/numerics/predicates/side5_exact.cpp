
#include "numerics/expansion.h"
#include "numerics/predicates.h"

namespace GEO
{
namespace PCK
{

Sign avro_side5_nd_exact_pck(const double* p0,const double* p1,const double* p2,const double* p3,const double* p4,const double* p5,
 const double* q0,const double* q1,const double* q2,const double* q3,const double* q4 ,unsigned short dim
){

 const expansion& l1 = expansion_sq_dist(p1,p0,dim);
 const expansion& l2 = expansion_sq_dist(p2,p0,dim);
 const expansion& l3 = expansion_sq_dist(p3,p0,dim);
 const expansion& l4 = expansion_sq_dist(p4,p0,dim);
 const expansion& l5 = expansion_sq_dist(p5,p0,dim);


 const expansion& a10 = expansion_dot_at(p1,q0,p0,dim).scale_fast(2.0);
 const expansion& a11 = expansion_dot_at(p1,q1,p0,dim).scale_fast(2.0);
 const expansion& a12 = expansion_dot_at(p1,q2,p0,dim).scale_fast(2.0);
 const expansion& a13 = expansion_dot_at(p1,q3,p0,dim).scale_fast(2.0);
 const expansion& a14 = expansion_dot_at(p1,q4,p0,dim).scale_fast(2.0);

 const expansion& a20 = expansion_dot_at(p2,q0,p0,dim).scale_fast(2.0);
 const expansion& a21 = expansion_dot_at(p2,q1,p0,dim).scale_fast(2.0);
 const expansion& a22 = expansion_dot_at(p2,q2,p0,dim).scale_fast(2.0);
 const expansion& a23 = expansion_dot_at(p2,q3,p0,dim).scale_fast(2.0);
 const expansion& a24 = expansion_dot_at(p2,q4,p0,dim).scale_fast(2.0);

 const expansion& a30 = expansion_dot_at(p3,q0,p0,dim).scale_fast(2.0);
 const expansion& a31 = expansion_dot_at(p3,q1,p0,dim).scale_fast(2.0);
 const expansion& a32 = expansion_dot_at(p3,q2,p0,dim).scale_fast(2.0);
 const expansion& a33 = expansion_dot_at(p3,q3,p0,dim).scale_fast(2.0);
 const expansion& a34 = expansion_dot_at(p3,q4,p0,dim).scale_fast(2.0);

 const expansion& a40 = expansion_dot_at(p4,q0,p0,dim).scale_fast(2.0);
 const expansion& a41 = expansion_dot_at(p4,q1,p0,dim).scale_fast(2.0);
 const expansion& a42 = expansion_dot_at(p4,q2,p0,dim).scale_fast(2.0);
 const expansion& a43 = expansion_dot_at(p4,q3,p0,dim).scale_fast(2.0);
 const expansion& a44 = expansion_dot_at(p4,q4,p0,dim).scale_fast(2.0);

 const expansion& a50 = expansion_dot_at(p5,q0,p0,dim).scale_fast(2.0);
 const expansion& a51 = expansion_dot_at(p5,q1,p0,dim).scale_fast(2.0);
 const expansion& a52 = expansion_dot_at(p5,q2,p0,dim).scale_fast(2.0);
 const expansion& a53 = expansion_dot_at(p5,q3,p0,dim).scale_fast(2.0);
 const expansion& a54 = expansion_dot_at(p5,q4,p0,dim).scale_fast(2.0);

 // [ b00 b01 b02 b03 b04 ]           [  1   1   1   1   1  ] -1
 // [ b10 b11 b12 b13 b14 ]           [ a10 a11 a12 a13 a14 ]
 // [ b20 b21 b22 b23 b24 ] = Delta * [ a20 a21 a22 a23 a24 ]
 // [ b30 b31 b32 b33 b34 ]           [ a30 a31 a32 a33 a34 ]
 // [ b40 b41 b42 b43 b44 ]           [ a40 a41 a42 a43 a44 ]

 const expansion& b00 = expansion_det4x4(a11,a12,a13,a14,a21,a22,a23,a24,a31,a32,a33,a34,a41,a42,a43,a44);
 const expansion& b01 = expansion_det_1111_3x4(a21,a22,a23,a24,a31,a32,a33,a34,a41,a42,a43,a44).negate();
 const expansion& b02 = expansion_det_1111_3x4(a11,a12,a13,a14,a31,a32,a33,a34,a41,a42,a43,a44);
 const expansion& b03 = expansion_det_1111_3x4(a11,a12,a13,a14,a21,a22,a23,a24,a41,a42,a43,a44).negate();
 const expansion& b04 = expansion_det_1111_3x4(a11,a12,a13,a14,a21,a22,a23,a24,a31,a32,a33,a34);

 const expansion& b10 = expansion_det4x4(a10,a12,a13,a14,a20,a22,a23,a24,a30,a32,a33,a34,a40,a42,a43,a44).negate();
 const expansion& b11 = expansion_det_1111_3x4(a20,a22,a23,a24,a30,a32,a33,a34,a40,a42,a43,a44);
 const expansion& b12 = expansion_det_1111_3x4(a10,a12,a13,a14,a30,a32,a33,a34,a40,a42,a43,a44).negate();
 const expansion& b13 = expansion_det_1111_3x4(a10,a12,a13,a14,a20,a22,a23,a24,a40,a42,a43,a44);
 const expansion& b14 = expansion_det_1111_3x4(a10,a12,a13,a14,a20,a22,a23,a24,a30,a32,a33,a34).negate();

 const expansion& b20 = expansion_det4x4(a10,a11,a13,a14,a20,a21,a23,a24,a30,a31,a33,a34,a40,a41,a43,a44);
 const expansion& b21 = expansion_det_1111_3x4(a20,a21,a23,a24,a30,a31,a33,a34,a40,a41,a43,a44).negate();
 const expansion& b22 = expansion_det_1111_3x4(a10,a11,a13,a14,a30,a31,a33,a34,a40,a41,a43,a44);
 const expansion& b23 = expansion_det_1111_3x4(a10,a11,a13,a14,a20,a21,a23,a24,a40,a41,a43,a44).negate();
 const expansion& b24 = expansion_det_1111_3x4(a10,a11,a13,a14,a20,a21,a23,a24,a30,a31,a33,a34);

 const expansion& b30 = expansion_det4x4(a10,a11,a12,a14,a20,a21,a22,a24,a30,a31,a32,a34,a40,a41,a42,a44).negate();
 const expansion& b31 = expansion_det_1111_3x4(a20,a21,a22,a24,a30,a31,a32,a34,a40,a41,a42,a44);
 const expansion& b32 = expansion_det_1111_3x4(a10,a11,a12,a14,a30,a31,a32,a34,a40,a41,a42,a44).negate();
 const expansion& b33 = expansion_det_1111_3x4(a10,a11,a12,a14,a20,a21,a22,a24,a40,a41,a42,a44);
 const expansion& b34 = expansion_det_1111_3x4(a10,a11,a12,a14,a20,a21,a22,a24,a30,a31,a32,a34).negate();

 const expansion& b40 = expansion_det4x4(a10,a11,a12,a13,a20,a21,a22,a23,a30,a31,a32,a33,a40,a41,a42,a43);
 const expansion& b41 = expansion_det_1111_3x4(a20,a21,a22,a23,a30,a31,a32,a33,a40,a41,a42,a43).negate();
 const expansion& b42 = expansion_det_1111_3x4(a10,a11,a12,a13,a30,a31,a32,a33,a40,a41,a42,a43);
 const expansion& b43 = expansion_det_1111_3x4(a10,a11,a12,a13,a20,a21,a22,a23,a40,a41,a42,a43).negate();
 const expansion& b44 = expansion_det_1111_3x4(a10,a11,a12,a13,a20,a21,a22,a23,a30,a31,a32,a33);

// expansion Delta = b00+b10+b20+b30+b40;
 expansion& Delta = expansion_sum( expansion_sum4(b00,b10,b20,b30),b40);
 Sign Delta_sign = Delta.sign();
 avro_assert( Delta_sign!=GEO::ZERO );

 //       [ Lambda0 ]   [ b01 b02 b03 b04 ]   [ l1 ]   [ b00 ]
 //       [ Lambda1 ]   [ b11 b12 b13 b14 ]   [ l2 ]   [ b10 ]
 // Delta [ Lambda2 ] = [ b21 b22 b23 b24 ] * [ l3 ] + [ b20 ]
 //       [ Lambda3 ]   [ b31 b32 b33 b34 ]   [ l4 ]   [ b30 ]
 //       [ Lambda4 ]   [ b41 b42 b43 b44 ]   [ l5 ]   [ b40 ]

 const expansion& b01_l1 = expansion_product(b01,l1);
 const expansion& b02_l2 = expansion_product(b02,l2);
 const expansion& b03_l3 = expansion_product(b03,l3);
 const expansion& b04_l4 = expansion_product(b04,l4);
 const expansion& DeltaLambda0 = expansion_sum( expansion_sum4(b01_l1,b02_l2,b03_l3,b04_l4),b00);

 const expansion& b11_l1 = expansion_product(b11,l1);
 const expansion& b12_l2 = expansion_product(b12,l2);
 const expansion& b13_l3 = expansion_product(b13,l3);
 const expansion& b14_l4 = expansion_product(b14,l4);
 const expansion& DeltaLambda1 = expansion_sum( expansion_sum4(b11_l1,b12_l2,b13_l3,b14_l4),b10);

 const expansion& b21_l1 = expansion_product(b21,l1);
 const expansion& b22_l2 = expansion_product(b22,l2);
 const expansion& b23_l3 = expansion_product(b23,l3);
 const expansion& b24_l4 = expansion_product(b24,l4);
 const expansion& DeltaLambda2 = expansion_sum( expansion_sum4(b21_l1,b22_l2,b23_l3,b24_l4),b20);

 const expansion& b31_l1 = expansion_product(b31,l1);
 const expansion& b32_l2 = expansion_product(b32,l2);
 const expansion& b33_l3 = expansion_product(b33,l3);
 const expansion& b34_l4 = expansion_product(b34,l4);
 const expansion& DeltaLambda3 = expansion_sum( expansion_sum4(b31_l1,b32_l2,b33_l3,b34_l4),b30);

 const expansion& b41_l1 = expansion_product(b41,l1);
 const expansion& b42_l2 = expansion_product(b42,l2);
 const expansion& b43_l3 = expansion_product(b43,l3);
 const expansion& b44_l4 = expansion_product(b44,l4);
 const expansion& DeltaLambda4 = expansion_sum( expansion_sum4(b41_l1,b42_l2,b43_l3,b44_l4),b40);

 // r = Delta*l5 - ( a50*DeltaLambda0 +a51*DeltaLambda1 +a52*DeltaLambda2 + a53*DeltaLambda3 + a54*DeltaLambda4 )
 const expansion& r0 = expansion_product(Delta,l5);
 const expansion& r1 = expansion_product(a50, DeltaLambda0);
 const expansion& r2 = expansion_product(a51, DeltaLambda1);
 const expansion& r3 = expansion_product(a52, DeltaLambda2);
 const expansion& r4 = expansion_product(a53, DeltaLambda3);
 const expansion& r5 = expansion_product(a54, DeltaLambda4);
 const expansion& r12345 = expansion_sum( r1 , expansion_sum4(r2,r3,r4,r5) );
 const expansion& r = expansion_diff( r0 , r12345 );
 Sign r_sign = r.sign();

 if (r_sign==ZERO)
 {
   const double* p_sort[6];
   p_sort[0] = p0;
   p_sort[1] = p1;
   p_sort[2] = p2;
   p_sort[3] = p3;
   p_sort[4] = p4;
   p_sort[5] = p5;
   std::sort(p_sort,p_sort+6);
   for (index_t i=0;i<6;++i) {
     if(p_sort[i] == p0) {
      //Sign result = Sign( Delta_sign*(Delta -( (b01+b02+b03+b04)*a50+ (b11+b12+b13+b14)*a51+ (b21+b22+b23+b24)*a52+ (b31+b32+b33+b34)*a53+ (b41+b42+b43+b44)*a54 )).sign());

      const expansion& z0 = expansion_product(a50,expansion_sum4(b01,b02,b03,b04));
      const expansion& z1 = expansion_product(a51,expansion_sum4(b11,b12,b13,b14));
      const expansion& z2 = expansion_product(a52,expansion_sum4(b21,b22,b23,b24));
      const expansion& z3 = expansion_product(a53,expansion_sum4(b31,b32,b33,b34));
      const expansion& z4 = expansion_product(a54,expansion_sum4(b41,b42,b43,b44));
      const expansion& z =  expansion_diff( Delta, expansion_sum4( z0,z1,z2, expansion_sum(z3,z4) ) );

      Sign z_sign = z.sign();
      if (z_sign!=ZERO) return Sign(Delta_sign*z_sign);

      //Sign result = Sign( Delta_sign*(z.sign() ));
      //if(result != ZERO) { return result ; }
    }
    if(p_sort[i] == p1) {
      //Sign result = Sign( Delta_sign*(a50*b01+a51*b11+a52*b21+a53*b31+a54*b41).sign());

      const expansion& z0 = expansion_product(a50,b01);
      const expansion& z1 = expansion_product(a51,b11);
      const expansion& z2 = expansion_product(a52,b21);
      const expansion& z3 = expansion_product(a53,b31);
      const expansion& z4 = expansion_product(a54,b41);
      const expansion& z = expansion_sum4( z0,z1,z2, expansion_sum(z3,z4) );

      Sign z_sign = z.sign();
      if (z_sign!=ZERO) return Sign(Delta_sign*z_sign);

      //Sign result = Sign( Delta_sign*( z.sign() ) );
      //if(result != ZERO) { return result ; }
    }
    if(p_sort[i] == p2) {
      //Sign result = Sign( Delta_sign*(a50*b02+a51*b12+a52*b22+a53*b32+a54*b42).sign());

      const expansion& z0 = expansion_product(a50,b02);
      const expansion& z1 = expansion_product(a51,b12);
      const expansion& z2 = expansion_product(a52,b22);
      const expansion& z3 = expansion_product(a53,b32);
      const expansion& z4 = expansion_product(a54,b42);
      const expansion& z = expansion_sum4( z0,z1,z2, expansion_sum(z3,z4) );

      Sign z_sign = z.sign();
      if (z_sign!=ZERO) return Sign(Delta_sign*z_sign);

      //Sign result = Sign( Delta_sign*( z.sign() ) );
      //if(result != ZERO) { return result ; }
    }
    if(p_sort[i] == p3) {
      //Sign result = Sign( Delta_sign*(a50*b03+a51*b13+a52*b23+a53*b33+a54*b43).sign());

      const expansion& z0 = expansion_product(a50,b03);
      const expansion& z1 = expansion_product(a51,b13);
      const expansion& z2 = expansion_product(a52,b23);
      const expansion& z3 = expansion_product(a53,b33);
      const expansion& z4 = expansion_product(a54,b43);
      const expansion& z = expansion_sum4( z0,z1,z2, expansion_sum(z3,z4) );

      Sign z_sign = z.sign();
      if (z_sign!=ZERO) return Sign(Delta_sign*z_sign);

      //Sign result = Sign( Delta_sign*( z.sign() ) );
      //if(result != ZERO) { return result ; }
    }
    if(p_sort[i] == p4) {
      //Sign result = Sign( Delta_sign*(a50*b04+a51*b14+a52*b24+a53*b34+a54*b44).sign());

      const expansion& z0 = expansion_product(a50,b04);
      const expansion& z1 = expansion_product(a51,b14);
      const expansion& z2 = expansion_product(a52,b24);
      const expansion& z3 = expansion_product(a53,b34);
      const expansion& z4 = expansion_product(a54,b44);
      const expansion& z = expansion_sum4( z0,z1,z2, expansion_sum(z3,z4) );

      Sign z_sign = z.sign();
      if (z_sign!=ZERO) return Sign(Delta_sign*z_sign);

      //Sign result = Sign( Delta_sign*( z.sign() ) );
      //if(result != ZERO) { return result ; }
    }

    if(p_sort[i] == p5) {
      return NEGATIVE;
    }
  } // loop over all p_sort
 }
 return Sign( r_sign*Delta_sign );
}

} // PCK

} // GEO
