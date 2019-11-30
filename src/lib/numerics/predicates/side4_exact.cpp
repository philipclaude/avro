// avro: Adaptive Voronoi Remesher
// Copyright 2017-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "numerics/expansion.h"
#include "numerics/predicates.h"

#include "numerics/predicates/geo_side_exact.h"

namespace GEO {

namespace PCK {


Sign avro_side4_nd_exact_pck(const double* p0,const double* p1,const double* p2,const double* p3,const double* p4,
 const double* q0,const double* q1,const double* q2,const double* q3 ,unsigned short dim )
{

    const expansion& l1 = expansion_sq_dist(p1, p0, dim);
    const expansion& l2 = expansion_sq_dist(p2, p0, dim);
    const expansion& l3 = expansion_sq_dist(p3, p0, dim);
    const expansion& l4 = expansion_sq_dist(p4, p0, dim);

    const expansion& a10 = expansion_dot_at(p1, q0, p0, dim).scale_fast(2.0);
    const expansion& a11 = expansion_dot_at(p1, q1, p0, dim).scale_fast(2.0);
    const expansion& a12 = expansion_dot_at(p1, q2, p0, dim).scale_fast(2.0);
    const expansion& a13 = expansion_dot_at(p1, q3, p0, dim).scale_fast(2.0);

    const expansion& a20 = expansion_dot_at(p2, q0, p0, dim).scale_fast(2.0);
    const expansion& a21 = expansion_dot_at(p2, q1, p0, dim).scale_fast(2.0);
    const expansion& a22 = expansion_dot_at(p2, q2, p0, dim).scale_fast(2.0);
    const expansion& a23 = expansion_dot_at(p2, q3, p0, dim).scale_fast(2.0);

    const expansion& a30 = expansion_dot_at(p3, q0, p0, dim).scale_fast(2.0);
    const expansion& a31 = expansion_dot_at(p3, q1, p0, dim).scale_fast(2.0);
    const expansion& a32 = expansion_dot_at(p3, q2, p0, dim).scale_fast(2.0);
    const expansion& a33 = expansion_dot_at(p3, q3, p0, dim).scale_fast(2.0);

    const expansion& a40 = expansion_dot_at(p4, q0, p0, dim).scale_fast(2.0);
    const expansion& a41 = expansion_dot_at(p4, q1, p0, dim).scale_fast(2.0);
    const expansion& a42 = expansion_dot_at(p4, q2, p0, dim).scale_fast(2.0);
    const expansion& a43 = expansion_dot_at(p4, q3, p0, dim).scale_fast(2.0);

    // [ b00 b01 b02 b03 ]           [  1   1   1   1  ]-1
    // [ b10 b11 b12 b13 ]           [ a10 a11 a12 a13 ]
    // [ b20 b21 b22 b23 ] = Delta * [ a20 a21 a22 a23 ]
    // [ b30 b31 b32 b33 ]           [ a30 a31 a32 a33 ]

    // Note: we could probably reuse some of the co-factors
    // (but for now I'd rather keep this form that is easier to
    //  read ... and to debug if need be !)

    const expansion& b00 = expansion_det3x3(a11, a12, a13, a21, a22, a23, a31, a32, a33);
    const expansion& b01 = expansion_det_111_2x3(a21, a22, a23, a31, a32, a33).negate();
    const expansion& b02 = expansion_det_111_2x3(a11, a12, a13, a31, a32, a33);
    const expansion& b03 = expansion_det_111_2x3(a11, a12, a13, a21, a22, a23).negate();

    const expansion& b10 = expansion_det3x3(a10, a12, a13, a20, a22, a23, a30, a32, a33).negate();
    const expansion& b11 = expansion_det_111_2x3(a20, a22, a23, a30, a32, a33);
    const expansion& b12 = expansion_det_111_2x3(a10, a12, a13, a30, a32, a33).negate();
    const expansion& b13 = expansion_det_111_2x3(a10, a12, a13, a20, a22, a23);

    const expansion& b20 = expansion_det3x3(a10, a11, a13, a20, a21, a23, a30, a31, a33);
    const expansion& b21 = expansion_det_111_2x3(a20, a21, a23, a30, a31, a33).negate();
    const expansion& b22 = expansion_det_111_2x3(a10, a11, a13, a30, a31, a33);
    const expansion& b23 = expansion_det_111_2x3(a10, a11, a13, a20, a21, a23).negate();

    const expansion& b30 = expansion_det3x3(a10, a11, a12, a20, a21, a22, a30, a31, a32).negate();
    const expansion& b31 = expansion_det_111_2x3(a20, a21, a22, a30, a31, a32);
    const expansion& b32 = expansion_det_111_2x3(a10, a11, a12, a30, a31, a32).negate();
    const expansion& b33 = expansion_det_111_2x3(a10, a11, a12, a20, a21, a22);

    const expansion& Delta = expansion_sum4(b00, b10, b20, b30);
    Sign Delta_sign = Delta.sign();
    geo_assert(Delta_sign != ZERO);

    //       [ Lambda0 ]   [ b01 b02 b03 ]   [ l1 ]   [ b00 ]
    //       [ Lambda1 ]   [ b11 b12 b13 ]   [ l2 ]   [ b10 ]
    // Delta [ Lambda2 ] = [ b21 b22 b23 ] * [ l3 ] + [ b20 ]
    //       [ Lambda3 ]   [ b31 b32 b33 ]   [ l4 ]   [ b30 ]

    const expansion& b01_l1 = expansion_product(b01, l1);
    const expansion& b02_l2 = expansion_product(b02, l2);
    const expansion& b03_l3 = expansion_product(b03, l3);
    const expansion& DeltaLambda0 = expansion_sum4(b01_l1, b02_l2, b03_l3, b00);

    const expansion& b11_l1 = expansion_product(b11, l1);
    const expansion& b12_l2 = expansion_product(b12, l2);
    const expansion& b13_l3 = expansion_product(b13, l3);
    const expansion& DeltaLambda1 = expansion_sum4(b11_l1, b12_l2, b13_l3, b10);

    const expansion& b21_l1 = expansion_product(b21, l1);
    const expansion& b22_l2 = expansion_product(b22, l2);
    const expansion& b23_l3 = expansion_product(b23, l3);
    const expansion& DeltaLambda2 = expansion_sum4(b21_l1, b22_l2, b23_l3, b20);

    const expansion& b31_l1 = expansion_product(b31, l1);
    const expansion& b32_l2 = expansion_product(b32, l2);
    const expansion& b33_l3 = expansion_product(b33, l3);
    const expansion& DeltaLambda3 = expansion_sum4(b31_l1, b32_l2, b33_l3, b30);

    const expansion& r0 = expansion_product(Delta, l4);
    const expansion& r1 = expansion_product(a40, DeltaLambda0);
    const expansion& r2 = expansion_product(a41, DeltaLambda1);
    const expansion& r3 = expansion_product(a42, DeltaLambda2);
    const expansion& r4 = expansion_product(a43, DeltaLambda3);
    const expansion& r1234 = expansion_sum4(r1, r2, r3, r4);
    const expansion& r = expansion_diff(r0, r1234);
    Sign r_sign = r.sign();

    // Simulation of Simplicity (symbolic perturbation)
    if(r_sign == ZERO) {
        const double* p_sort[5];
        p_sort[0] = p0;
        p_sort[1] = p1;
        p_sort[2] = p2;
        p_sort[3] = p3;
        p_sort[4] = p4;
        std::sort(p_sort, p_sort + 5);
        for(index_t i = 0; i < 5; ++i) {
            if(p_sort[i] == p0) {
                const expansion& z1_0 = expansion_sum3(b01, b02, b03);
                const expansion& z1 = expansion_product(a40, z1_0);
                const expansion& z2_0 = expansion_sum3(b11, b12, b13);
                const expansion& z2 = expansion_product(a41, z2_0);
                const expansion& z3_0 = expansion_sum3(b21, b22, b23);
                const expansion& z3 = expansion_product(a42, z3_0);
                const expansion& z4_0 = expansion_sum3(b31, b32, b33);
                const expansion& z4 = expansion_product(a43, z4_0);
                const expansion& z1234 = expansion_sum4(z1, z2, z3, z4);
                const expansion& z = expansion_diff(Delta, z1234);
                Sign z_sign = z.sign();
                len_side4_SOS = geo_max(len_side4_SOS, z.length());
                if(z_sign != ZERO) {
                    return Sign(Delta_sign * z_sign);
                }
            } else if(p_sort[i] == p1) {
                const expansion& z1 = expansion_product(a40, b01);
                const expansion& z2 = expansion_product(a41, b11);
                const expansion& z3 = expansion_product(a42, b21);
                const expansion& z4 = expansion_product(a43, b31);
                const expansion& z = expansion_sum4(z1, z2, z3, z4);
                Sign z_sign = z.sign();
                len_side4_SOS = geo_max(len_side4_SOS, z.length());
                if(z_sign != ZERO) {
                    return Sign(Delta_sign * z_sign);
                }
            } else if(p_sort[i] == p2) {
                const expansion& z1 = expansion_product(a40, b02);
                const expansion& z2 = expansion_product(a41, b12);
                const expansion& z3 = expansion_product(a42, b22);
                const expansion& z4 = expansion_product(a43, b32);
                const expansion& z = expansion_sum4(z1, z2, z3, z4);
                Sign z_sign = z.sign();
                len_side4_SOS = geo_max(len_side4_SOS, z.length());
                if(z_sign != ZERO) {
                    return Sign(Delta_sign * z_sign);
                }
            } else if(p_sort[i] == p3) {
                const expansion& z1 = expansion_product(a40, b03);
                const expansion& z2 = expansion_product(a41, b13);
                const expansion& z3 = expansion_product(a42, b23);
                const expansion& z4 = expansion_product(a43, b33);
                const expansion& z = expansion_sum4(z1, z2, z3, z4);
                Sign z_sign = z.sign();
                len_side4_SOS = geo_max(len_side4_SOS, z.length());
                if(z_sign != ZERO) {
                    return Sign(Delta_sign * z_sign);
                }
            } else if(p_sort[i] == p4) {
                return NEGATIVE;
            }
        }
        geo_assert_not_reached;
    }
    return Sign(r_sign * Delta_sign);
}

} // PCK

} // GEO
