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

#include "numerics/predicates/geo_side_exact.h"

namespace GEO {

namespace PCK {

Sign avro_side4_nd_exact_pck(const double* p0,const double* p1,const double* p2,const double* p3,const double* p4,
 const double* q0,const double* q1,const double* q2,const double* q3 ,unsigned short dim )
{

  printf("side 4 exact!\n");

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

/**
 * \brief Exact implementation of the side4_3d_SOS() predicate
 *  using low-level exact arithmetics API (expansion class).
 * \param[in] sos if true, applies symbolic perturbation when
 *  result is zero, else returns zero
 */
Sign avro_side4_3d_exact_pck(
    const double* p0, const double* p1, const double* p2, const double* p3,
    const double* p4, bool sos
) {

    cnt_side4_exact++;

    const expansion& a11 = expansion_diff(p1[0], p0[0]);
    const expansion& a12 = expansion_diff(p1[1], p0[1]);
    const expansion& a13 = expansion_diff(p1[2], p0[2]);
    const expansion& a14 = expansion_sq_dist(p1, p0, 3).negate();

    const expansion& a21 = expansion_diff(p2[0], p0[0]);
    const expansion& a22 = expansion_diff(p2[1], p0[1]);
    const expansion& a23 = expansion_diff(p2[2], p0[2]);
    const expansion& a24 = expansion_sq_dist(p2, p0, 3).negate();

    const expansion& a31 = expansion_diff(p3[0], p0[0]);
    const expansion& a32 = expansion_diff(p3[1], p0[1]);
    const expansion& a33 = expansion_diff(p3[2], p0[2]);
    const expansion& a34 = expansion_sq_dist(p3, p0, 3).negate();

    const expansion& a41 = expansion_diff(p4[0], p0[0]);
    const expansion& a42 = expansion_diff(p4[1], p0[1]);
    const expansion& a43 = expansion_diff(p4[2], p0[2]);
    const expansion& a44 = expansion_sq_dist(p4, p0, 3).negate();

    // This commented-out version does not reuse
    // the 2x2 minors.
/*
    const expansion& Delta1 = expansion_det3x3(
        a21, a22, a23,
        a31, a32, a33,
        a41, a42, a43
    );
    const expansion& Delta2 = expansion_det3x3(
        a11, a12, a13,
        a31, a32, a33,
        a41, a42, a43
    );
    const expansion& Delta3 = expansion_det3x3(
        a11, a12, a13,
        a21, a22, a23,
        a41, a42, a43
    );
    const expansion& Delta4 = expansion_det3x3(
        a11, a12, a13,
        a21, a22, a23,
        a31, a32, a33
    );
*/

    // Optimized version that reuses the 2x2 minors
    const expansion& m12 = expansion_det2x2(a12,a13,a22,a23);
    const expansion& m13 = expansion_det2x2(a12,a13,a32,a33);
    const expansion& m14 = expansion_det2x2(a12,a13,a42,a43);
    const expansion& m23 = expansion_det2x2(a22,a23,a32,a33);
    const expansion& m24 = expansion_det2x2(a22,a23,a42,a43);
    const expansion& m34 = expansion_det2x2(a32,a33,a42,a43);


    const expansion& z11 = expansion_product(a21,m34);
    const expansion& z12 = expansion_product(a31,m24).negate();
    const expansion& z13 = expansion_product(a41,m23);
    const expansion& Delta1 = expansion_sum3(z11,z12,z13);

    const expansion& z21 = expansion_product(a11,m34);
    const expansion& z22 = expansion_product(a31,m14).negate();
    const expansion& z23 = expansion_product(a41,m13);
    const expansion& Delta2 = expansion_sum3(z21,z22,z23);

    const expansion& z31 = expansion_product(a11,m24);
    const expansion& z32 = expansion_product(a21,m14).negate();
    const expansion& z33 = expansion_product(a41,m12);
    const expansion& Delta3 = expansion_sum3(z31,z32,z33);

    const expansion& z41 = expansion_product(a11,m23);
    const expansion& z42 = expansion_product(a21,m13).negate();
    const expansion& z43 = expansion_product(a31,m12);
    const expansion& Delta4 = expansion_sum3(z41,z42,z43);


    Sign Delta4_sign = Delta4.sign();
    geo_assert(Delta4_sign != ZERO);

    const expansion& r_1 = expansion_product(Delta1, a14);
    const expansion& r_2 = expansion_product(Delta2, a24).negate();
    const expansion& r_3 = expansion_product(Delta3, a34);
    const expansion& r_4 = expansion_product(Delta4, a44).negate();
    const expansion& r = expansion_sum4(r_1, r_2, r_3, r_4);
    Sign r_sign = r.sign();

    // Statistics
    len_side4_num = std::max(len_side4_num, r.length());
    len_side4_denom = std::max(len_side4_denom, Delta1.length());

    // Simulation of Simplicity (symbolic perturbation)
    if(sos && r_sign == ZERO) {
        cnt_side4_SOS++;
        const double* p_sort[5];
        p_sort[0] = p0;
        p_sort[1] = p1;
        p_sort[2] = p2;
        p_sort[3] = p3;
        p_sort[4] = p4;
        //SOS_sort(p_sort, p_sort + 5, 3);
        std::sort(p_sort,p_sort+5);
        for(index_t i = 0; i < 5; ++i) {
            if(p_sort[i] == p0) {
                const expansion& z1 = expansion_diff(Delta2, Delta1);
                const expansion& z2 = expansion_diff(Delta4, Delta3);
                const expansion& z = expansion_sum(z1, z2);
                Sign z_sign = z.sign();
                len_side4_SOS = std::max(len_side4_SOS, z.length());
                if(z_sign != ZERO) {
                    return Sign(Delta4_sign * z_sign);
                }
            } else if(p_sort[i] == p1) {
                Sign Delta1_sign = Delta1.sign();
                if(Delta1_sign != ZERO) {
                    len_side4_SOS = std::max(len_side4_SOS, Delta1.length());
                    return Sign(Delta4_sign * Delta1_sign);
                }
            } else if(p_sort[i] == p2) {
                Sign Delta2_sign = Delta2.sign();
                if(Delta2_sign != ZERO) {
                    len_side4_SOS = std::max(len_side4_SOS, Delta2.length());
                    return Sign(-Delta4_sign * Delta2_sign);
                }
            } else if(p_sort[i] == p3) {
                Sign Delta3_sign = Delta3.sign();
                if(Delta3_sign != ZERO) {
                    len_side4_SOS = std::max(len_side4_SOS, Delta3.length());
                    return Sign(Delta4_sign * Delta3_sign);
                }
            } else if(p_sort[i] == p4) {
                return NEGATIVE;
            }
        }
    }
    return Sign(Delta4_sign * r_sign);
}

} // PCK

} // GEO
