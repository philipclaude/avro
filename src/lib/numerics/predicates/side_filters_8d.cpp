// avro: Adaptive Voronoi Remesher
// Copyright 2017-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "numerics/predicates.h"

namespace GEO {

namespace PCK {

int avro_side1_8d_filter( const double* p0, const double* p1, const double* q0) {
    double p0_0_p1_0 = (p0[0] - p1[0]);
    double p0_1_p1_1 = (p0[1] - p1[1]);
    double p0_2_p1_2 = (p0[2] - p1[2]);
    double p0_3_p1_3 = (p0[3] - p1[3]);
    double p0_4_p1_4 = (p0[4] - p1[4]);
    double p0_5_p1_5 = (p0[5] - p1[5]);
    double p0_6_p1_6 = (p0[6] - p1[6]);
    double p0_7_p1_7 = (p0[7] - p1[7]);
    double r;
    r = (1 * ((((((((p0_0_p1_0 * p0_0_p1_0) + (p0_1_p1_1 * p0_1_p1_1)) + (p0_2_p1_2 * p0_2_p1_2)) + (p0_3_p1_3 * p0_3_p1_3)) + (p0_4_p1_4 * p0_4_p1_4)) + (p0_5_p1_5 * p0_5_p1_5)) + (p0_6_p1_6 * p0_6_p1_6)) + (p0_7_p1_7 * p0_7_p1_7)));
    double p1_0_p0_0 = (p1[0] - p0[0]);
    double q0_0_p0_0 = (q0[0] - p0[0]);
    double p1_1_p0_1 = (p1[1] - p0[1]);
    double q0_1_p0_1 = (q0[1] - p0[1]);
    double p1_2_p0_2 = (p1[2] - p0[2]);
    double q0_2_p0_2 = (q0[2] - p0[2]);
    double p1_3_p0_3 = (p1[3] - p0[3]);
    double q0_3_p0_3 = (q0[3] - p0[3]);
    double p1_4_p0_4 = (p1[4] - p0[4]);
    double q0_4_p0_4 = (q0[4] - p0[4]);
    double p1_5_p0_5 = (p1[5] - p0[5]);
    double q0_5_p0_5 = (q0[5] - p0[5]);
    double p1_6_p0_6 = (p1[6] - p0[6]);
    double q0_6_p0_6 = (q0[6] - p0[6]);
    double p1_7_p0_7 = (p1[7] - p0[7]);
    double q0_7_p0_7 = (q0[7] - p0[7]);
    r = (r - (2 * ((((((((p1_0_p0_0 * q0_0_p0_0) + (p1_1_p0_1 * q0_1_p0_1)) + (p1_2_p0_2 * q0_2_p0_2)) + (p1_3_p0_3 * q0_3_p0_3)) + (p1_4_p0_4 * q0_4_p0_4)) + (p1_5_p0_5 * q0_5_p0_5)) + (p1_6_p0_6 * q0_6_p0_6)) + (p1_7_p0_7 * q0_7_p0_7))));
    int int_tmp_result;
    double eps;
    double max1 = fabs(p0_2_p1_2);
    if( (max1 < fabs(p0_0_p1_0)) )
    {
        max1 = fabs(p0_0_p1_0);
    }
    if( (max1 < fabs(p0_1_p1_1)) )
    {
        max1 = fabs(p0_1_p1_1);
    }
    if( (max1 < fabs(p0_3_p1_3)) )
    {
        max1 = fabs(p0_3_p1_3);
    }
    if( (max1 < fabs(p0_4_p1_4)) )
    {
        max1 = fabs(p0_4_p1_4);
    }
    if( (max1 < fabs(p0_5_p1_5)) )
    {
        max1 = fabs(p0_5_p1_5);
    }
    if( (max1 < fabs(p0_6_p1_6)) )
    {
        max1 = fabs(p0_6_p1_6);
    }
    if( (max1 < fabs(p0_7_p1_7)) )
    {
        max1 = fabs(p0_7_p1_7);
    }
    if( (max1 < fabs(p1_0_p0_0)) )
    {
        max1 = fabs(p1_0_p0_0);
    }
    if( (max1 < fabs(p1_1_p0_1)) )
    {
        max1 = fabs(p1_1_p0_1);
    }
    if( (max1 < fabs(p1_2_p0_2)) )
    {
        max1 = fabs(p1_2_p0_2);
    }
    if( (max1 < fabs(p1_3_p0_3)) )
    {
        max1 = fabs(p1_3_p0_3);
    }
    if( (max1 < fabs(p1_4_p0_4)) )
    {
        max1 = fabs(p1_4_p0_4);
    }
    if( (max1 < fabs(p1_5_p0_5)) )
    {
        max1 = fabs(p1_5_p0_5);
    }
    if( (max1 < fabs(p1_6_p0_6)) )
    {
        max1 = fabs(p1_6_p0_6);
    }
    if( (max1 < fabs(p1_7_p0_7)) )
    {
        max1 = fabs(p1_7_p0_7);
    }
    double max2 = fabs(p0_2_p1_2);
    if( (max2 < fabs(p0_0_p1_0)) )
    {
        max2 = fabs(p0_0_p1_0);
    }
    if( (max2 < fabs(p0_1_p1_1)) )
    {
        max2 = fabs(p0_1_p1_1);
    }
    if( (max2 < fabs(p0_3_p1_3)) )
    {
        max2 = fabs(p0_3_p1_3);
    }
    if( (max2 < fabs(p0_4_p1_4)) )
    {
        max2 = fabs(p0_4_p1_4);
    }
    if( (max2 < fabs(p0_5_p1_5)) )
    {
        max2 = fabs(p0_5_p1_5);
    }
    if( (max2 < fabs(p0_6_p1_6)) )
    {
        max2 = fabs(p0_6_p1_6);
    }
    if( (max2 < fabs(p0_7_p1_7)) )
    {
        max2 = fabs(p0_7_p1_7);
    }
    if( (max2 < fabs(q0_0_p0_0)) )
    {
        max2 = fabs(q0_0_p0_0);
    }
    if( (max2 < fabs(q0_1_p0_1)) )
    {
        max2 = fabs(q0_1_p0_1);
    }
    if( (max2 < fabs(q0_2_p0_2)) )
    {
        max2 = fabs(q0_2_p0_2);
    }
    if( (max2 < fabs(q0_3_p0_3)) )
    {
        max2 = fabs(q0_3_p0_3);
    }
    if( (max2 < fabs(q0_4_p0_4)) )
    {
        max2 = fabs(q0_4_p0_4);
    }
    if( (max2 < fabs(q0_5_p0_5)) )
    {
        max2 = fabs(q0_5_p0_5);
    }
    if( (max2 < fabs(q0_6_p0_6)) )
    {
        max2 = fabs(q0_6_p0_6);
    }
    if( (max2 < fabs(q0_7_p0_7)) )
    {
        max2 = fabs(q0_7_p0_7);
    }
    double lower_bound_1;
    double upper_bound_1;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    }
    else
    {
        if( (max2 > upper_bound_1) )
        {
            upper_bound_1 = max2;
        }
    }
    if( (lower_bound_1 < 1.15542931091530087067e-147) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 5.59936185544450928309e+101) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (1.66670090166682227006e-14 * (max1 * max2));
        if( (r > eps) )
        {
            int_tmp_result = 1;
        }
        else
        {
            if( (r < -eps) )
            {
                int_tmp_result = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    return int_tmp_result;
}


int avro_side2_8d_filter( const double* p0, const double* p1, const double* p2, const double* q0, const double* q1) {
    double p1_0_p0_0 = (p1[0] - p0[0]);
    double p1_1_p0_1 = (p1[1] - p0[1]);
    double p1_2_p0_2 = (p1[2] - p0[2]);
    double p1_3_p0_3 = (p1[3] - p0[3]);
    double p1_4_p0_4 = (p1[4] - p0[4]);
    double p1_5_p0_5 = (p1[5] - p0[5]);
    double p1_6_p0_6 = (p1[6] - p0[6]);
    double p1_7_p0_7 = (p1[7] - p0[7]);
    double l1;
    l1 = (1 * ((((((((p1_0_p0_0 * p1_0_p0_0) + (p1_1_p0_1 * p1_1_p0_1)) + (p1_2_p0_2 * p1_2_p0_2)) + (p1_3_p0_3 * p1_3_p0_3)) + (p1_4_p0_4 * p1_4_p0_4)) + (p1_5_p0_5 * p1_5_p0_5)) + (p1_6_p0_6 * p1_6_p0_6)) + (p1_7_p0_7 * p1_7_p0_7)));
    double p2_0_p0_0 = (p2[0] - p0[0]);
    double p2_1_p0_1 = (p2[1] - p0[1]);
    double p2_2_p0_2 = (p2[2] - p0[2]);
    double p2_3_p0_3 = (p2[3] - p0[3]);
    double p2_4_p0_4 = (p2[4] - p0[4]);
    double p2_5_p0_5 = (p2[5] - p0[5]);
    double p2_6_p0_6 = (p2[6] - p0[6]);
    double p2_7_p0_7 = (p2[7] - p0[7]);
    double l2;
    l2 = (1 * ((((((((p2_0_p0_0 * p2_0_p0_0) + (p2_1_p0_1 * p2_1_p0_1)) + (p2_2_p0_2 * p2_2_p0_2)) + (p2_3_p0_3 * p2_3_p0_3)) + (p2_4_p0_4 * p2_4_p0_4)) + (p2_5_p0_5 * p2_5_p0_5)) + (p2_6_p0_6 * p2_6_p0_6)) + (p2_7_p0_7 * p2_7_p0_7)));
    double q0_0_p0_0 = (q0[0] - p0[0]);
    double q0_1_p0_1 = (q0[1] - p0[1]);
    double q0_2_p0_2 = (q0[2] - p0[2]);
    double q0_3_p0_3 = (q0[3] - p0[3]);
    double q0_4_p0_4 = (q0[4] - p0[4]);
    double q0_5_p0_5 = (q0[5] - p0[5]);
    double q0_6_p0_6 = (q0[6] - p0[6]);
    double q0_7_p0_7 = (q0[7] - p0[7]);
    double a10;
    a10 = (2 * ((((((((p1_0_p0_0 * q0_0_p0_0) + (p1_1_p0_1 * q0_1_p0_1)) + (p1_2_p0_2 * q0_2_p0_2)) + (p1_3_p0_3 * q0_3_p0_3)) + (p1_4_p0_4 * q0_4_p0_4)) + (p1_5_p0_5 * q0_5_p0_5)) + (p1_6_p0_6 * q0_6_p0_6)) + (p1_7_p0_7 * q0_7_p0_7)));
    double q1_0_p0_0 = (q1[0] - p0[0]);
    double q1_1_p0_1 = (q1[1] - p0[1]);
    double q1_2_p0_2 = (q1[2] - p0[2]);
    double q1_3_p0_3 = (q1[3] - p0[3]);
    double q1_4_p0_4 = (q1[4] - p0[4]);
    double q1_5_p0_5 = (q1[5] - p0[5]);
    double q1_6_p0_6 = (q1[6] - p0[6]);
    double q1_7_p0_7 = (q1[7] - p0[7]);
    double a11;
    a11 = (2 * ((((((((p1_0_p0_0 * q1_0_p0_0) + (p1_1_p0_1 * q1_1_p0_1)) + (p1_2_p0_2 * q1_2_p0_2)) + (p1_3_p0_3 * q1_3_p0_3)) + (p1_4_p0_4 * q1_4_p0_4)) + (p1_5_p0_5 * q1_5_p0_5)) + (p1_6_p0_6 * q1_6_p0_6)) + (p1_7_p0_7 * q1_7_p0_7)));
    double a20;
    a20 = (2 * ((((((((p2_0_p0_0 * q0_0_p0_0) + (p2_1_p0_1 * q0_1_p0_1)) + (p2_2_p0_2 * q0_2_p0_2)) + (p2_3_p0_3 * q0_3_p0_3)) + (p2_4_p0_4 * q0_4_p0_4)) + (p2_5_p0_5 * q0_5_p0_5)) + (p2_6_p0_6 * q0_6_p0_6)) + (p2_7_p0_7 * q0_7_p0_7)));
    double a21;
    a21 = (2 * ((((((((p2_0_p0_0 * q1_0_p0_0) + (p2_1_p0_1 * q1_1_p0_1)) + (p2_2_p0_2 * q1_2_p0_2)) + (p2_3_p0_3 * q1_3_p0_3)) + (p2_4_p0_4 * q1_4_p0_4)) + (p2_5_p0_5 * q1_5_p0_5)) + (p2_6_p0_6 * q1_6_p0_6)) + (p2_7_p0_7 * q1_7_p0_7)));
    double Delta;
    Delta = (a11 - a10);
    double DeltaLambda0;
    DeltaLambda0 = (a11 - l1);
    double DeltaLambda1;
    DeltaLambda1 = (l1 - a10);
    double r;
    r = (((Delta * l2) - (a20 * DeltaLambda0)) - (a21 * DeltaLambda1));
    double eps;
    double max1 = fabs(p1_2_p0_2);
    if( (max1 < fabs(p1_1_p0_1)) )
    {
        max1 = fabs(p1_1_p0_1);
    }
    if( (max1 < fabs(p1_0_p0_0)) )
    {
        max1 = fabs(p1_0_p0_0);
    }
    if( (max1 < fabs(p1_4_p0_4)) )
    {
        max1 = fabs(p1_4_p0_4);
    }
    if( (max1 < fabs(p1_3_p0_3)) )
    {
        max1 = fabs(p1_3_p0_3);
    }
    if( (max1 < fabs(p1_5_p0_5)) )
    {
        max1 = fabs(p1_5_p0_5);
    }
    if( (max1 < fabs(p1_6_p0_6)) )
    {
        max1 = fabs(p1_6_p0_6);
    }
    if( (max1 < fabs(p1_7_p0_7)) )
    {
        max1 = fabs(p1_7_p0_7);
    }
    double max2 = fabs(q0_0_p0_0);
    if( (max2 < fabs(q0_1_p0_1)) )
    {
        max2 = fabs(q0_1_p0_1);
    }
    if( (max2 < fabs(q0_2_p0_2)) )
    {
        max2 = fabs(q0_2_p0_2);
    }
    if( (max2 < fabs(q0_3_p0_3)) )
    {
        max2 = fabs(q0_3_p0_3);
    }
    if( (max2 < fabs(q0_4_p0_4)) )
    {
        max2 = fabs(q0_4_p0_4);
    }
    if( (max2 < fabs(q0_5_p0_5)) )
    {
        max2 = fabs(q0_5_p0_5);
    }
    if( (max2 < fabs(q0_6_p0_6)) )
    {
        max2 = fabs(q0_6_p0_6);
    }
    if( (max2 < fabs(q0_7_p0_7)) )
    {
        max2 = fabs(q0_7_p0_7);
    }
    if( (max2 < fabs(q1_0_p0_0)) )
    {
        max2 = fabs(q1_0_p0_0);
    }
    if( (max2 < fabs(q1_1_p0_1)) )
    {
        max2 = fabs(q1_1_p0_1);
    }
    if( (max2 < fabs(q1_2_p0_2)) )
    {
        max2 = fabs(q1_2_p0_2);
    }
    if( (max2 < fabs(q1_3_p0_3)) )
    {
        max2 = fabs(q1_3_p0_3);
    }
    if( (max2 < fabs(q1_4_p0_4)) )
    {
        max2 = fabs(q1_4_p0_4);
    }
    if( (max2 < fabs(q1_5_p0_5)) )
    {
        max2 = fabs(q1_5_p0_5);
    }
    if( (max2 < fabs(q1_6_p0_6)) )
    {
        max2 = fabs(q1_6_p0_6);
    }
    if( (max2 < fabs(q1_7_p0_7)) )
    {
        max2 = fabs(q1_7_p0_7);
    }
    double lower_bound_1;
    double upper_bound_1;
    double Delta_sign;
    int int_tmp_result;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    }
    else
    {
        if( (max2 > upper_bound_1) )
        {
            upper_bound_1 = max2;
        }
    }
    if( (lower_bound_1 < 1.15542931091530087067e-147) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 1.87072209578355531992e+50) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (1.66670090166682227006e-14 * (max1 * max2));
        if( (Delta > eps) )
        {
            int_tmp_result = 1;
        }
        else
        {
            if( (Delta < -eps) )
            {
                int_tmp_result = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    Delta_sign = int_tmp_result;
    double max3 = max1;
    if( (max3 < max2) )
    {
        max3 = max2;
    }
    double max4 = max2;
    if( (max4 < fabs(p2_3_p0_3)) )
    {
        max4 = fabs(p2_3_p0_3);
    }
    if( (max4 < fabs(p2_2_p0_2)) )
    {
        max4 = fabs(p2_2_p0_2);
    }
    if( (max4 < fabs(p2_0_p0_0)) )
    {
        max4 = fabs(p2_0_p0_0);
    }
    if( (max4 < fabs(p2_1_p0_1)) )
    {
        max4 = fabs(p2_1_p0_1);
    }
    if( (max4 < fabs(p2_4_p0_4)) )
    {
        max4 = fabs(p2_4_p0_4);
    }
    if( (max4 < fabs(p2_5_p0_5)) )
    {
        max4 = fabs(p2_5_p0_5);
    }
    if( (max4 < fabs(p2_6_p0_6)) )
    {
        max4 = fabs(p2_6_p0_6);
    }
    if( (max4 < fabs(p2_7_p0_7)) )
    {
        max4 = fabs(p2_7_p0_7);
    }
    if( (max3 < max4) )
    {
        max3 = max4;
    }
    double r_sign;
    int int_tmp_result_FFWKCAA;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if( (max3 < lower_bound_1) )
    {
        lower_bound_1 = max3;
    }
    else
    {
        if( (max3 > upper_bound_1) )
        {
            upper_bound_1 = max3;
        }
    }
    if( (max4 < lower_bound_1) )
    {
        lower_bound_1 = max4;
    }
    if( (lower_bound_1 < 1.26419510663115923609e-74) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 1.87072209578355531992e+50) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (8.71140112255785451890e-13 * (((max1 * max4) * max4) * max3));
        if( (r > eps) )
        {
            int_tmp_result_FFWKCAA = 1;
        }
        else
        {
            if( (r < -eps) )
            {
                int_tmp_result_FFWKCAA = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    r_sign = int_tmp_result_FFWKCAA;
    int int_tmp_result_k60Ocge;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if( (max3 < lower_bound_1) )
    {
        lower_bound_1 = max3;
    }
    else
    {
        if( (max3 > upper_bound_1) )
        {
            upper_bound_1 = max3;
        }
    }
    if( (max4 < lower_bound_1) )
    {
        lower_bound_1 = max4;
    }
    if( (lower_bound_1 < 1.26419510663115923609e-74) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 1.87072209578355531992e+50) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (8.71140112255785451890e-13 * (((max1 * max4) * max4) * max3));
        if( (r > eps) )
        {
            int_tmp_result_k60Ocge = 1;
        }
        else
        {
            if( (r < -eps) )
            {
                int_tmp_result_k60Ocge = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    //return int_tmp_result_k60Ocge;
    r_sign = int_tmp_result_k60Ocge;
    return (Delta_sign*r_sign);
}


int avro_side3_8d_filter( const double* p0, const double* p1, const double* p2, const double* p3, const double* q0, const double* q1, const double* q2) {
    double p1_0_p0_0 = (p1[0] - p0[0]);
    double p1_1_p0_1 = (p1[1] - p0[1]);
    double p1_2_p0_2 = (p1[2] - p0[2]);
    double p1_3_p0_3 = (p1[3] - p0[3]);
    double p1_4_p0_4 = (p1[4] - p0[4]);
    double p1_5_p0_5 = (p1[5] - p0[5]);
    double p1_6_p0_6 = (p1[6] - p0[6]);
    double p1_7_p0_7 = (p1[7] - p0[7]);
    double l1;
    l1 = (1 * ((((((((p1_0_p0_0 * p1_0_p0_0) + (p1_1_p0_1 * p1_1_p0_1)) + (p1_2_p0_2 * p1_2_p0_2)) + (p1_3_p0_3 * p1_3_p0_3)) + (p1_4_p0_4 * p1_4_p0_4)) + (p1_5_p0_5 * p1_5_p0_5)) + (p1_6_p0_6 * p1_6_p0_6)) + (p1_7_p0_7 * p1_7_p0_7)));
    double p2_0_p0_0 = (p2[0] - p0[0]);
    double p2_1_p0_1 = (p2[1] - p0[1]);
    double p2_2_p0_2 = (p2[2] - p0[2]);
    double p2_3_p0_3 = (p2[3] - p0[3]);
    double p2_4_p0_4 = (p2[4] - p0[4]);
    double p2_5_p0_5 = (p2[5] - p0[5]);
    double p2_6_p0_6 = (p2[6] - p0[6]);
    double p2_7_p0_7 = (p2[7] - p0[7]);
    double l2;
    l2 = (1 * ((((((((p2_0_p0_0 * p2_0_p0_0) + (p2_1_p0_1 * p2_1_p0_1)) + (p2_2_p0_2 * p2_2_p0_2)) + (p2_3_p0_3 * p2_3_p0_3)) + (p2_4_p0_4 * p2_4_p0_4)) + (p2_5_p0_5 * p2_5_p0_5)) + (p2_6_p0_6 * p2_6_p0_6)) + (p2_7_p0_7 * p2_7_p0_7)));
    double p3_0_p0_0 = (p3[0] - p0[0]);
    double p3_1_p0_1 = (p3[1] - p0[1]);
    double p3_2_p0_2 = (p3[2] - p0[2]);
    double p3_3_p0_3 = (p3[3] - p0[3]);
    double p3_4_p0_4 = (p3[4] - p0[4]);
    double p3_5_p0_5 = (p3[5] - p0[5]);
    double p3_6_p0_6 = (p3[6] - p0[6]);
    double p3_7_p0_7 = (p3[7] - p0[7]);
    double l3;
    l3 = (1 * ((((((((p3_0_p0_0 * p3_0_p0_0) + (p3_1_p0_1 * p3_1_p0_1)) + (p3_2_p0_2 * p3_2_p0_2)) + (p3_3_p0_3 * p3_3_p0_3)) + (p3_4_p0_4 * p3_4_p0_4)) + (p3_5_p0_5 * p3_5_p0_5)) + (p3_6_p0_6 * p3_6_p0_6)) + (p3_7_p0_7 * p3_7_p0_7)));
    double q0_0_p0_0 = (q0[0] - p0[0]);
    double q0_1_p0_1 = (q0[1] - p0[1]);
    double q0_2_p0_2 = (q0[2] - p0[2]);
    double q0_3_p0_3 = (q0[3] - p0[3]);
    double q0_4_p0_4 = (q0[4] - p0[4]);
    double q0_5_p0_5 = (q0[5] - p0[5]);
    double q0_6_p0_6 = (q0[6] - p0[6]);
    double q0_7_p0_7 = (q0[7] - p0[7]);
    double a10;
    a10 = (2 * ((((((((p1_0_p0_0 * q0_0_p0_0) + (p1_1_p0_1 * q0_1_p0_1)) + (p1_2_p0_2 * q0_2_p0_2)) + (p1_3_p0_3 * q0_3_p0_3)) + (p1_4_p0_4 * q0_4_p0_4)) + (p1_5_p0_5 * q0_5_p0_5)) + (p1_6_p0_6 * q0_6_p0_6)) + (p1_7_p0_7 * q0_7_p0_7)));
    double q1_0_p0_0 = (q1[0] - p0[0]);
    double q1_1_p0_1 = (q1[1] - p0[1]);
    double q1_2_p0_2 = (q1[2] - p0[2]);
    double q1_3_p0_3 = (q1[3] - p0[3]);
    double q1_4_p0_4 = (q1[4] - p0[4]);
    double q1_5_p0_5 = (q1[5] - p0[5]);
    double q1_6_p0_6 = (q1[6] - p0[6]);
    double q1_7_p0_7 = (q1[7] - p0[7]);
    double a11;
    a11 = (2 * ((((((((p1_0_p0_0 * q1_0_p0_0) + (p1_1_p0_1 * q1_1_p0_1)) + (p1_2_p0_2 * q1_2_p0_2)) + (p1_3_p0_3 * q1_3_p0_3)) + (p1_4_p0_4 * q1_4_p0_4)) + (p1_5_p0_5 * q1_5_p0_5)) + (p1_6_p0_6 * q1_6_p0_6)) + (p1_7_p0_7 * q1_7_p0_7)));
    double q2_0_p0_0 = (q2[0] - p0[0]);
    double q2_1_p0_1 = (q2[1] - p0[1]);
    double q2_2_p0_2 = (q2[2] - p0[2]);
    double q2_3_p0_3 = (q2[3] - p0[3]);
    double q2_4_p0_4 = (q2[4] - p0[4]);
    double q2_5_p0_5 = (q2[5] - p0[5]);
    double q2_6_p0_6 = (q2[6] - p0[6]);
    double q2_7_p0_7 = (q2[7] - p0[7]);
    double a12;
    a12 = (2 * ((((((((p1_0_p0_0 * q2_0_p0_0) + (p1_1_p0_1 * q2_1_p0_1)) + (p1_2_p0_2 * q2_2_p0_2)) + (p1_3_p0_3 * q2_3_p0_3)) + (p1_4_p0_4 * q2_4_p0_4)) + (p1_5_p0_5 * q2_5_p0_5)) + (p1_6_p0_6 * q2_6_p0_6)) + (p1_7_p0_7 * q2_7_p0_7)));
    double a20;
    a20 = (2 * ((((((((p2_0_p0_0 * q0_0_p0_0) + (p2_1_p0_1 * q0_1_p0_1)) + (p2_2_p0_2 * q0_2_p0_2)) + (p2_3_p0_3 * q0_3_p0_3)) + (p2_4_p0_4 * q0_4_p0_4)) + (p2_5_p0_5 * q0_5_p0_5)) + (p2_6_p0_6 * q0_6_p0_6)) + (p2_7_p0_7 * q0_7_p0_7)));
    double a21;
    a21 = (2 * ((((((((p2_0_p0_0 * q1_0_p0_0) + (p2_1_p0_1 * q1_1_p0_1)) + (p2_2_p0_2 * q1_2_p0_2)) + (p2_3_p0_3 * q1_3_p0_3)) + (p2_4_p0_4 * q1_4_p0_4)) + (p2_5_p0_5 * q1_5_p0_5)) + (p2_6_p0_6 * q1_6_p0_6)) + (p2_7_p0_7 * q1_7_p0_7)));
    double a22;
    a22 = (2 * ((((((((p2_0_p0_0 * q2_0_p0_0) + (p2_1_p0_1 * q2_1_p0_1)) + (p2_2_p0_2 * q2_2_p0_2)) + (p2_3_p0_3 * q2_3_p0_3)) + (p2_4_p0_4 * q2_4_p0_4)) + (p2_5_p0_5 * q2_5_p0_5)) + (p2_6_p0_6 * q2_6_p0_6)) + (p2_7_p0_7 * q2_7_p0_7)));
    double a30;
    a30 = (2 * ((((((((p3_0_p0_0 * q0_0_p0_0) + (p3_1_p0_1 * q0_1_p0_1)) + (p3_2_p0_2 * q0_2_p0_2)) + (p3_3_p0_3 * q0_3_p0_3)) + (p3_4_p0_4 * q0_4_p0_4)) + (p3_5_p0_5 * q0_5_p0_5)) + (p3_6_p0_6 * q0_6_p0_6)) + (p3_7_p0_7 * q0_7_p0_7)));
    double a31;
    a31 = (2 * ((((((((p3_0_p0_0 * q1_0_p0_0) + (p3_1_p0_1 * q1_1_p0_1)) + (p3_2_p0_2 * q1_2_p0_2)) + (p3_3_p0_3 * q1_3_p0_3)) + (p3_4_p0_4 * q1_4_p0_4)) + (p3_5_p0_5 * q1_5_p0_5)) + (p3_6_p0_6 * q1_6_p0_6)) + (p3_7_p0_7 * q1_7_p0_7)));
    double a32;
    a32 = (2 * ((((((((p3_0_p0_0 * q2_0_p0_0) + (p3_1_p0_1 * q2_1_p0_1)) + (p3_2_p0_2 * q2_2_p0_2)) + (p3_3_p0_3 * q2_3_p0_3)) + (p3_4_p0_4 * q2_4_p0_4)) + (p3_5_p0_5 * q2_5_p0_5)) + (p3_6_p0_6 * q2_6_p0_6)) + (p3_7_p0_7 * q2_7_p0_7)));
    double b00;
    b00 = ((a11 * a22) - (a12 * a21));
    double b01;
    b01 = -(a22 - a21);
    double b02;
    b02 = (a12 - a11);
    double b10;
    b10 = -((a10 * a22) - (a12 * a20));
    double b11;
    b11 = (a22 - a20);
    double b12;
    b12 = -(a12 - a10);
    double b20;
    b20 = ((a10 * a21) - (a11 * a20));
    double b21;
    b21 = -(a21 - a20);
    double b22;
    b22 = (a11 - a10);
    double Delta;
    Delta = ((b00 + b10) + b20);
    double DeltaLambda0;
    DeltaLambda0 = ((b00 + (b01 * l1)) + (b02 * l2));
    double DeltaLambda1;
    DeltaLambda1 = ((b10 + (b11 * l1)) + (b12 * l2));
    double DeltaLambda2;
    DeltaLambda2 = ((b20 + (b21 * l1)) + (b22 * l2));
    double r;
    r = ((Delta * l3) - (((a30 * DeltaLambda0) + (a31 * DeltaLambda1)) + (a32 * DeltaLambda2)));
    double eps;
    double max1 = fabs(p1_0_p0_0);
    if( (max1 < fabs(p1_3_p0_3)) )
    {
        max1 = fabs(p1_3_p0_3);
    }
    if( (max1 < fabs(p1_1_p0_1)) )
    {
        max1 = fabs(p1_1_p0_1);
    }
    if( (max1 < fabs(p1_4_p0_4)) )
    {
        max1 = fabs(p1_4_p0_4);
    }
    if( (max1 < fabs(p1_5_p0_5)) )
    {
        max1 = fabs(p1_5_p0_5);
    }
    if( (max1 < fabs(p1_6_p0_6)) )
    {
        max1 = fabs(p1_6_p0_6);
    }
    if( (max1 < fabs(p1_2_p0_2)) )
    {
        max1 = fabs(p1_2_p0_2);
    }
    if( (max1 < fabs(p1_7_p0_7)) )
    {
        max1 = fabs(p1_7_p0_7);
    }
    double max2 = fabs(p2_0_p0_0);
    if( (max2 < fabs(p2_1_p0_1)) )
    {
        max2 = fabs(p2_1_p0_1);
    }
    if( (max2 < fabs(p2_2_p0_2)) )
    {
        max2 = fabs(p2_2_p0_2);
    }
    if( (max2 < fabs(p2_3_p0_3)) )
    {
        max2 = fabs(p2_3_p0_3);
    }
    if( (max2 < fabs(p2_4_p0_4)) )
    {
        max2 = fabs(p2_4_p0_4);
    }
    if( (max2 < fabs(p2_5_p0_5)) )
    {
        max2 = fabs(p2_5_p0_5);
    }
    if( (max2 < fabs(p2_6_p0_6)) )
    {
        max2 = fabs(p2_6_p0_6);
    }
    if( (max2 < fabs(p2_7_p0_7)) )
    {
        max2 = fabs(p2_7_p0_7);
    }
    double max3 = fabs(q0_0_p0_0);
    if( (max3 < fabs(q0_1_p0_1)) )
    {
        max3 = fabs(q0_1_p0_1);
    }
    if( (max3 < fabs(q0_2_p0_2)) )
    {
        max3 = fabs(q0_2_p0_2);
    }
    if( (max3 < fabs(q0_3_p0_3)) )
    {
        max3 = fabs(q0_3_p0_3);
    }
    if( (max3 < fabs(q0_4_p0_4)) )
    {
        max3 = fabs(q0_4_p0_4);
    }
    if( (max3 < fabs(q0_5_p0_5)) )
    {
        max3 = fabs(q0_5_p0_5);
    }
    if( (max3 < fabs(q0_6_p0_6)) )
    {
        max3 = fabs(q0_6_p0_6);
    }
    if( (max3 < fabs(q0_7_p0_7)) )
    {
        max3 = fabs(q0_7_p0_7);
    }
    if( (max3 < fabs(q1_0_p0_0)) )
    {
        max3 = fabs(q1_0_p0_0);
    }
    if( (max3 < fabs(q1_1_p0_1)) )
    {
        max3 = fabs(q1_1_p0_1);
    }
    if( (max3 < fabs(q1_2_p0_2)) )
    {
        max3 = fabs(q1_2_p0_2);
    }
    if( (max3 < fabs(q1_3_p0_3)) )
    {
        max3 = fabs(q1_3_p0_3);
    }
    if( (max3 < fabs(q1_4_p0_4)) )
    {
        max3 = fabs(q1_4_p0_4);
    }
    if( (max3 < fabs(q1_5_p0_5)) )
    {
        max3 = fabs(q1_5_p0_5);
    }
    if( (max3 < fabs(q1_6_p0_6)) )
    {
        max3 = fabs(q1_6_p0_6);
    }
    if( (max3 < fabs(q1_7_p0_7)) )
    {
        max3 = fabs(q1_7_p0_7);
    }
    double max4 = fabs(q1_0_p0_0);
    if( (max4 < fabs(q1_1_p0_1)) )
    {
        max4 = fabs(q1_1_p0_1);
    }
    if( (max4 < fabs(q1_2_p0_2)) )
    {
        max4 = fabs(q1_2_p0_2);
    }
    if( (max4 < fabs(q1_3_p0_3)) )
    {
        max4 = fabs(q1_3_p0_3);
    }
    if( (max4 < fabs(q1_4_p0_4)) )
    {
        max4 = fabs(q1_4_p0_4);
    }
    if( (max4 < fabs(q1_5_p0_5)) )
    {
        max4 = fabs(q1_5_p0_5);
    }
    if( (max4 < fabs(q1_6_p0_6)) )
    {
        max4 = fabs(q1_6_p0_6);
    }
    if( (max4 < fabs(q1_7_p0_7)) )
    {
        max4 = fabs(q1_7_p0_7);
    }
    if( (max4 < fabs(q2_0_p0_0)) )
    {
        max4 = fabs(q2_0_p0_0);
    }
    if( (max4 < fabs(q2_1_p0_1)) )
    {
        max4 = fabs(q2_1_p0_1);
    }
    if( (max4 < fabs(q2_2_p0_2)) )
    {
        max4 = fabs(q2_2_p0_2);
    }
    if( (max4 < fabs(q2_3_p0_3)) )
    {
        max4 = fabs(q2_3_p0_3);
    }
    if( (max4 < fabs(q2_4_p0_4)) )
    {
        max4 = fabs(q2_4_p0_4);
    }
    if( (max4 < fabs(q2_5_p0_5)) )
    {
        max4 = fabs(q2_5_p0_5);
    }
    if( (max4 < fabs(q2_6_p0_6)) )
    {
        max4 = fabs(q2_6_p0_6);
    }
    if( (max4 < fabs(q2_7_p0_7)) )
    {
        max4 = fabs(q2_7_p0_7);
    }
    double lower_bound_1;
    double upper_bound_1;
    int Delta_sign;
    int int_tmp_result;
    lower_bound_1 = max2;
    upper_bound_1 = max2;
    if( (max1 < lower_bound_1) )
    {
        lower_bound_1 = max1;
    }
    else
    {
        if( (max1 > upper_bound_1) )
        {
            upper_bound_1 = max1;
        }
    }
    if( (max4 < lower_bound_1) )
    {
        lower_bound_1 = max4;
    }
    else
    {
        if( (max4 > upper_bound_1) )
        {
            upper_bound_1 = max4;
        }
    }
    if( (max3 < lower_bound_1) )
    {
        lower_bound_1 = max3;
    }
    else
    {
        if( (max3 > upper_bound_1) )
        {
            upper_bound_1 = max3;
        }
    }
    if( (lower_bound_1 < 1.26419510663115923609e-74) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 1.29807421463370647479e+33) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (8.71140112255785451890e-13 * (((max1 * max3) * max2) * max4));
        if( (Delta > eps) )
        {
            int_tmp_result = 1;
        }
        else
        {
            if( (Delta < -eps) )
            {
                int_tmp_result = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    Delta_sign = int_tmp_result;
    int int_tmp_result_3SPBwDj;
    double max5 = max2;
    if( (max5 < max1) )
    {
        max5 = max1;
    }
    if( (max5 < max4) )
    {
        max5 = max4;
    }
    if( (max5 < fabs(p3_0_p0_0)) )
    {
        max5 = fabs(p3_0_p0_0);
    }
    if( (max5 < fabs(p3_1_p0_1)) )
    {
        max5 = fabs(p3_1_p0_1);
    }
    if( (max5 < fabs(p3_2_p0_2)) )
    {
        max5 = fabs(p3_2_p0_2);
    }
    if( (max5 < fabs(p3_3_p0_3)) )
    {
        max5 = fabs(p3_3_p0_3);
    }
    if( (max5 < fabs(p3_4_p0_4)) )
    {
        max5 = fabs(p3_4_p0_4);
    }
    if( (max5 < fabs(p3_5_p0_5)) )
    {
        max5 = fabs(p3_5_p0_5);
    }
    if( (max5 < fabs(p3_6_p0_6)) )
    {
        max5 = fabs(p3_6_p0_6);
    }
    if( (max5 < fabs(p3_7_p0_7)) )
    {
        max5 = fabs(p3_7_p0_7);
    }
    double max6 = max2;
    if( (max6 < max4) )
    {
        max6 = max4;
    }
    if( (max6 < max3) )
    {
        max6 = max3;
    }
    double max7 = max2;
    if( (max7 < max4) )
    {
        max7 = max4;
    }
    if( (max7 < max3) )
    {
        max7 = max3;
    }
    if( (max6 < max7) )
    {
        max6 = max7;
    }
    double max8 = max3;
    if( (max8 < fabs(p3_0_p0_0)) )
    {
        max8 = fabs(p3_0_p0_0);
    }
    if( (max8 < fabs(p3_1_p0_1)) )
    {
        max8 = fabs(p3_1_p0_1);
    }
    if( (max8 < fabs(p3_2_p0_2)) )
    {
        max8 = fabs(p3_2_p0_2);
    }
    if( (max8 < fabs(p3_3_p0_3)) )
    {
        max8 = fabs(p3_3_p0_3);
    }
    if( (max8 < fabs(p3_4_p0_4)) )
    {
        max8 = fabs(p3_4_p0_4);
    }
    if( (max8 < fabs(p3_5_p0_5)) )
    {
        max8 = fabs(p3_5_p0_5);
    }
    if( (max8 < fabs(p3_6_p0_6)) )
    {
        max8 = fabs(p3_6_p0_6);
    }
    if( (max8 < fabs(p3_7_p0_7)) )
    {
        max8 = fabs(p3_7_p0_7);
    }
    if( (max6 < max8) )
    {
        max6 = max8;
    }
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if( (max5 < lower_bound_1) )
    {
        lower_bound_1 = max5;
    }
    else
    {
        if( (max5 > upper_bound_1) )
        {
            upper_bound_1 = max5;
        }
    }
    if( (max6 < lower_bound_1) )
    {
        lower_bound_1 = max6;
    }
    else
    {
        if( (max6 > upper_bound_1) )
        {
            upper_bound_1 = max6;
        }
    }
    if( (max7 < lower_bound_1) )
    {
        lower_bound_1 = max7;
    }
    if( (max8 < lower_bound_1) )
    {
        lower_bound_1 = max8;
    }
    if( (lower_bound_1 < 2.82528483194754087282e-50) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 1.29807421463370647479e+33) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (4.37492894694731169807e-11 * (((((max1 * max8) * max7) * max7) * max6) * max5));
        if( (r > eps) )
        {
            int_tmp_result_3SPBwDj = 1;
        }
        else
        {
            if( (r < -eps) )
            {
                int_tmp_result_3SPBwDj = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    return (Delta_sign * int_tmp_result_3SPBwDj);
}


int avro_side4_8d_filter( const double* p0, const double* p1, const double* p2, const double* p3, const double* p4, const double* q0, const double* q1, const double* q2, const double* q3) {
    double p1_0_p0_0 = (p1[0] - p0[0]);
    double p1_1_p0_1 = (p1[1] - p0[1]);
    double p1_2_p0_2 = (p1[2] - p0[2]);
    double p1_3_p0_3 = (p1[3] - p0[3]);
    double p1_4_p0_4 = (p1[4] - p0[4]);
    double p1_5_p0_5 = (p1[5] - p0[5]);
    double p1_6_p0_6 = (p1[6] - p0[6]);
    double p1_7_p0_7 = (p1[7] - p0[7]);
    double l1;
    l1 = (1 * ((((((((p1_0_p0_0 * p1_0_p0_0) + (p1_1_p0_1 * p1_1_p0_1)) + (p1_2_p0_2 * p1_2_p0_2)) + (p1_3_p0_3 * p1_3_p0_3)) + (p1_4_p0_4 * p1_4_p0_4)) + (p1_5_p0_5 * p1_5_p0_5)) + (p1_6_p0_6 * p1_6_p0_6)) + (p1_7_p0_7 * p1_7_p0_7)));
    double p2_0_p0_0 = (p2[0] - p0[0]);
    double p2_1_p0_1 = (p2[1] - p0[1]);
    double p2_2_p0_2 = (p2[2] - p0[2]);
    double p2_3_p0_3 = (p2[3] - p0[3]);
    double p2_4_p0_4 = (p2[4] - p0[4]);
    double p2_5_p0_5 = (p2[5] - p0[5]);
    double p2_6_p0_6 = (p2[6] - p0[6]);
    double p2_7_p0_7 = (p2[7] - p0[7]);
    double l2;
    l2 = (1 * ((((((((p2_0_p0_0 * p2_0_p0_0) + (p2_1_p0_1 * p2_1_p0_1)) + (p2_2_p0_2 * p2_2_p0_2)) + (p2_3_p0_3 * p2_3_p0_3)) + (p2_4_p0_4 * p2_4_p0_4)) + (p2_5_p0_5 * p2_5_p0_5)) + (p2_6_p0_6 * p2_6_p0_6)) + (p2_7_p0_7 * p2_7_p0_7)));
    double p3_0_p0_0 = (p3[0] - p0[0]);
    double p3_1_p0_1 = (p3[1] - p0[1]);
    double p3_2_p0_2 = (p3[2] - p0[2]);
    double p3_3_p0_3 = (p3[3] - p0[3]);
    double p3_4_p0_4 = (p3[4] - p0[4]);
    double p3_5_p0_5 = (p3[5] - p0[5]);
    double p3_6_p0_6 = (p3[6] - p0[6]);
    double p3_7_p0_7 = (p3[7] - p0[7]);
    double l3;
    l3 = (1 * ((((((((p3_0_p0_0 * p3_0_p0_0) + (p3_1_p0_1 * p3_1_p0_1)) + (p3_2_p0_2 * p3_2_p0_2)) + (p3_3_p0_3 * p3_3_p0_3)) + (p3_4_p0_4 * p3_4_p0_4)) + (p3_5_p0_5 * p3_5_p0_5)) + (p3_6_p0_6 * p3_6_p0_6)) + (p3_7_p0_7 * p3_7_p0_7)));
    double p4_0_p0_0 = (p4[0] - p0[0]);
    double p4_1_p0_1 = (p4[1] - p0[1]);
    double p4_2_p0_2 = (p4[2] - p0[2]);
    double p4_3_p0_3 = (p4[3] - p0[3]);
    double p4_4_p0_4 = (p4[4] - p0[4]);
    double p4_5_p0_5 = (p4[5] - p0[5]);
    double p4_6_p0_6 = (p4[6] - p0[6]);
    double p4_7_p0_7 = (p4[7] - p0[7]);
    double l4;
    l4 = (1 * ((((((((p4_0_p0_0 * p4_0_p0_0) + (p4_1_p0_1 * p4_1_p0_1)) + (p4_2_p0_2 * p4_2_p0_2)) + (p4_3_p0_3 * p4_3_p0_3)) + (p4_4_p0_4 * p4_4_p0_4)) + (p4_5_p0_5 * p4_5_p0_5)) + (p4_6_p0_6 * p4_6_p0_6)) + (p4_7_p0_7 * p4_7_p0_7)));
    double q0_0_p0_0 = (q0[0] - p0[0]);
    double q0_1_p0_1 = (q0[1] - p0[1]);
    double q0_2_p0_2 = (q0[2] - p0[2]);
    double q0_3_p0_3 = (q0[3] - p0[3]);
    double q0_4_p0_4 = (q0[4] - p0[4]);
    double q0_5_p0_5 = (q0[5] - p0[5]);
    double q0_6_p0_6 = (q0[6] - p0[6]);
    double q0_7_p0_7 = (q0[7] - p0[7]);
    double a10;
    a10 = (2 * ((((((((p1_0_p0_0 * q0_0_p0_0) + (p1_1_p0_1 * q0_1_p0_1)) + (p1_2_p0_2 * q0_2_p0_2)) + (p1_3_p0_3 * q0_3_p0_3)) + (p1_4_p0_4 * q0_4_p0_4)) + (p1_5_p0_5 * q0_5_p0_5)) + (p1_6_p0_6 * q0_6_p0_6)) + (p1_7_p0_7 * q0_7_p0_7)));
    double q1_0_p0_0 = (q1[0] - p0[0]);
    double q1_1_p0_1 = (q1[1] - p0[1]);
    double q1_2_p0_2 = (q1[2] - p0[2]);
    double q1_3_p0_3 = (q1[3] - p0[3]);
    double q1_4_p0_4 = (q1[4] - p0[4]);
    double q1_5_p0_5 = (q1[5] - p0[5]);
    double q1_6_p0_6 = (q1[6] - p0[6]);
    double q1_7_p0_7 = (q1[7] - p0[7]);
    double a11;
    a11 = (2 * ((((((((p1_0_p0_0 * q1_0_p0_0) + (p1_1_p0_1 * q1_1_p0_1)) + (p1_2_p0_2 * q1_2_p0_2)) + (p1_3_p0_3 * q1_3_p0_3)) + (p1_4_p0_4 * q1_4_p0_4)) + (p1_5_p0_5 * q1_5_p0_5)) + (p1_6_p0_6 * q1_6_p0_6)) + (p1_7_p0_7 * q1_7_p0_7)));
    double q2_0_p0_0 = (q2[0] - p0[0]);
    double q2_1_p0_1 = (q2[1] - p0[1]);
    double q2_2_p0_2 = (q2[2] - p0[2]);
    double q2_3_p0_3 = (q2[3] - p0[3]);
    double q2_4_p0_4 = (q2[4] - p0[4]);
    double q2_5_p0_5 = (q2[5] - p0[5]);
    double q2_6_p0_6 = (q2[6] - p0[6]);
    double q2_7_p0_7 = (q2[7] - p0[7]);
    double a12;
    a12 = (2 * ((((((((p1_0_p0_0 * q2_0_p0_0) + (p1_1_p0_1 * q2_1_p0_1)) + (p1_2_p0_2 * q2_2_p0_2)) + (p1_3_p0_3 * q2_3_p0_3)) + (p1_4_p0_4 * q2_4_p0_4)) + (p1_5_p0_5 * q2_5_p0_5)) + (p1_6_p0_6 * q2_6_p0_6)) + (p1_7_p0_7 * q2_7_p0_7)));
    double q3_0_p0_0 = (q3[0] - p0[0]);
    double q3_1_p0_1 = (q3[1] - p0[1]);
    double q3_2_p0_2 = (q3[2] - p0[2]);
    double q3_3_p0_3 = (q3[3] - p0[3]);
    double q3_4_p0_4 = (q3[4] - p0[4]);
    double q3_5_p0_5 = (q3[5] - p0[5]);
    double q3_6_p0_6 = (q3[6] - p0[6]);
    double q3_7_p0_7 = (q3[7] - p0[7]);
    double a13;
    a13 = (2 * ((((((((p1_0_p0_0 * q3_0_p0_0) + (p1_1_p0_1 * q3_1_p0_1)) + (p1_2_p0_2 * q3_2_p0_2)) + (p1_3_p0_3 * q3_3_p0_3)) + (p1_4_p0_4 * q3_4_p0_4)) + (p1_5_p0_5 * q3_5_p0_5)) + (p1_6_p0_6 * q3_6_p0_6)) + (p1_7_p0_7 * q3_7_p0_7)));
    double a20;
    a20 = (2 * ((((((((p2_0_p0_0 * q0_0_p0_0) + (p2_1_p0_1 * q0_1_p0_1)) + (p2_2_p0_2 * q0_2_p0_2)) + (p2_3_p0_3 * q0_3_p0_3)) + (p2_4_p0_4 * q0_4_p0_4)) + (p2_5_p0_5 * q0_5_p0_5)) + (p2_6_p0_6 * q0_6_p0_6)) + (p2_7_p0_7 * q0_7_p0_7)));
    double a21;
    a21 = (2 * ((((((((p2_0_p0_0 * q1_0_p0_0) + (p2_1_p0_1 * q1_1_p0_1)) + (p2_2_p0_2 * q1_2_p0_2)) + (p2_3_p0_3 * q1_3_p0_3)) + (p2_4_p0_4 * q1_4_p0_4)) + (p2_5_p0_5 * q1_5_p0_5)) + (p2_6_p0_6 * q1_6_p0_6)) + (p2_7_p0_7 * q1_7_p0_7)));
    double a22;
    a22 = (2 * ((((((((p2_0_p0_0 * q2_0_p0_0) + (p2_1_p0_1 * q2_1_p0_1)) + (p2_2_p0_2 * q2_2_p0_2)) + (p2_3_p0_3 * q2_3_p0_3)) + (p2_4_p0_4 * q2_4_p0_4)) + (p2_5_p0_5 * q2_5_p0_5)) + (p2_6_p0_6 * q2_6_p0_6)) + (p2_7_p0_7 * q2_7_p0_7)));
    double a23;
    a23 = (2 * ((((((((p2_0_p0_0 * q3_0_p0_0) + (p2_1_p0_1 * q3_1_p0_1)) + (p2_2_p0_2 * q3_2_p0_2)) + (p2_3_p0_3 * q3_3_p0_3)) + (p2_4_p0_4 * q3_4_p0_4)) + (p2_5_p0_5 * q3_5_p0_5)) + (p2_6_p0_6 * q3_6_p0_6)) + (p2_7_p0_7 * q3_7_p0_7)));
    double a30;
    a30 = (2 * ((((((((p3_0_p0_0 * q0_0_p0_0) + (p3_1_p0_1 * q0_1_p0_1)) + (p3_2_p0_2 * q0_2_p0_2)) + (p3_3_p0_3 * q0_3_p0_3)) + (p3_4_p0_4 * q0_4_p0_4)) + (p3_5_p0_5 * q0_5_p0_5)) + (p3_6_p0_6 * q0_6_p0_6)) + (p3_7_p0_7 * q0_7_p0_7)));
    double a31;
    a31 = (2 * ((((((((p3_0_p0_0 * q1_0_p0_0) + (p3_1_p0_1 * q1_1_p0_1)) + (p3_2_p0_2 * q1_2_p0_2)) + (p3_3_p0_3 * q1_3_p0_3)) + (p3_4_p0_4 * q1_4_p0_4)) + (p3_5_p0_5 * q1_5_p0_5)) + (p3_6_p0_6 * q1_6_p0_6)) + (p3_7_p0_7 * q1_7_p0_7)));
    double a32;
    a32 = (2 * ((((((((p3_0_p0_0 * q2_0_p0_0) + (p3_1_p0_1 * q2_1_p0_1)) + (p3_2_p0_2 * q2_2_p0_2)) + (p3_3_p0_3 * q2_3_p0_3)) + (p3_4_p0_4 * q2_4_p0_4)) + (p3_5_p0_5 * q2_5_p0_5)) + (p3_6_p0_6 * q2_6_p0_6)) + (p3_7_p0_7 * q2_7_p0_7)));
    double a33;
    a33 = (2 * ((((((((p3_0_p0_0 * q3_0_p0_0) + (p3_1_p0_1 * q3_1_p0_1)) + (p3_2_p0_2 * q3_2_p0_2)) + (p3_3_p0_3 * q3_3_p0_3)) + (p3_4_p0_4 * q3_4_p0_4)) + (p3_5_p0_5 * q3_5_p0_5)) + (p3_6_p0_6 * q3_6_p0_6)) + (p3_7_p0_7 * q3_7_p0_7)));
    double a40;
    a40 = (2 * ((((((((p4_0_p0_0 * q0_0_p0_0) + (p4_1_p0_1 * q0_1_p0_1)) + (p4_2_p0_2 * q0_2_p0_2)) + (p4_3_p0_3 * q0_3_p0_3)) + (p4_4_p0_4 * q0_4_p0_4)) + (p4_5_p0_5 * q0_5_p0_5)) + (p4_6_p0_6 * q0_6_p0_6)) + (p4_7_p0_7 * q0_7_p0_7)));
    double a41;
    a41 = (2 * ((((((((p4_0_p0_0 * q1_0_p0_0) + (p4_1_p0_1 * q1_1_p0_1)) + (p4_2_p0_2 * q1_2_p0_2)) + (p4_3_p0_3 * q1_3_p0_3)) + (p4_4_p0_4 * q1_4_p0_4)) + (p4_5_p0_5 * q1_5_p0_5)) + (p4_6_p0_6 * q1_6_p0_6)) + (p4_7_p0_7 * q1_7_p0_7)));
    double a42;
    a42 = (2 * ((((((((p4_0_p0_0 * q2_0_p0_0) + (p4_1_p0_1 * q2_1_p0_1)) + (p4_2_p0_2 * q2_2_p0_2)) + (p4_3_p0_3 * q2_3_p0_3)) + (p4_4_p0_4 * q2_4_p0_4)) + (p4_5_p0_5 * q2_5_p0_5)) + (p4_6_p0_6 * q2_6_p0_6)) + (p4_7_p0_7 * q2_7_p0_7)));
    double a43;
    a43 = (2 * ((((((((p4_0_p0_0 * q3_0_p0_0) + (p4_1_p0_1 * q3_1_p0_1)) + (p4_2_p0_2 * q3_2_p0_2)) + (p4_3_p0_3 * q3_3_p0_3)) + (p4_4_p0_4 * q3_4_p0_4)) + (p4_5_p0_5 * q3_5_p0_5)) + (p4_6_p0_6 * q3_6_p0_6)) + (p4_7_p0_7 * q3_7_p0_7)));
    double b00;
    b00 = (((((((a11 * a22) * a33) - ((a11 * a23) * a32)) - ((a12 * a21) * a33)) + ((a12 * a23) * a31)) + ((a13 * a21) * a32)) - ((a13 * a22) * a31));
    double b01;
    b01 = -((((((a21 * a32) - (a22 * a31)) - (a21 * a33)) + (a23 * a31)) + (a22 * a33)) - (a23 * a32));
    double b02;
    b02 = ((((((a11 * a32) - (a12 * a31)) - (a11 * a33)) + (a13 * a31)) + (a12 * a33)) - (a13 * a32));
    double b03;
    b03 = -((((((a11 * a22) - (a12 * a21)) - (a11 * a23)) + (a13 * a21)) + (a12 * a23)) - (a13 * a22));
    double b10;
    b10 = -(((((((a10 * a22) * a33) - ((a10 * a23) * a32)) - ((a12 * a20) * a33)) + ((a12 * a23) * a30)) + ((a13 * a20) * a32)) - ((a13 * a22) * a30));
    double b11;
    b11 = ((((((a20 * a32) - (a22 * a30)) - (a20 * a33)) + (a23 * a30)) + (a22 * a33)) - (a23 * a32));
    double b12;
    b12 = -((((((a10 * a32) - (a12 * a30)) - (a10 * a33)) + (a13 * a30)) + (a12 * a33)) - (a13 * a32));
    double b13;
    b13 = ((((((a10 * a22) - (a12 * a20)) - (a10 * a23)) + (a13 * a20)) + (a12 * a23)) - (a13 * a22));
    double b20;
    b20 = (((((((a10 * a21) * a33) - ((a10 * a23) * a31)) - ((a11 * a20) * a33)) + ((a11 * a23) * a30)) + ((a13 * a20) * a31)) - ((a13 * a21) * a30));
    double b21;
    b21 = -((((((a20 * a31) - (a21 * a30)) - (a20 * a33)) + (a23 * a30)) + (a21 * a33)) - (a23 * a31));
    double b22;
    b22 = ((((((a10 * a31) - (a11 * a30)) - (a10 * a33)) + (a13 * a30)) + (a11 * a33)) - (a13 * a31));
    double b23;
    b23 = -((((((a10 * a21) - (a11 * a20)) - (a10 * a23)) + (a13 * a20)) + (a11 * a23)) - (a13 * a21));
    double b30;
    b30 = -(((((((a10 * a21) * a32) - ((a10 * a22) * a31)) - ((a11 * a20) * a32)) + ((a11 * a22) * a30)) + ((a12 * a20) * a31)) - ((a12 * a21) * a30));
    double b31;
    b31 = ((((((a20 * a31) - (a21 * a30)) - (a20 * a32)) + (a22 * a30)) + (a21 * a32)) - (a22 * a31));
    double b32;
    b32 = -((((((a10 * a31) - (a11 * a30)) - (a10 * a32)) + (a12 * a30)) + (a11 * a32)) - (a12 * a31));
    double b33;
    b33 = ((((((a10 * a21) - (a11 * a20)) - (a10 * a22)) + (a12 * a20)) + (a11 * a22)) - (a12 * a21));
    double Delta;
    Delta = (((b00 + b10) + b20) + b30);
    double DeltaLambda0;
    DeltaLambda0 = (((b00 + (b01 * l1)) + (b02 * l2)) + (b03 * l3));
    double DeltaLambda1;
    DeltaLambda1 = (((b10 + (b11 * l1)) + (b12 * l2)) + (b13 * l3));
    double DeltaLambda2;
    DeltaLambda2 = (((b20 + (b21 * l1)) + (b22 * l2)) + (b23 * l3));
    double DeltaLambda3;
    DeltaLambda3 = (((b30 + (b31 * l1)) + (b32 * l2)) + (b33 * l3));
    double r;
    r = ((Delta * l4) - ((((a40 * DeltaLambda0) + (a41 * DeltaLambda1)) + (a42 * DeltaLambda2)) + (a43 * DeltaLambda3)));
    double eps;
    double max1 = fabs(p1_0_p0_0);
    if( (max1 < fabs(p1_2_p0_2)) )
    {
        max1 = fabs(p1_2_p0_2);
    }
    if( (max1 < fabs(p1_7_p0_7)) )
    {
        max1 = fabs(p1_7_p0_7);
    }
    if( (max1 < fabs(p1_3_p0_3)) )
    {
        max1 = fabs(p1_3_p0_3);
    }
    if( (max1 < fabs(p1_1_p0_1)) )
    {
        max1 = fabs(p1_1_p0_1);
    }
    if( (max1 < fabs(p1_4_p0_4)) )
    {
        max1 = fabs(p1_4_p0_4);
    }
    if( (max1 < fabs(p1_5_p0_5)) )
    {
        max1 = fabs(p1_5_p0_5);
    }
    if( (max1 < fabs(p1_6_p0_6)) )
    {
        max1 = fabs(p1_6_p0_6);
    }
    double max2 = fabs(p2_0_p0_0);
    if( (max2 < fabs(p2_3_p0_3)) )
    {
        max2 = fabs(p2_3_p0_3);
    }
    if( (max2 < fabs(p2_1_p0_1)) )
    {
        max2 = fabs(p2_1_p0_1);
    }
    if( (max2 < fabs(p2_5_p0_5)) )
    {
        max2 = fabs(p2_5_p0_5);
    }
    if( (max2 < fabs(p2_6_p0_6)) )
    {
        max2 = fabs(p2_6_p0_6);
    }
    if( (max2 < fabs(p2_7_p0_7)) )
    {
        max2 = fabs(p2_7_p0_7);
    }
    if( (max2 < fabs(p2_2_p0_2)) )
    {
        max2 = fabs(p2_2_p0_2);
    }
    if( (max2 < fabs(p2_4_p0_4)) )
    {
        max2 = fabs(p2_4_p0_4);
    }
    double max3 = fabs(p3_5_p0_5);
    if( (max3 < fabs(p3_0_p0_0)) )
    {
        max3 = fabs(p3_0_p0_0);
    }
    if( (max3 < fabs(p3_1_p0_1)) )
    {
        max3 = fabs(p3_1_p0_1);
    }
    if( (max3 < fabs(p3_2_p0_2)) )
    {
        max3 = fabs(p3_2_p0_2);
    }
    if( (max3 < fabs(p3_3_p0_3)) )
    {
        max3 = fabs(p3_3_p0_3);
    }
    if( (max3 < fabs(p3_4_p0_4)) )
    {
        max3 = fabs(p3_4_p0_4);
    }
    if( (max3 < fabs(p3_6_p0_6)) )
    {
        max3 = fabs(p3_6_p0_6);
    }
    if( (max3 < fabs(p3_7_p0_7)) )
    {
        max3 = fabs(p3_7_p0_7);
    }
    double max4 = fabs(q0_0_p0_0);
    if( (max4 < fabs(q0_1_p0_1)) )
    {
        max4 = fabs(q0_1_p0_1);
    }
    if( (max4 < fabs(q0_2_p0_2)) )
    {
        max4 = fabs(q0_2_p0_2);
    }
    if( (max4 < fabs(q0_3_p0_3)) )
    {
        max4 = fabs(q0_3_p0_3);
    }
    if( (max4 < fabs(q0_4_p0_4)) )
    {
        max4 = fabs(q0_4_p0_4);
    }
    if( (max4 < fabs(q0_5_p0_5)) )
    {
        max4 = fabs(q0_5_p0_5);
    }
    if( (max4 < fabs(q0_6_p0_6)) )
    {
        max4 = fabs(q0_6_p0_6);
    }
    if( (max4 < fabs(q0_7_p0_7)) )
    {
        max4 = fabs(q0_7_p0_7);
    }
    if( (max4 < fabs(q1_0_p0_0)) )
    {
        max4 = fabs(q1_0_p0_0);
    }
    if( (max4 < fabs(q1_1_p0_1)) )
    {
        max4 = fabs(q1_1_p0_1);
    }
    if( (max4 < fabs(q1_2_p0_2)) )
    {
        max4 = fabs(q1_2_p0_2);
    }
    if( (max4 < fabs(q1_3_p0_3)) )
    {
        max4 = fabs(q1_3_p0_3);
    }
    if( (max4 < fabs(q1_4_p0_4)) )
    {
        max4 = fabs(q1_4_p0_4);
    }
    if( (max4 < fabs(q1_5_p0_5)) )
    {
        max4 = fabs(q1_5_p0_5);
    }
    if( (max4 < fabs(q1_6_p0_6)) )
    {
        max4 = fabs(q1_6_p0_6);
    }
    if( (max4 < fabs(q1_7_p0_7)) )
    {
        max4 = fabs(q1_7_p0_7);
    }
    double max5 = fabs(q1_0_p0_0);
    if( (max5 < fabs(q1_1_p0_1)) )
    {
        max5 = fabs(q1_1_p0_1);
    }
    if( (max5 < fabs(q1_2_p0_2)) )
    {
        max5 = fabs(q1_2_p0_2);
    }
    if( (max5 < fabs(q1_3_p0_3)) )
    {
        max5 = fabs(q1_3_p0_3);
    }
    if( (max5 < fabs(q1_4_p0_4)) )
    {
        max5 = fabs(q1_4_p0_4);
    }
    if( (max5 < fabs(q1_5_p0_5)) )
    {
        max5 = fabs(q1_5_p0_5);
    }
    if( (max5 < fabs(q1_6_p0_6)) )
    {
        max5 = fabs(q1_6_p0_6);
    }
    if( (max5 < fabs(q1_7_p0_7)) )
    {
        max5 = fabs(q1_7_p0_7);
    }
    if( (max5 < fabs(q2_0_p0_0)) )
    {
        max5 = fabs(q2_0_p0_0);
    }
    if( (max5 < fabs(q2_1_p0_1)) )
    {
        max5 = fabs(q2_1_p0_1);
    }
    if( (max5 < fabs(q2_2_p0_2)) )
    {
        max5 = fabs(q2_2_p0_2);
    }
    if( (max5 < fabs(q2_3_p0_3)) )
    {
        max5 = fabs(q2_3_p0_3);
    }
    if( (max5 < fabs(q2_4_p0_4)) )
    {
        max5 = fabs(q2_4_p0_4);
    }
    if( (max5 < fabs(q2_5_p0_5)) )
    {
        max5 = fabs(q2_5_p0_5);
    }
    if( (max5 < fabs(q2_6_p0_6)) )
    {
        max5 = fabs(q2_6_p0_6);
    }
    if( (max5 < fabs(q2_7_p0_7)) )
    {
        max5 = fabs(q2_7_p0_7);
    }
    double max6 = fabs(q2_0_p0_0);
    if( (max6 < fabs(q2_1_p0_1)) )
    {
        max6 = fabs(q2_1_p0_1);
    }
    if( (max6 < fabs(q2_2_p0_2)) )
    {
        max6 = fabs(q2_2_p0_2);
    }
    if( (max6 < fabs(q2_3_p0_3)) )
    {
        max6 = fabs(q2_3_p0_3);
    }
    if( (max6 < fabs(q2_4_p0_4)) )
    {
        max6 = fabs(q2_4_p0_4);
    }
    if( (max6 < fabs(q2_5_p0_5)) )
    {
        max6 = fabs(q2_5_p0_5);
    }
    if( (max6 < fabs(q2_6_p0_6)) )
    {
        max6 = fabs(q2_6_p0_6);
    }
    if( (max6 < fabs(q2_7_p0_7)) )
    {
        max6 = fabs(q2_7_p0_7);
    }
    if( (max6 < fabs(q3_0_p0_0)) )
    {
        max6 = fabs(q3_0_p0_0);
    }
    if( (max6 < fabs(q3_1_p0_1)) )
    {
        max6 = fabs(q3_1_p0_1);
    }
    if( (max6 < fabs(q3_2_p0_2)) )
    {
        max6 = fabs(q3_2_p0_2);
    }
    if( (max6 < fabs(q3_3_p0_3)) )
    {
        max6 = fabs(q3_3_p0_3);
    }
    if( (max6 < fabs(q3_4_p0_4)) )
    {
        max6 = fabs(q3_4_p0_4);
    }
    if( (max6 < fabs(q3_5_p0_5)) )
    {
        max6 = fabs(q3_5_p0_5);
    }
    if( (max6 < fabs(q3_6_p0_6)) )
    {
        max6 = fabs(q3_6_p0_6);
    }
    if( (max6 < fabs(q3_7_p0_7)) )
    {
        max6 = fabs(q3_7_p0_7);
    }
    double lower_bound_1;
    double upper_bound_1;
    int Delta_sign;
    int int_tmp_result;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    }
    else
    {
        if( (max2 > upper_bound_1) )
        {
            upper_bound_1 = max2;
        }
    }
    if( (max3 < lower_bound_1) )
    {
        lower_bound_1 = max3;
    }
    else
    {
        if( (max3 > upper_bound_1) )
        {
            upper_bound_1 = max3;
        }
    }
    if( (max4 < lower_bound_1) )
    {
        lower_bound_1 = max4;
    }
    else
    {
        if( (max4 > upper_bound_1) )
        {
            upper_bound_1 = max4;
        }
    }
    if( (max5 < lower_bound_1) )
    {
        lower_bound_1 = max5;
    }
    else
    {
        if( (max5 > upper_bound_1) )
        {
            upper_bound_1 = max5;
        }
    }
    if( (max6 < lower_bound_1) )
    {
        lower_bound_1 = max6;
    }
    else
    {
        if( (max6 > upper_bound_1) )
        {
            upper_bound_1 = max6;
        }
    }
    if( (lower_bound_1 < 2.81560341188832891396e-50) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 4.83570327845851562508e+24) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (4.46596725665061863033e-11 * (((((max1 * max4) * max2) * max5) * max3) * max6));
        if( (Delta > eps) )
        {
            int_tmp_result = 1;
        }
        else
        {
            if( (Delta < -eps) )
            {
                int_tmp_result = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    Delta_sign = int_tmp_result;
    int int_tmp_result_k3Lzf6g;
    double max7 = max1;
    if( (max7 < max2) )
    {
        max7 = max2;
    }
    if( (max7 < max3) )
    {
        max7 = max3;
    }
    if( (max7 < max6) )
    {
        max7 = max6;
    }
    if( (max7 < fabs(p4_0_p0_0)) )
    {
        max7 = fabs(p4_0_p0_0);
    }
    if( (max7 < fabs(p4_1_p0_1)) )
    {
        max7 = fabs(p4_1_p0_1);
    }
    if( (max7 < fabs(p4_2_p0_2)) )
    {
        max7 = fabs(p4_2_p0_2);
    }
    if( (max7 < fabs(p4_3_p0_3)) )
    {
        max7 = fabs(p4_3_p0_3);
    }
    if( (max7 < fabs(p4_4_p0_4)) )
    {
        max7 = fabs(p4_4_p0_4);
    }
    if( (max7 < fabs(p4_5_p0_5)) )
    {
        max7 = fabs(p4_5_p0_5);
    }
    if( (max7 < fabs(p4_6_p0_6)) )
    {
        max7 = fabs(p4_6_p0_6);
    }
    if( (max7 < fabs(p4_7_p0_7)) )
    {
        max7 = fabs(p4_7_p0_7);
    }
    double max8 = max2;
    if( (max8 < max3) )
    {
        max8 = max3;
    }
    if( (max8 < max5) )
    {
        max8 = max5;
    }
    if( (max8 < max6) )
    {
        max8 = max6;
    }
    double max11 = max3;
    if( (max11 < max5) )
    {
        max11 = max5;
    }
    if( (max11 < max6) )
    {
        max11 = max6;
    }
    if( (max8 < max11) )
    {
        max8 = max11;
    }
    if( (max8 < fabs(p4_0_p0_0)) )
    {
        max8 = fabs(p4_0_p0_0);
    }
    if( (max8 < fabs(p4_1_p0_1)) )
    {
        max8 = fabs(p4_1_p0_1);
    }
    if( (max8 < fabs(p4_2_p0_2)) )
    {
        max8 = fabs(p4_2_p0_2);
    }
    if( (max8 < fabs(p4_3_p0_3)) )
    {
        max8 = fabs(p4_3_p0_3);
    }
    if( (max8 < fabs(p4_4_p0_4)) )
    {
        max8 = fabs(p4_4_p0_4);
    }
    if( (max8 < fabs(p4_5_p0_5)) )
    {
        max8 = fabs(p4_5_p0_5);
    }
    if( (max8 < fabs(p4_6_p0_6)) )
    {
        max8 = fabs(p4_6_p0_6);
    }
    if( (max8 < fabs(p4_7_p0_7)) )
    {
        max8 = fabs(p4_7_p0_7);
    }
    double max9 = max2;
    if( (max9 < max3) )
    {
        max9 = max3;
    }
    if( (max9 < max4) )
    {
        max9 = max4;
    }
    if( (max9 < max5) )
    {
        max9 = max5;
    }
    if( (max9 < max6) )
    {
        max9 = max6;
    }
    double max10 = max2;
    if( (max10 < max4) )
    {
        max10 = max4;
    }
    if( (max10 < max5) )
    {
        max10 = max5;
    }
    if( (max10 < max6) )
    {
        max10 = max6;
    }
    if( (max9 < max10) )
    {
        max9 = max10;
    }
    if( (max9 < max11) )
    {
        max9 = max11;
    }
    double max12 = max4;
    if( (max12 < fabs(p4_0_p0_0)) )
    {
        max12 = fabs(p4_0_p0_0);
    }
    if( (max12 < fabs(p4_1_p0_1)) )
    {
        max12 = fabs(p4_1_p0_1);
    }
    if( (max12 < fabs(p4_2_p0_2)) )
    {
        max12 = fabs(p4_2_p0_2);
    }
    if( (max12 < fabs(p4_3_p0_3)) )
    {
        max12 = fabs(p4_3_p0_3);
    }
    if( (max12 < fabs(p4_4_p0_4)) )
    {
        max12 = fabs(p4_4_p0_4);
    }
    if( (max12 < fabs(p4_5_p0_5)) )
    {
        max12 = fabs(p4_5_p0_5);
    }
    if( (max12 < fabs(p4_6_p0_6)) )
    {
        max12 = fabs(p4_6_p0_6);
    }
    if( (max12 < fabs(p4_7_p0_7)) )
    {
        max12 = fabs(p4_7_p0_7);
    }
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if( (max7 < lower_bound_1) )
    {
        lower_bound_1 = max7;
    }
    else
    {
        if( (max7 > upper_bound_1) )
        {
            upper_bound_1 = max7;
        }
    }
    if( (max8 < lower_bound_1) )
    {
        lower_bound_1 = max8;
    }
    else
    {
        if( (max8 > upper_bound_1) )
        {
            upper_bound_1 = max8;
        }
    }
    if( (max9 < lower_bound_1) )
    {
        lower_bound_1 = max9;
    }
    else
    {
        if( (max9 > upper_bound_1) )
        {
            upper_bound_1 = max9;
        }
    }
    if( (max10 < lower_bound_1) )
    {
        lower_bound_1 = max10;
    }
    if( (max11 < lower_bound_1) )
    {
        lower_bound_1 = max11;
    }
    if( (max12 < lower_bound_1) )
    {
        lower_bound_1 = max12;
    }
    else
    {
        if( (max12 > upper_bound_1) )
        {
            upper_bound_1 = max12;
        }
    }
    if( (lower_bound_1 < 4.16622355992800490501e-38) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 4.83570327845851562508e+24) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (2.45134177995386317815e-09 * (((((((max1 * max12) * max10) * max10) * max9) * max11) * max8) * max7));
        if( (r > eps) )
        {
            int_tmp_result_k3Lzf6g = 1;
        }
        else
        {
            if( (r < -eps) )
            {
                int_tmp_result_k3Lzf6g = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    return (Delta_sign * int_tmp_result_k3Lzf6g);
}

int avro_side5_8d_filter( const double* p0, const double* p1, const double* p2, const double* p3, const double* p4, const double* p5, const double* q0, const double* q1, const double* q2, const double* q3, const double* q4) {
    double p1_0_p0_0 = (p1[0] - p0[0]);
    double p1_1_p0_1 = (p1[1] - p0[1]);
    double p1_2_p0_2 = (p1[2] - p0[2]);
    double p1_3_p0_3 = (p1[3] - p0[3]);
    double p1_4_p0_4 = (p1[4] - p0[4]);
    double p1_5_p0_5 = (p1[5] - p0[5]);
    double p1_6_p0_6 = (p1[6] - p0[6]);
    double p1_7_p0_7 = (p1[7] - p0[7]);
    double l1;
    l1 = (1 * ((((((((p1_0_p0_0 * p1_0_p0_0) + (p1_1_p0_1 * p1_1_p0_1)) + (p1_2_p0_2 * p1_2_p0_2)) + (p1_3_p0_3 * p1_3_p0_3)) + (p1_4_p0_4 * p1_4_p0_4)) + (p1_5_p0_5 * p1_5_p0_5)) + (p1_6_p0_6 * p1_6_p0_6)) + (p1_7_p0_7 * p1_7_p0_7)));
    double p2_0_p0_0 = (p2[0] - p0[0]);
    double p2_1_p0_1 = (p2[1] - p0[1]);
    double p2_2_p0_2 = (p2[2] - p0[2]);
    double p2_3_p0_3 = (p2[3] - p0[3]);
    double p2_4_p0_4 = (p2[4] - p0[4]);
    double p2_5_p0_5 = (p2[5] - p0[5]);
    double p2_6_p0_6 = (p2[6] - p0[6]);
    double p2_7_p0_7 = (p2[7] - p0[7]);
    double l2;
    l2 = (1 * ((((((((p2_0_p0_0 * p2_0_p0_0) + (p2_1_p0_1 * p2_1_p0_1)) + (p2_2_p0_2 * p2_2_p0_2)) + (p2_3_p0_3 * p2_3_p0_3)) + (p2_4_p0_4 * p2_4_p0_4)) + (p2_5_p0_5 * p2_5_p0_5)) + (p2_6_p0_6 * p2_6_p0_6)) + (p2_7_p0_7 * p2_7_p0_7)));
    double p3_0_p0_0 = (p3[0] - p0[0]);
    double p3_1_p0_1 = (p3[1] - p0[1]);
    double p3_2_p0_2 = (p3[2] - p0[2]);
    double p3_3_p0_3 = (p3[3] - p0[3]);
    double p3_4_p0_4 = (p3[4] - p0[4]);
    double p3_5_p0_5 = (p3[5] - p0[5]);
    double p3_6_p0_6 = (p3[6] - p0[6]);
    double p3_7_p0_7 = (p3[7] - p0[7]);
    double l3;
    l3 = (1 * ((((((((p3_0_p0_0 * p3_0_p0_0) + (p3_1_p0_1 * p3_1_p0_1)) + (p3_2_p0_2 * p3_2_p0_2)) + (p3_3_p0_3 * p3_3_p0_3)) + (p3_4_p0_4 * p3_4_p0_4)) + (p3_5_p0_5 * p3_5_p0_5)) + (p3_6_p0_6 * p3_6_p0_6)) + (p3_7_p0_7 * p3_7_p0_7)));
    double p4_0_p0_0 = (p4[0] - p0[0]);
    double p4_1_p0_1 = (p4[1] - p0[1]);
    double p4_2_p0_2 = (p4[2] - p0[2]);
    double p4_3_p0_3 = (p4[3] - p0[3]);
    double p4_4_p0_4 = (p4[4] - p0[4]);
    double p4_5_p0_5 = (p4[5] - p0[5]);
    double p4_6_p0_6 = (p4[6] - p0[6]);
    double p4_7_p0_7 = (p4[7] - p0[7]);
    double l4;
    l4 = (1 * ((((((((p4_0_p0_0 * p4_0_p0_0) + (p4_1_p0_1 * p4_1_p0_1)) + (p4_2_p0_2 * p4_2_p0_2)) + (p4_3_p0_3 * p4_3_p0_3)) + (p4_4_p0_4 * p4_4_p0_4)) + (p4_5_p0_5 * p4_5_p0_5)) + (p4_6_p0_6 * p4_6_p0_6)) + (p4_7_p0_7 * p4_7_p0_7)));
    double p5_0_p0_0 = (p5[0] - p0[0]);
    double p5_1_p0_1 = (p5[1] - p0[1]);
    double p5_2_p0_2 = (p5[2] - p0[2]);
    double p5_3_p0_3 = (p5[3] - p0[3]);
    double p5_4_p0_4 = (p5[4] - p0[4]);
    double p5_5_p0_5 = (p5[5] - p0[5]);
    double p5_6_p0_6 = (p5[6] - p0[6]);
    double p5_7_p0_7 = (p5[7] - p0[7]);
    double l5;
    l5 = (1 * ((((((((p5_0_p0_0 * p5_0_p0_0) + (p5_1_p0_1 * p5_1_p0_1)) + (p5_2_p0_2 * p5_2_p0_2)) + (p5_3_p0_3 * p5_3_p0_3)) + (p5_4_p0_4 * p5_4_p0_4)) + (p5_5_p0_5 * p5_5_p0_5)) + (p5_6_p0_6 * p5_6_p0_6)) + (p5_7_p0_7 * p5_7_p0_7)));
    double q0_0_p0_0 = (q0[0] - p0[0]);
    double q0_1_p0_1 = (q0[1] - p0[1]);
    double q0_2_p0_2 = (q0[2] - p0[2]);
    double q0_3_p0_3 = (q0[3] - p0[3]);
    double q0_4_p0_4 = (q0[4] - p0[4]);
    double q0_5_p0_5 = (q0[5] - p0[5]);
    double q0_6_p0_6 = (q0[6] - p0[6]);
    double q0_7_p0_7 = (q0[7] - p0[7]);
    double a10;
    a10 = (2 * ((((((((p1_0_p0_0 * q0_0_p0_0) + (p1_1_p0_1 * q0_1_p0_1)) + (p1_2_p0_2 * q0_2_p0_2)) + (p1_3_p0_3 * q0_3_p0_3)) + (p1_4_p0_4 * q0_4_p0_4)) + (p1_5_p0_5 * q0_5_p0_5)) + (p1_6_p0_6 * q0_6_p0_6)) + (p1_7_p0_7 * q0_7_p0_7)));
    double q1_0_p0_0 = (q1[0] - p0[0]);
    double q1_1_p0_1 = (q1[1] - p0[1]);
    double q1_2_p0_2 = (q1[2] - p0[2]);
    double q1_3_p0_3 = (q1[3] - p0[3]);
    double q1_4_p0_4 = (q1[4] - p0[4]);
    double q1_5_p0_5 = (q1[5] - p0[5]);
    double q1_6_p0_6 = (q1[6] - p0[6]);
    double q1_7_p0_7 = (q1[7] - p0[7]);
    double a11;
    a11 = (2 * ((((((((p1_0_p0_0 * q1_0_p0_0) + (p1_1_p0_1 * q1_1_p0_1)) + (p1_2_p0_2 * q1_2_p0_2)) + (p1_3_p0_3 * q1_3_p0_3)) + (p1_4_p0_4 * q1_4_p0_4)) + (p1_5_p0_5 * q1_5_p0_5)) + (p1_6_p0_6 * q1_6_p0_6)) + (p1_7_p0_7 * q1_7_p0_7)));
    double q2_0_p0_0 = (q2[0] - p0[0]);
    double q2_1_p0_1 = (q2[1] - p0[1]);
    double q2_2_p0_2 = (q2[2] - p0[2]);
    double q2_3_p0_3 = (q2[3] - p0[3]);
    double q2_4_p0_4 = (q2[4] - p0[4]);
    double q2_5_p0_5 = (q2[5] - p0[5]);
    double q2_6_p0_6 = (q2[6] - p0[6]);
    double q2_7_p0_7 = (q2[7] - p0[7]);
    double a12;
    a12 = (2 * ((((((((p1_0_p0_0 * q2_0_p0_0) + (p1_1_p0_1 * q2_1_p0_1)) + (p1_2_p0_2 * q2_2_p0_2)) + (p1_3_p0_3 * q2_3_p0_3)) + (p1_4_p0_4 * q2_4_p0_4)) + (p1_5_p0_5 * q2_5_p0_5)) + (p1_6_p0_6 * q2_6_p0_6)) + (p1_7_p0_7 * q2_7_p0_7)));
    double q3_0_p0_0 = (q3[0] - p0[0]);
    double q3_1_p0_1 = (q3[1] - p0[1]);
    double q3_2_p0_2 = (q3[2] - p0[2]);
    double q3_3_p0_3 = (q3[3] - p0[3]);
    double q3_4_p0_4 = (q3[4] - p0[4]);
    double q3_5_p0_5 = (q3[5] - p0[5]);
    double q3_6_p0_6 = (q3[6] - p0[6]);
    double q3_7_p0_7 = (q3[7] - p0[7]);
    double a13;
    a13 = (2 * ((((((((p1_0_p0_0 * q3_0_p0_0) + (p1_1_p0_1 * q3_1_p0_1)) + (p1_2_p0_2 * q3_2_p0_2)) + (p1_3_p0_3 * q3_3_p0_3)) + (p1_4_p0_4 * q3_4_p0_4)) + (p1_5_p0_5 * q3_5_p0_5)) + (p1_6_p0_6 * q3_6_p0_6)) + (p1_7_p0_7 * q3_7_p0_7)));
    double q4_0_p0_0 = (q4[0] - p0[0]);
    double q4_1_p0_1 = (q4[1] - p0[1]);
    double q4_2_p0_2 = (q4[2] - p0[2]);
    double q4_3_p0_3 = (q4[3] - p0[3]);
    double q4_4_p0_4 = (q4[4] - p0[4]);
    double q4_5_p0_5 = (q4[5] - p0[5]);
    double q4_6_p0_6 = (q4[6] - p0[6]);
    double q4_7_p0_7 = (q4[7] - p0[7]);
    double a14;
    a14 = (2 * ((((((((p1_0_p0_0 * q4_0_p0_0) + (p1_1_p0_1 * q4_1_p0_1)) + (p1_2_p0_2 * q4_2_p0_2)) + (p1_3_p0_3 * q4_3_p0_3)) + (p1_4_p0_4 * q4_4_p0_4)) + (p1_5_p0_5 * q4_5_p0_5)) + (p1_6_p0_6 * q4_6_p0_6)) + (p1_7_p0_7 * q4_7_p0_7)));
    double a20;
    a20 = (2 * ((((((((p2_0_p0_0 * q0_0_p0_0) + (p2_1_p0_1 * q0_1_p0_1)) + (p2_2_p0_2 * q0_2_p0_2)) + (p2_3_p0_3 * q0_3_p0_3)) + (p2_4_p0_4 * q0_4_p0_4)) + (p2_5_p0_5 * q0_5_p0_5)) + (p2_6_p0_6 * q0_6_p0_6)) + (p2_7_p0_7 * q0_7_p0_7)));
    double a21;
    a21 = (2 * ((((((((p2_0_p0_0 * q1_0_p0_0) + (p2_1_p0_1 * q1_1_p0_1)) + (p2_2_p0_2 * q1_2_p0_2)) + (p2_3_p0_3 * q1_3_p0_3)) + (p2_4_p0_4 * q1_4_p0_4)) + (p2_5_p0_5 * q1_5_p0_5)) + (p2_6_p0_6 * q1_6_p0_6)) + (p2_7_p0_7 * q1_7_p0_7)));
    double a22;
    a22 = (2 * ((((((((p2_0_p0_0 * q2_0_p0_0) + (p2_1_p0_1 * q2_1_p0_1)) + (p2_2_p0_2 * q2_2_p0_2)) + (p2_3_p0_3 * q2_3_p0_3)) + (p2_4_p0_4 * q2_4_p0_4)) + (p2_5_p0_5 * q2_5_p0_5)) + (p2_6_p0_6 * q2_6_p0_6)) + (p2_7_p0_7 * q2_7_p0_7)));
    double a23;
    a23 = (2 * ((((((((p2_0_p0_0 * q3_0_p0_0) + (p2_1_p0_1 * q3_1_p0_1)) + (p2_2_p0_2 * q3_2_p0_2)) + (p2_3_p0_3 * q3_3_p0_3)) + (p2_4_p0_4 * q3_4_p0_4)) + (p2_5_p0_5 * q3_5_p0_5)) + (p2_6_p0_6 * q3_6_p0_6)) + (p2_7_p0_7 * q3_7_p0_7)));
    double a24;
    a24 = (2 * ((((((((p2_0_p0_0 * q4_0_p0_0) + (p2_1_p0_1 * q4_1_p0_1)) + (p2_2_p0_2 * q4_2_p0_2)) + (p2_3_p0_3 * q4_3_p0_3)) + (p2_4_p0_4 * q4_4_p0_4)) + (p2_5_p0_5 * q4_5_p0_5)) + (p2_6_p0_6 * q4_6_p0_6)) + (p2_7_p0_7 * q4_7_p0_7)));
    double a30;
    a30 = (2 * ((((((((p3_0_p0_0 * q0_0_p0_0) + (p3_1_p0_1 * q0_1_p0_1)) + (p3_2_p0_2 * q0_2_p0_2)) + (p3_3_p0_3 * q0_3_p0_3)) + (p3_4_p0_4 * q0_4_p0_4)) + (p3_5_p0_5 * q0_5_p0_5)) + (p3_6_p0_6 * q0_6_p0_6)) + (p3_7_p0_7 * q0_7_p0_7)));
    double a31;
    a31 = (2 * ((((((((p3_0_p0_0 * q1_0_p0_0) + (p3_1_p0_1 * q1_1_p0_1)) + (p3_2_p0_2 * q1_2_p0_2)) + (p3_3_p0_3 * q1_3_p0_3)) + (p3_4_p0_4 * q1_4_p0_4)) + (p3_5_p0_5 * q1_5_p0_5)) + (p3_6_p0_6 * q1_6_p0_6)) + (p3_7_p0_7 * q1_7_p0_7)));
    double a32;
    a32 = (2 * ((((((((p3_0_p0_0 * q2_0_p0_0) + (p3_1_p0_1 * q2_1_p0_1)) + (p3_2_p0_2 * q2_2_p0_2)) + (p3_3_p0_3 * q2_3_p0_3)) + (p3_4_p0_4 * q2_4_p0_4)) + (p3_5_p0_5 * q2_5_p0_5)) + (p3_6_p0_6 * q2_6_p0_6)) + (p3_7_p0_7 * q2_7_p0_7)));
    double a33;
    a33 = (2 * ((((((((p3_0_p0_0 * q3_0_p0_0) + (p3_1_p0_1 * q3_1_p0_1)) + (p3_2_p0_2 * q3_2_p0_2)) + (p3_3_p0_3 * q3_3_p0_3)) + (p3_4_p0_4 * q3_4_p0_4)) + (p3_5_p0_5 * q3_5_p0_5)) + (p3_6_p0_6 * q3_6_p0_6)) + (p3_7_p0_7 * q3_7_p0_7)));
    double a34;
    a34 = (2 * ((((((((p3_0_p0_0 * q4_0_p0_0) + (p3_1_p0_1 * q4_1_p0_1)) + (p3_2_p0_2 * q4_2_p0_2)) + (p3_3_p0_3 * q4_3_p0_3)) + (p3_4_p0_4 * q4_4_p0_4)) + (p3_5_p0_5 * q4_5_p0_5)) + (p3_6_p0_6 * q4_6_p0_6)) + (p3_7_p0_7 * q4_7_p0_7)));
    double a40;
    a40 = (2 * ((((((((p4_0_p0_0 * q0_0_p0_0) + (p4_1_p0_1 * q0_1_p0_1)) + (p4_2_p0_2 * q0_2_p0_2)) + (p4_3_p0_3 * q0_3_p0_3)) + (p4_4_p0_4 * q0_4_p0_4)) + (p4_5_p0_5 * q0_5_p0_5)) + (p4_6_p0_6 * q0_6_p0_6)) + (p4_7_p0_7 * q0_7_p0_7)));
    double a41;
    a41 = (2 * ((((((((p4_0_p0_0 * q1_0_p0_0) + (p4_1_p0_1 * q1_1_p0_1)) + (p4_2_p0_2 * q1_2_p0_2)) + (p4_3_p0_3 * q1_3_p0_3)) + (p4_4_p0_4 * q1_4_p0_4)) + (p4_5_p0_5 * q1_5_p0_5)) + (p4_6_p0_6 * q1_6_p0_6)) + (p4_7_p0_7 * q1_7_p0_7)));
    double a42;
    a42 = (2 * ((((((((p4_0_p0_0 * q2_0_p0_0) + (p4_1_p0_1 * q2_1_p0_1)) + (p4_2_p0_2 * q2_2_p0_2)) + (p4_3_p0_3 * q2_3_p0_3)) + (p4_4_p0_4 * q2_4_p0_4)) + (p4_5_p0_5 * q2_5_p0_5)) + (p4_6_p0_6 * q2_6_p0_6)) + (p4_7_p0_7 * q2_7_p0_7)));
    double a43;
    a43 = (2 * ((((((((p4_0_p0_0 * q3_0_p0_0) + (p4_1_p0_1 * q3_1_p0_1)) + (p4_2_p0_2 * q3_2_p0_2)) + (p4_3_p0_3 * q3_3_p0_3)) + (p4_4_p0_4 * q3_4_p0_4)) + (p4_5_p0_5 * q3_5_p0_5)) + (p4_6_p0_6 * q3_6_p0_6)) + (p4_7_p0_7 * q3_7_p0_7)));
    double a44;
    a44 = (2 * ((((((((p4_0_p0_0 * q4_0_p0_0) + (p4_1_p0_1 * q4_1_p0_1)) + (p4_2_p0_2 * q4_2_p0_2)) + (p4_3_p0_3 * q4_3_p0_3)) + (p4_4_p0_4 * q4_4_p0_4)) + (p4_5_p0_5 * q4_5_p0_5)) + (p4_6_p0_6 * q4_6_p0_6)) + (p4_7_p0_7 * q4_7_p0_7)));
    double a50;
    a50 = (2 * ((((((((p5_0_p0_0 * q0_0_p0_0) + (p5_1_p0_1 * q0_1_p0_1)) + (p5_2_p0_2 * q0_2_p0_2)) + (p5_3_p0_3 * q0_3_p0_3)) + (p5_4_p0_4 * q0_4_p0_4)) + (p5_5_p0_5 * q0_5_p0_5)) + (p5_6_p0_6 * q0_6_p0_6)) + (p5_7_p0_7 * q0_7_p0_7)));
    double a51;
    a51 = (2 * ((((((((p5_0_p0_0 * q1_0_p0_0) + (p5_1_p0_1 * q1_1_p0_1)) + (p5_2_p0_2 * q1_2_p0_2)) + (p5_3_p0_3 * q1_3_p0_3)) + (p5_4_p0_4 * q1_4_p0_4)) + (p5_5_p0_5 * q1_5_p0_5)) + (p5_6_p0_6 * q1_6_p0_6)) + (p5_7_p0_7 * q1_7_p0_7)));
    double a52;
    a52 = (2 * ((((((((p5_0_p0_0 * q2_0_p0_0) + (p5_1_p0_1 * q2_1_p0_1)) + (p5_2_p0_2 * q2_2_p0_2)) + (p5_3_p0_3 * q2_3_p0_3)) + (p5_4_p0_4 * q2_4_p0_4)) + (p5_5_p0_5 * q2_5_p0_5)) + (p5_6_p0_6 * q2_6_p0_6)) + (p5_7_p0_7 * q2_7_p0_7)));
    double a53;
    a53 = (2 * ((((((((p5_0_p0_0 * q3_0_p0_0) + (p5_1_p0_1 * q3_1_p0_1)) + (p5_2_p0_2 * q3_2_p0_2)) + (p5_3_p0_3 * q3_3_p0_3)) + (p5_4_p0_4 * q3_4_p0_4)) + (p5_5_p0_5 * q3_5_p0_5)) + (p5_6_p0_6 * q3_6_p0_6)) + (p5_7_p0_7 * q3_7_p0_7)));
    double a54;
    a54 = (2 * ((((((((p5_0_p0_0 * q4_0_p0_0) + (p5_1_p0_1 * q4_1_p0_1)) + (p5_2_p0_2 * q4_2_p0_2)) + (p5_3_p0_3 * q4_3_p0_3)) + (p5_4_p0_4 * q4_4_p0_4)) + (p5_5_p0_5 * q4_5_p0_5)) + (p5_6_p0_6 * q4_6_p0_6)) + (p5_7_p0_7 * q4_7_p0_7)));
    double b00;
    b00 = ((((((((((((((((((((((((((a11 * a22) * a33) * a44) - (((a11 * a22) * a34) * a43)) - (((a11 * a23) * a32) * a44)) + (((a11 * a23) * a34) * a42)) + (((a11 * a24) * a32) * a43)) - (((a11 * a24) * a33) * a42)) - (((a12 * a21) * a33) * a44)) + (((a12 * a21) * a34) * a43)) + (((a12 * a23) * a31) * a44)) - (((a12 * a23) * a34) * a41)) - (((a12 * a24) * a31) * a43)) + (((a12 * a24) * a33) * a41)) + (((a13 * a21) * a32) * a44)) - (((a13 * a21) * a34) * a42)) - (((a13 * a22) * a31) * a44)) + (((a13 * a22) * a34) * a41)) + (((a13 * a24) * a31) * a42)) - (((a13 * a24) * a32) * a41)) - (((a14 * a21) * a32) * a43)) + (((a14 * a21) * a33) * a42)) + (((a14 * a22) * a31) * a43)) - (((a14 * a22) * a33) * a41)) - (((a14 * a23) * a31) * a42)) + (((a14 * a23) * a32) * a41));
    double b01;
    b01 = -(((((((((((((((((((((((((a21 * a33) * a42) - ((a21 * a32) * a43)) + ((a22 * a31) * a43)) - ((a22 * a33) * a41)) - ((a23 * a31) * a42)) + ((a23 * a32) * a41)) + ((a21 * a32) * a44)) - ((a21 * a34) * a42)) - ((a22 * a31) * a44)) + ((a22 * a34) * a41)) + ((a24 * a31) * a42)) - ((a24 * a32) * a41)) - ((a21 * a33) * a44)) + ((a21 * a34) * a43)) + ((a23 * a31) * a44)) - ((a23 * a34) * a41)) - ((a24 * a31) * a43)) + ((a24 * a33) * a41)) + ((a22 * a33) * a44)) - ((a22 * a34) * a43)) - ((a23 * a32) * a44)) + ((a23 * a34) * a42)) + ((a24 * a32) * a43)) - ((a24 * a33) * a42));
    double b02;
    b02 = (((((((((((((((((((((((((a11 * a33) * a42) - ((a11 * a32) * a43)) + ((a12 * a31) * a43)) - ((a12 * a33) * a41)) - ((a13 * a31) * a42)) + ((a13 * a32) * a41)) + ((a11 * a32) * a44)) - ((a11 * a34) * a42)) - ((a12 * a31) * a44)) + ((a12 * a34) * a41)) + ((a14 * a31) * a42)) - ((a14 * a32) * a41)) - ((a11 * a33) * a44)) + ((a11 * a34) * a43)) + ((a13 * a31) * a44)) - ((a13 * a34) * a41)) - ((a14 * a31) * a43)) + ((a14 * a33) * a41)) + ((a12 * a33) * a44)) - ((a12 * a34) * a43)) - ((a13 * a32) * a44)) + ((a13 * a34) * a42)) + ((a14 * a32) * a43)) - ((a14 * a33) * a42));
    double b03;
    b03 = -(((((((((((((((((((((((((a11 * a23) * a42) - ((a11 * a22) * a43)) + ((a12 * a21) * a43)) - ((a12 * a23) * a41)) - ((a13 * a21) * a42)) + ((a13 * a22) * a41)) + ((a11 * a22) * a44)) - ((a11 * a24) * a42)) - ((a12 * a21) * a44)) + ((a12 * a24) * a41)) + ((a14 * a21) * a42)) - ((a14 * a22) * a41)) - ((a11 * a23) * a44)) + ((a11 * a24) * a43)) + ((a13 * a21) * a44)) - ((a13 * a24) * a41)) - ((a14 * a21) * a43)) + ((a14 * a23) * a41)) + ((a12 * a23) * a44)) - ((a12 * a24) * a43)) - ((a13 * a22) * a44)) + ((a13 * a24) * a42)) + ((a14 * a22) * a43)) - ((a14 * a23) * a42));
    double b04;
    b04 = (((((((((((((((((((((((((a11 * a23) * a32) - ((a11 * a22) * a33)) + ((a12 * a21) * a33)) - ((a12 * a23) * a31)) - ((a13 * a21) * a32)) + ((a13 * a22) * a31)) + ((a11 * a22) * a34)) - ((a11 * a24) * a32)) - ((a12 * a21) * a34)) + ((a12 * a24) * a31)) + ((a14 * a21) * a32)) - ((a14 * a22) * a31)) - ((a11 * a23) * a34)) + ((a11 * a24) * a33)) + ((a13 * a21) * a34)) - ((a13 * a24) * a31)) - ((a14 * a21) * a33)) + ((a14 * a23) * a31)) + ((a12 * a23) * a34)) - ((a12 * a24) * a33)) - ((a13 * a22) * a34)) + ((a13 * a24) * a32)) + ((a14 * a22) * a33)) - ((a14 * a23) * a32));
    double b10;
    b10 = -((((((((((((((((((((((((((a10 * a22) * a33) * a44) - (((a10 * a22) * a34) * a43)) - (((a10 * a23) * a32) * a44)) + (((a10 * a23) * a34) * a42)) + (((a10 * a24) * a32) * a43)) - (((a10 * a24) * a33) * a42)) - (((a12 * a20) * a33) * a44)) + (((a12 * a20) * a34) * a43)) + (((a12 * a23) * a30) * a44)) - (((a12 * a23) * a34) * a40)) - (((a12 * a24) * a30) * a43)) + (((a12 * a24) * a33) * a40)) + (((a13 * a20) * a32) * a44)) - (((a13 * a20) * a34) * a42)) - (((a13 * a22) * a30) * a44)) + (((a13 * a22) * a34) * a40)) + (((a13 * a24) * a30) * a42)) - (((a13 * a24) * a32) * a40)) - (((a14 * a20) * a32) * a43)) + (((a14 * a20) * a33) * a42)) + (((a14 * a22) * a30) * a43)) - (((a14 * a22) * a33) * a40)) - (((a14 * a23) * a30) * a42)) + (((a14 * a23) * a32) * a40));
    double b11;
    b11 = (((((((((((((((((((((((((a20 * a33) * a42) - ((a20 * a32) * a43)) + ((a22 * a30) * a43)) - ((a22 * a33) * a40)) - ((a23 * a30) * a42)) + ((a23 * a32) * a40)) + ((a20 * a32) * a44)) - ((a20 * a34) * a42)) - ((a22 * a30) * a44)) + ((a22 * a34) * a40)) + ((a24 * a30) * a42)) - ((a24 * a32) * a40)) - ((a20 * a33) * a44)) + ((a20 * a34) * a43)) + ((a23 * a30) * a44)) - ((a23 * a34) * a40)) - ((a24 * a30) * a43)) + ((a24 * a33) * a40)) + ((a22 * a33) * a44)) - ((a22 * a34) * a43)) - ((a23 * a32) * a44)) + ((a23 * a34) * a42)) + ((a24 * a32) * a43)) - ((a24 * a33) * a42));
    double b12;
    b12 = -(((((((((((((((((((((((((a10 * a33) * a42) - ((a10 * a32) * a43)) + ((a12 * a30) * a43)) - ((a12 * a33) * a40)) - ((a13 * a30) * a42)) + ((a13 * a32) * a40)) + ((a10 * a32) * a44)) - ((a10 * a34) * a42)) - ((a12 * a30) * a44)) + ((a12 * a34) * a40)) + ((a14 * a30) * a42)) - ((a14 * a32) * a40)) - ((a10 * a33) * a44)) + ((a10 * a34) * a43)) + ((a13 * a30) * a44)) - ((a13 * a34) * a40)) - ((a14 * a30) * a43)) + ((a14 * a33) * a40)) + ((a12 * a33) * a44)) - ((a12 * a34) * a43)) - ((a13 * a32) * a44)) + ((a13 * a34) * a42)) + ((a14 * a32) * a43)) - ((a14 * a33) * a42));
    double b13;
    b13 = (((((((((((((((((((((((((a10 * a23) * a42) - ((a10 * a22) * a43)) + ((a12 * a20) * a43)) - ((a12 * a23) * a40)) - ((a13 * a20) * a42)) + ((a13 * a22) * a40)) + ((a10 * a22) * a44)) - ((a10 * a24) * a42)) - ((a12 * a20) * a44)) + ((a12 * a24) * a40)) + ((a14 * a20) * a42)) - ((a14 * a22) * a40)) - ((a10 * a23) * a44)) + ((a10 * a24) * a43)) + ((a13 * a20) * a44)) - ((a13 * a24) * a40)) - ((a14 * a20) * a43)) + ((a14 * a23) * a40)) + ((a12 * a23) * a44)) - ((a12 * a24) * a43)) - ((a13 * a22) * a44)) + ((a13 * a24) * a42)) + ((a14 * a22) * a43)) - ((a14 * a23) * a42));
    double b14;
    b14 = -(((((((((((((((((((((((((a10 * a23) * a32) - ((a10 * a22) * a33)) + ((a12 * a20) * a33)) - ((a12 * a23) * a30)) - ((a13 * a20) * a32)) + ((a13 * a22) * a30)) + ((a10 * a22) * a34)) - ((a10 * a24) * a32)) - ((a12 * a20) * a34)) + ((a12 * a24) * a30)) + ((a14 * a20) * a32)) - ((a14 * a22) * a30)) - ((a10 * a23) * a34)) + ((a10 * a24) * a33)) + ((a13 * a20) * a34)) - ((a13 * a24) * a30)) - ((a14 * a20) * a33)) + ((a14 * a23) * a30)) + ((a12 * a23) * a34)) - ((a12 * a24) * a33)) - ((a13 * a22) * a34)) + ((a13 * a24) * a32)) + ((a14 * a22) * a33)) - ((a14 * a23) * a32));
    double b20;
    b20 = ((((((((((((((((((((((((((a10 * a21) * a33) * a44) - (((a10 * a21) * a34) * a43)) - (((a10 * a23) * a31) * a44)) + (((a10 * a23) * a34) * a41)) + (((a10 * a24) * a31) * a43)) - (((a10 * a24) * a33) * a41)) - (((a11 * a20) * a33) * a44)) + (((a11 * a20) * a34) * a43)) + (((a11 * a23) * a30) * a44)) - (((a11 * a23) * a34) * a40)) - (((a11 * a24) * a30) * a43)) + (((a11 * a24) * a33) * a40)) + (((a13 * a20) * a31) * a44)) - (((a13 * a20) * a34) * a41)) - (((a13 * a21) * a30) * a44)) + (((a13 * a21) * a34) * a40)) + (((a13 * a24) * a30) * a41)) - (((a13 * a24) * a31) * a40)) - (((a14 * a20) * a31) * a43)) + (((a14 * a20) * a33) * a41)) + (((a14 * a21) * a30) * a43)) - (((a14 * a21) * a33) * a40)) - (((a14 * a23) * a30) * a41)) + (((a14 * a23) * a31) * a40));
    double b21;
    b21 = -(((((((((((((((((((((((((a20 * a33) * a41) - ((a20 * a31) * a43)) + ((a21 * a30) * a43)) - ((a21 * a33) * a40)) - ((a23 * a30) * a41)) + ((a23 * a31) * a40)) + ((a20 * a31) * a44)) - ((a20 * a34) * a41)) - ((a21 * a30) * a44)) + ((a21 * a34) * a40)) + ((a24 * a30) * a41)) - ((a24 * a31) * a40)) - ((a20 * a33) * a44)) + ((a20 * a34) * a43)) + ((a23 * a30) * a44)) - ((a23 * a34) * a40)) - ((a24 * a30) * a43)) + ((a24 * a33) * a40)) + ((a21 * a33) * a44)) - ((a21 * a34) * a43)) - ((a23 * a31) * a44)) + ((a23 * a34) * a41)) + ((a24 * a31) * a43)) - ((a24 * a33) * a41));
    double b22;
    b22 = (((((((((((((((((((((((((a10 * a33) * a41) - ((a10 * a31) * a43)) + ((a11 * a30) * a43)) - ((a11 * a33) * a40)) - ((a13 * a30) * a41)) + ((a13 * a31) * a40)) + ((a10 * a31) * a44)) - ((a10 * a34) * a41)) - ((a11 * a30) * a44)) + ((a11 * a34) * a40)) + ((a14 * a30) * a41)) - ((a14 * a31) * a40)) - ((a10 * a33) * a44)) + ((a10 * a34) * a43)) + ((a13 * a30) * a44)) - ((a13 * a34) * a40)) - ((a14 * a30) * a43)) + ((a14 * a33) * a40)) + ((a11 * a33) * a44)) - ((a11 * a34) * a43)) - ((a13 * a31) * a44)) + ((a13 * a34) * a41)) + ((a14 * a31) * a43)) - ((a14 * a33) * a41));
    double b23;
    b23 = -(((((((((((((((((((((((((a10 * a23) * a41) - ((a10 * a21) * a43)) + ((a11 * a20) * a43)) - ((a11 * a23) * a40)) - ((a13 * a20) * a41)) + ((a13 * a21) * a40)) + ((a10 * a21) * a44)) - ((a10 * a24) * a41)) - ((a11 * a20) * a44)) + ((a11 * a24) * a40)) + ((a14 * a20) * a41)) - ((a14 * a21) * a40)) - ((a10 * a23) * a44)) + ((a10 * a24) * a43)) + ((a13 * a20) * a44)) - ((a13 * a24) * a40)) - ((a14 * a20) * a43)) + ((a14 * a23) * a40)) + ((a11 * a23) * a44)) - ((a11 * a24) * a43)) - ((a13 * a21) * a44)) + ((a13 * a24) * a41)) + ((a14 * a21) * a43)) - ((a14 * a23) * a41));
    double b24;
    b24 = (((((((((((((((((((((((((a10 * a23) * a31) - ((a10 * a21) * a33)) + ((a11 * a20) * a33)) - ((a11 * a23) * a30)) - ((a13 * a20) * a31)) + ((a13 * a21) * a30)) + ((a10 * a21) * a34)) - ((a10 * a24) * a31)) - ((a11 * a20) * a34)) + ((a11 * a24) * a30)) + ((a14 * a20) * a31)) - ((a14 * a21) * a30)) - ((a10 * a23) * a34)) + ((a10 * a24) * a33)) + ((a13 * a20) * a34)) - ((a13 * a24) * a30)) - ((a14 * a20) * a33)) + ((a14 * a23) * a30)) + ((a11 * a23) * a34)) - ((a11 * a24) * a33)) - ((a13 * a21) * a34)) + ((a13 * a24) * a31)) + ((a14 * a21) * a33)) - ((a14 * a23) * a31));
    double b30;
    b30 = -((((((((((((((((((((((((((a10 * a21) * a32) * a44) - (((a10 * a21) * a34) * a42)) - (((a10 * a22) * a31) * a44)) + (((a10 * a22) * a34) * a41)) + (((a10 * a24) * a31) * a42)) - (((a10 * a24) * a32) * a41)) - (((a11 * a20) * a32) * a44)) + (((a11 * a20) * a34) * a42)) + (((a11 * a22) * a30) * a44)) - (((a11 * a22) * a34) * a40)) - (((a11 * a24) * a30) * a42)) + (((a11 * a24) * a32) * a40)) + (((a12 * a20) * a31) * a44)) - (((a12 * a20) * a34) * a41)) - (((a12 * a21) * a30) * a44)) + (((a12 * a21) * a34) * a40)) + (((a12 * a24) * a30) * a41)) - (((a12 * a24) * a31) * a40)) - (((a14 * a20) * a31) * a42)) + (((a14 * a20) * a32) * a41)) + (((a14 * a21) * a30) * a42)) - (((a14 * a21) * a32) * a40)) - (((a14 * a22) * a30) * a41)) + (((a14 * a22) * a31) * a40));
    double b31;
    b31 = (((((((((((((((((((((((((a20 * a32) * a41) - ((a20 * a31) * a42)) + ((a21 * a30) * a42)) - ((a21 * a32) * a40)) - ((a22 * a30) * a41)) + ((a22 * a31) * a40)) + ((a20 * a31) * a44)) - ((a20 * a34) * a41)) - ((a21 * a30) * a44)) + ((a21 * a34) * a40)) + ((a24 * a30) * a41)) - ((a24 * a31) * a40)) - ((a20 * a32) * a44)) + ((a20 * a34) * a42)) + ((a22 * a30) * a44)) - ((a22 * a34) * a40)) - ((a24 * a30) * a42)) + ((a24 * a32) * a40)) + ((a21 * a32) * a44)) - ((a21 * a34) * a42)) - ((a22 * a31) * a44)) + ((a22 * a34) * a41)) + ((a24 * a31) * a42)) - ((a24 * a32) * a41));
    double b32;
    b32 = -(((((((((((((((((((((((((a10 * a32) * a41) - ((a10 * a31) * a42)) + ((a11 * a30) * a42)) - ((a11 * a32) * a40)) - ((a12 * a30) * a41)) + ((a12 * a31) * a40)) + ((a10 * a31) * a44)) - ((a10 * a34) * a41)) - ((a11 * a30) * a44)) + ((a11 * a34) * a40)) + ((a14 * a30) * a41)) - ((a14 * a31) * a40)) - ((a10 * a32) * a44)) + ((a10 * a34) * a42)) + ((a12 * a30) * a44)) - ((a12 * a34) * a40)) - ((a14 * a30) * a42)) + ((a14 * a32) * a40)) + ((a11 * a32) * a44)) - ((a11 * a34) * a42)) - ((a12 * a31) * a44)) + ((a12 * a34) * a41)) + ((a14 * a31) * a42)) - ((a14 * a32) * a41));
    double b33;
    b33 = (((((((((((((((((((((((((a10 * a22) * a41) - ((a10 * a21) * a42)) + ((a11 * a20) * a42)) - ((a11 * a22) * a40)) - ((a12 * a20) * a41)) + ((a12 * a21) * a40)) + ((a10 * a21) * a44)) - ((a10 * a24) * a41)) - ((a11 * a20) * a44)) + ((a11 * a24) * a40)) + ((a14 * a20) * a41)) - ((a14 * a21) * a40)) - ((a10 * a22) * a44)) + ((a10 * a24) * a42)) + ((a12 * a20) * a44)) - ((a12 * a24) * a40)) - ((a14 * a20) * a42)) + ((a14 * a22) * a40)) + ((a11 * a22) * a44)) - ((a11 * a24) * a42)) - ((a12 * a21) * a44)) + ((a12 * a24) * a41)) + ((a14 * a21) * a42)) - ((a14 * a22) * a41));
    double b34;
    b34 = -(((((((((((((((((((((((((a10 * a22) * a31) - ((a10 * a21) * a32)) + ((a11 * a20) * a32)) - ((a11 * a22) * a30)) - ((a12 * a20) * a31)) + ((a12 * a21) * a30)) + ((a10 * a21) * a34)) - ((a10 * a24) * a31)) - ((a11 * a20) * a34)) + ((a11 * a24) * a30)) + ((a14 * a20) * a31)) - ((a14 * a21) * a30)) - ((a10 * a22) * a34)) + ((a10 * a24) * a32)) + ((a12 * a20) * a34)) - ((a12 * a24) * a30)) - ((a14 * a20) * a32)) + ((a14 * a22) * a30)) + ((a11 * a22) * a34)) - ((a11 * a24) * a32)) - ((a12 * a21) * a34)) + ((a12 * a24) * a31)) + ((a14 * a21) * a32)) - ((a14 * a22) * a31));
    double b40;
    b40 = ((((((((((((((((((((((((((a10 * a21) * a32) * a43) - (((a10 * a21) * a33) * a42)) - (((a10 * a22) * a31) * a43)) + (((a10 * a22) * a33) * a41)) + (((a10 * a23) * a31) * a42)) - (((a10 * a23) * a32) * a41)) - (((a11 * a20) * a32) * a43)) + (((a11 * a20) * a33) * a42)) + (((a11 * a22) * a30) * a43)) - (((a11 * a22) * a33) * a40)) - (((a11 * a23) * a30) * a42)) + (((a11 * a23) * a32) * a40)) + (((a12 * a20) * a31) * a43)) - (((a12 * a20) * a33) * a41)) - (((a12 * a21) * a30) * a43)) + (((a12 * a21) * a33) * a40)) + (((a12 * a23) * a30) * a41)) - (((a12 * a23) * a31) * a40)) - (((a13 * a20) * a31) * a42)) + (((a13 * a20) * a32) * a41)) + (((a13 * a21) * a30) * a42)) - (((a13 * a21) * a32) * a40)) - (((a13 * a22) * a30) * a41)) + (((a13 * a22) * a31) * a40));
    double b41;
    b41 = -(((((((((((((((((((((((((a20 * a32) * a41) - ((a20 * a31) * a42)) + ((a21 * a30) * a42)) - ((a21 * a32) * a40)) - ((a22 * a30) * a41)) + ((a22 * a31) * a40)) + ((a20 * a31) * a43)) - ((a20 * a33) * a41)) - ((a21 * a30) * a43)) + ((a21 * a33) * a40)) + ((a23 * a30) * a41)) - ((a23 * a31) * a40)) - ((a20 * a32) * a43)) + ((a20 * a33) * a42)) + ((a22 * a30) * a43)) - ((a22 * a33) * a40)) - ((a23 * a30) * a42)) + ((a23 * a32) * a40)) + ((a21 * a32) * a43)) - ((a21 * a33) * a42)) - ((a22 * a31) * a43)) + ((a22 * a33) * a41)) + ((a23 * a31) * a42)) - ((a23 * a32) * a41));
    double b42;
    b42 = (((((((((((((((((((((((((a10 * a32) * a41) - ((a10 * a31) * a42)) + ((a11 * a30) * a42)) - ((a11 * a32) * a40)) - ((a12 * a30) * a41)) + ((a12 * a31) * a40)) + ((a10 * a31) * a43)) - ((a10 * a33) * a41)) - ((a11 * a30) * a43)) + ((a11 * a33) * a40)) + ((a13 * a30) * a41)) - ((a13 * a31) * a40)) - ((a10 * a32) * a43)) + ((a10 * a33) * a42)) + ((a12 * a30) * a43)) - ((a12 * a33) * a40)) - ((a13 * a30) * a42)) + ((a13 * a32) * a40)) + ((a11 * a32) * a43)) - ((a11 * a33) * a42)) - ((a12 * a31) * a43)) + ((a12 * a33) * a41)) + ((a13 * a31) * a42)) - ((a13 * a32) * a41));
    double b43;
    b43 = -(((((((((((((((((((((((((a10 * a22) * a41) - ((a10 * a21) * a42)) + ((a11 * a20) * a42)) - ((a11 * a22) * a40)) - ((a12 * a20) * a41)) + ((a12 * a21) * a40)) + ((a10 * a21) * a43)) - ((a10 * a23) * a41)) - ((a11 * a20) * a43)) + ((a11 * a23) * a40)) + ((a13 * a20) * a41)) - ((a13 * a21) * a40)) - ((a10 * a22) * a43)) + ((a10 * a23) * a42)) + ((a12 * a20) * a43)) - ((a12 * a23) * a40)) - ((a13 * a20) * a42)) + ((a13 * a22) * a40)) + ((a11 * a22) * a43)) - ((a11 * a23) * a42)) - ((a12 * a21) * a43)) + ((a12 * a23) * a41)) + ((a13 * a21) * a42)) - ((a13 * a22) * a41));
    double b44;
    b44 = (((((((((((((((((((((((((a10 * a22) * a31) - ((a10 * a21) * a32)) + ((a11 * a20) * a32)) - ((a11 * a22) * a30)) - ((a12 * a20) * a31)) + ((a12 * a21) * a30)) + ((a10 * a21) * a33)) - ((a10 * a23) * a31)) - ((a11 * a20) * a33)) + ((a11 * a23) * a30)) + ((a13 * a20) * a31)) - ((a13 * a21) * a30)) - ((a10 * a22) * a33)) + ((a10 * a23) * a32)) + ((a12 * a20) * a33)) - ((a12 * a23) * a30)) - ((a13 * a20) * a32)) + ((a13 * a22) * a30)) + ((a11 * a22) * a33)) - ((a11 * a23) * a32)) - ((a12 * a21) * a33)) + ((a12 * a23) * a31)) + ((a13 * a21) * a32)) - ((a13 * a22) * a31));
    double Delta;
    Delta = ((((b00 + b10) + b20) + b30) + b40);
    double DeltaLambda0;
    DeltaLambda0 = ((((b00 + (b01 * l1)) + (b02 * l2)) + (b03 * l3)) + (b04 * l4));
    double DeltaLambda1;
    DeltaLambda1 = ((((b10 + (b11 * l1)) + (b12 * l2)) + (b13 * l3)) + (b14 * l4));
    double DeltaLambda2;
    DeltaLambda2 = ((((b20 + (b21 * l1)) + (b22 * l2)) + (b23 * l3)) + (b24 * l4));
    double DeltaLambda3;
    DeltaLambda3 = ((((b30 + (b31 * l1)) + (b32 * l2)) + (b33 * l3)) + (b34 * l4));
    double DeltaLambda4;
    DeltaLambda4 = ((((b40 + (b41 * l1)) + (b42 * l2)) + (b43 * l3)) + (b44 * l4));
    double r;
    r = ((Delta * l5) - (((((a50 * DeltaLambda0) + (a51 * DeltaLambda1)) + (a52 * DeltaLambda2)) + (a53 * DeltaLambda3)) + (a54 * DeltaLambda4)));
    double eps;
    double max1 = fabs(p2_2_p0_2);
    if( (max1 < fabs(p2_4_p0_4)) )
    {
        max1 = fabs(p2_4_p0_4);
    }
    if( (max1 < fabs(p2_3_p0_3)) )
    {
        max1 = fabs(p2_3_p0_3);
    }
    if( (max1 < fabs(p2_1_p0_1)) )
    {
        max1 = fabs(p2_1_p0_1);
    }
    if( (max1 < fabs(p2_6_p0_6)) )
    {
        max1 = fabs(p2_6_p0_6);
    }
    if( (max1 < fabs(p2_7_p0_7)) )
    {
        max1 = fabs(p2_7_p0_7);
    }
    if( (max1 < fabs(p2_0_p0_0)) )
    {
        max1 = fabs(p2_0_p0_0);
    }
    if( (max1 < fabs(p2_5_p0_5)) )
    {
        max1 = fabs(p2_5_p0_5);
    }
    double max2 = fabs(p3_3_p0_3);
    if( (max2 < fabs(p3_1_p0_1)) )
    {
        max2 = fabs(p3_1_p0_1);
    }
    if( (max2 < fabs(p3_0_p0_0)) )
    {
        max2 = fabs(p3_0_p0_0);
    }
    if( (max2 < fabs(p3_2_p0_2)) )
    {
        max2 = fabs(p3_2_p0_2);
    }
    if( (max2 < fabs(p3_4_p0_4)) )
    {
        max2 = fabs(p3_4_p0_4);
    }
    if( (max2 < fabs(p3_5_p0_5)) )
    {
        max2 = fabs(p3_5_p0_5);
    }
    if( (max2 < fabs(p3_6_p0_6)) )
    {
        max2 = fabs(p3_6_p0_6);
    }
    if( (max2 < fabs(p3_7_p0_7)) )
    {
        max2 = fabs(p3_7_p0_7);
    }
    double max3 = fabs(p4_0_p0_0);
    if( (max3 < fabs(p4_7_p0_7)) )
    {
        max3 = fabs(p4_7_p0_7);
    }
    if( (max3 < fabs(p4_6_p0_6)) )
    {
        max3 = fabs(p4_6_p0_6);
    }
    if( (max3 < fabs(p4_1_p0_1)) )
    {
        max3 = fabs(p4_1_p0_1);
    }
    if( (max3 < fabs(p4_2_p0_2)) )
    {
        max3 = fabs(p4_2_p0_2);
    }
    if( (max3 < fabs(p4_3_p0_3)) )
    {
        max3 = fabs(p4_3_p0_3);
    }
    if( (max3 < fabs(p4_4_p0_4)) )
    {
        max3 = fabs(p4_4_p0_4);
    }
    if( (max3 < fabs(p4_5_p0_5)) )
    {
        max3 = fabs(p4_5_p0_5);
    }
    double max4 = fabs(p1_2_p0_2);
    if( (max4 < fabs(p1_0_p0_0)) )
    {
        max4 = fabs(p1_0_p0_0);
    }
    if( (max4 < fabs(p1_1_p0_1)) )
    {
        max4 = fabs(p1_1_p0_1);
    }
    if( (max4 < fabs(p1_3_p0_3)) )
    {
        max4 = fabs(p1_3_p0_3);
    }
    if( (max4 < fabs(p1_4_p0_4)) )
    {
        max4 = fabs(p1_4_p0_4);
    }
    if( (max4 < fabs(p1_5_p0_5)) )
    {
        max4 = fabs(p1_5_p0_5);
    }
    if( (max4 < fabs(p1_6_p0_6)) )
    {
        max4 = fabs(p1_6_p0_6);
    }
    if( (max4 < fabs(p1_7_p0_7)) )
    {
        max4 = fabs(p1_7_p0_7);
    }
    double max5 = fabs(q2_5_p0_5);
    if( (max5 < fabs(q2_6_p0_6)) )
    {
        max5 = fabs(q2_6_p0_6);
    }
    if( (max5 < fabs(q2_7_p0_7)) )
    {
        max5 = fabs(q2_7_p0_7);
    }
    if( (max5 < fabs(q3_0_p0_0)) )
    {
        max5 = fabs(q3_0_p0_0);
    }
    if( (max5 < fabs(q3_2_p0_2)) )
    {
        max5 = fabs(q3_2_p0_2);
    }
    if( (max5 < fabs(q3_3_p0_3)) )
    {
        max5 = fabs(q3_3_p0_3);
    }
    if( (max5 < fabs(q3_7_p0_7)) )
    {
        max5 = fabs(q3_7_p0_7);
    }
    if( (max5 < fabs(q3_4_p0_4)) )
    {
        max5 = fabs(q3_4_p0_4);
    }
    if( (max5 < fabs(q3_5_p0_5)) )
    {
        max5 = fabs(q3_5_p0_5);
    }
    if( (max5 < fabs(q2_0_p0_0)) )
    {
        max5 = fabs(q2_0_p0_0);
    }
    if( (max5 < fabs(q3_6_p0_6)) )
    {
        max5 = fabs(q3_6_p0_6);
    }
    if( (max5 < fabs(q2_1_p0_1)) )
    {
        max5 = fabs(q2_1_p0_1);
    }
    if( (max5 < fabs(q2_2_p0_2)) )
    {
        max5 = fabs(q2_2_p0_2);
    }
    if( (max5 < fabs(q2_3_p0_3)) )
    {
        max5 = fabs(q2_3_p0_3);
    }
    if( (max5 < fabs(q2_4_p0_4)) )
    {
        max5 = fabs(q2_4_p0_4);
    }
    if( (max5 < fabs(q3_1_p0_1)) )
    {
        max5 = fabs(q3_1_p0_1);
    }
    double max6 = fabs(q2_5_p0_5);
    if( (max6 < fabs(q2_6_p0_6)) )
    {
        max6 = fabs(q2_6_p0_6);
    }
    if( (max6 < fabs(q2_7_p0_7)) )
    {
        max6 = fabs(q2_7_p0_7);
    }
    if( (max6 < fabs(q1_4_p0_4)) )
    {
        max6 = fabs(q1_4_p0_4);
    }
    if( (max6 < fabs(q1_5_p0_5)) )
    {
        max6 = fabs(q1_5_p0_5);
    }
    if( (max6 < fabs(q1_6_p0_6)) )
    {
        max6 = fabs(q1_6_p0_6);
    }
    if( (max6 < fabs(q1_0_p0_0)) )
    {
        max6 = fabs(q1_0_p0_0);
    }
    if( (max6 < fabs(q1_7_p0_7)) )
    {
        max6 = fabs(q1_7_p0_7);
    }
    if( (max6 < fabs(q1_1_p0_1)) )
    {
        max6 = fabs(q1_1_p0_1);
    }
    if( (max6 < fabs(q2_0_p0_0)) )
    {
        max6 = fabs(q2_0_p0_0);
    }
    if( (max6 < fabs(q2_1_p0_1)) )
    {
        max6 = fabs(q2_1_p0_1);
    }
    if( (max6 < fabs(q1_2_p0_2)) )
    {
        max6 = fabs(q1_2_p0_2);
    }
    if( (max6 < fabs(q2_2_p0_2)) )
    {
        max6 = fabs(q2_2_p0_2);
    }
    if( (max6 < fabs(q2_3_p0_3)) )
    {
        max6 = fabs(q2_3_p0_3);
    }
    if( (max6 < fabs(q1_3_p0_3)) )
    {
        max6 = fabs(q1_3_p0_3);
    }
    if( (max6 < fabs(q2_4_p0_4)) )
    {
        max6 = fabs(q2_4_p0_4);
    }
    double max7 = fabs(q4_6_p0_6);
    if( (max7 < fabs(q3_0_p0_0)) )
    {
        max7 = fabs(q3_0_p0_0);
    }
    if( (max7 < fabs(q3_2_p0_2)) )
    {
        max7 = fabs(q3_2_p0_2);
    }
    if( (max7 < fabs(q3_3_p0_3)) )
    {
        max7 = fabs(q3_3_p0_3);
    }
    if( (max7 < fabs(q4_7_p0_7)) )
    {
        max7 = fabs(q4_7_p0_7);
    }
    if( (max7 < fabs(q3_7_p0_7)) )
    {
        max7 = fabs(q3_7_p0_7);
    }
    if( (max7 < fabs(q3_4_p0_4)) )
    {
        max7 = fabs(q3_4_p0_4);
    }
    if( (max7 < fabs(q4_0_p0_0)) )
    {
        max7 = fabs(q4_0_p0_0);
    }
    if( (max7 < fabs(q3_5_p0_5)) )
    {
        max7 = fabs(q3_5_p0_5);
    }
    if( (max7 < fabs(q3_6_p0_6)) )
    {
        max7 = fabs(q3_6_p0_6);
    }
    if( (max7 < fabs(q4_1_p0_1)) )
    {
        max7 = fabs(q4_1_p0_1);
    }
    if( (max7 < fabs(q4_2_p0_2)) )
    {
        max7 = fabs(q4_2_p0_2);
    }
    if( (max7 < fabs(q4_3_p0_3)) )
    {
        max7 = fabs(q4_3_p0_3);
    }
    if( (max7 < fabs(q4_4_p0_4)) )
    {
        max7 = fabs(q4_4_p0_4);
    }
    if( (max7 < fabs(q4_5_p0_5)) )
    {
        max7 = fabs(q4_5_p0_5);
    }
    if( (max7 < fabs(q3_1_p0_1)) )
    {
        max7 = fabs(q3_1_p0_1);
    }
    double max8 = fabs(q0_0_p0_0);
    if( (max8 < fabs(q0_1_p0_1)) )
    {
        max8 = fabs(q0_1_p0_1);
    }
    if( (max8 < fabs(q0_2_p0_2)) )
    {
        max8 = fabs(q0_2_p0_2);
    }
    if( (max8 < fabs(q0_3_p0_3)) )
    {
        max8 = fabs(q0_3_p0_3);
    }
    if( (max8 < fabs(q0_4_p0_4)) )
    {
        max8 = fabs(q0_4_p0_4);
    }
    if( (max8 < fabs(q1_4_p0_4)) )
    {
        max8 = fabs(q1_4_p0_4);
    }
    if( (max8 < fabs(q0_5_p0_5)) )
    {
        max8 = fabs(q0_5_p0_5);
    }
    if( (max8 < fabs(q1_5_p0_5)) )
    {
        max8 = fabs(q1_5_p0_5);
    }
    if( (max8 < fabs(q0_6_p0_6)) )
    {
        max8 = fabs(q0_6_p0_6);
    }
    if( (max8 < fabs(q0_7_p0_7)) )
    {
        max8 = fabs(q0_7_p0_7);
    }
    if( (max8 < fabs(q1_6_p0_6)) )
    {
        max8 = fabs(q1_6_p0_6);
    }
    if( (max8 < fabs(q1_0_p0_0)) )
    {
        max8 = fabs(q1_0_p0_0);
    }
    if( (max8 < fabs(q1_7_p0_7)) )
    {
        max8 = fabs(q1_7_p0_7);
    }
    if( (max8 < fabs(q1_1_p0_1)) )
    {
        max8 = fabs(q1_1_p0_1);
    }
    if( (max8 < fabs(q1_2_p0_2)) )
    {
        max8 = fabs(q1_2_p0_2);
    }
    if( (max8 < fabs(q1_3_p0_3)) )
    {
        max8 = fabs(q1_3_p0_3);
    }
    double lower_bound_1;
    double upper_bound_1;
    int Delta_sign;
    int int_tmp_result;
    lower_bound_1 = max3;
    upper_bound_1 = max3;
    if( (max5 < lower_bound_1) )
    {
        lower_bound_1 = max5;
    }
    else
    {
        if( (max5 > upper_bound_1) )
        {
            upper_bound_1 = max5;
        }
    }
    if( (max6 < lower_bound_1) )
    {
        lower_bound_1 = max6;
    }
    else
    {
        if( (max6 > upper_bound_1) )
        {
            upper_bound_1 = max6;
        }
    }
    if( (max7 < lower_bound_1) )
    {
        lower_bound_1 = max7;
    }
    else
    {
        if( (max7 > upper_bound_1) )
        {
            upper_bound_1 = max7;
        }
    }
    if( (max4 < lower_bound_1) )
    {
        lower_bound_1 = max4;
    }
    else
    {
        if( (max4 > upper_bound_1) )
        {
            upper_bound_1 = max4;
        }
    }
    if( (max8 < lower_bound_1) )
    {
        lower_bound_1 = max8;
    }
    else
    {
        if( (max8 > upper_bound_1) )
        {
            upper_bound_1 = max8;
        }
    }
    if( (max1 < lower_bound_1) )
    {
        lower_bound_1 = max1;
    }
    else
    {
        if( (max1 > upper_bound_1) )
        {
            upper_bound_1 = max1;
        }
    }
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    }
    else
    {
        if( (max2 > upper_bound_1) )
        {
            upper_bound_1 = max2;
        }
    }
    if( (lower_bound_1 < 4.09265118977974013133e-38) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 1.84467440737095475200e+19) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (2.82687480748000781186e-09 * (((((((max4 * max8) * max1) * max6) * max2) * max5) * max3) * max7));
        if( (Delta > eps) )
        {
            int_tmp_result = 1;
        }
        else
        {
            if( (Delta < -eps) )
            {
                int_tmp_result = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    Delta_sign = int_tmp_result;
    int int_tmp_result_Gx4H;
    double max9 = max3;
    if( (max9 < max7) )
    {
        max9 = max7;
    }
    if( (max9 < max4) )
    {
        max9 = max4;
    }
    if( (max9 < max1) )
    {
        max9 = max1;
    }
    if( (max9 < max2) )
    {
        max9 = max2;
    }
    if( (max9 < fabs(p5_1_p0_1)) )
    {
        max9 = fabs(p5_1_p0_1);
    }
    if( (max9 < fabs(p5_2_p0_2)) )
    {
        max9 = fabs(p5_2_p0_2);
    }
    if( (max9 < fabs(p5_3_p0_3)) )
    {
        max9 = fabs(p5_3_p0_3);
    }
    if( (max9 < fabs(p5_7_p0_7)) )
    {
        max9 = fabs(p5_7_p0_7);
    }
    if( (max9 < fabs(p5_4_p0_4)) )
    {
        max9 = fabs(p5_4_p0_4);
    }
    if( (max9 < fabs(p5_5_p0_5)) )
    {
        max9 = fabs(p5_5_p0_5);
    }
    if( (max9 < fabs(p5_6_p0_6)) )
    {
        max9 = fabs(p5_6_p0_6);
    }
    if( (max9 < fabs(p5_0_p0_0)) )
    {
        max9 = fabs(p5_0_p0_0);
    }
    double max10 = max3;
    if( (max10 < max5) )
    {
        max10 = max5;
    }
    if( (max10 < max6) )
    {
        max10 = max6;
    }
    if( (max10 < max7) )
    {
        max10 = max7;
    }
    double max13 = max3;
    if( (max13 < max5) )
    {
        max13 = max5;
    }
    if( (max13 < max6) )
    {
        max13 = max6;
    }
    if( (max13 < max7) )
    {
        max13 = max7;
    }
    double max14 = max5;
    if( (max14 < max6) )
    {
        max14 = max6;
    }
    if( (max14 < max7) )
    {
        max14 = max7;
    }
    if( (max14 < max2) )
    {
        max14 = max2;
    }
    if( (max13 < max14) )
    {
        max13 = max14;
    }
    double max15 = max3;
    if( (max15 < max5) )
    {
        max15 = max5;
    }
    if( (max15 < max6) )
    {
        max15 = max6;
    }
    if( (max15 < max7) )
    {
        max15 = max7;
    }
    if( (max13 < max15) )
    {
        max13 = max15;
    }
    if( (max13 < max2) )
    {
        max13 = max2;
    }
    if( (max10 < max13) )
    {
        max10 = max13;
    }
    if( (max10 < max14) )
    {
        max10 = max14;
    }
    if( (max10 < max15) )
    {
        max10 = max15;
    }
    if( (max10 < max1) )
    {
        max10 = max1;
    }
    if( (max10 < max2) )
    {
        max10 = max2;
    }
    if( (max10 < fabs(p5_1_p0_1)) )
    {
        max10 = fabs(p5_1_p0_1);
    }
    if( (max10 < fabs(p5_2_p0_2)) )
    {
        max10 = fabs(p5_2_p0_2);
    }
    if( (max10 < fabs(p5_3_p0_3)) )
    {
        max10 = fabs(p5_3_p0_3);
    }
    if( (max10 < fabs(p5_7_p0_7)) )
    {
        max10 = fabs(p5_7_p0_7);
    }
    if( (max10 < fabs(p5_4_p0_4)) )
    {
        max10 = fabs(p5_4_p0_4);
    }
    if( (max10 < fabs(p5_5_p0_5)) )
    {
        max10 = fabs(p5_5_p0_5);
    }
    if( (max10 < fabs(p5_6_p0_6)) )
    {
        max10 = fabs(p5_6_p0_6);
    }
    if( (max10 < fabs(p5_0_p0_0)) )
    {
        max10 = fabs(p5_0_p0_0);
    }
    double max11 = max5;
    if( (max11 < max6) )
    {
        max11 = max6;
    }
    if( (max11 < max7) )
    {
        max11 = max7;
    }
    if( (max11 < max8) )
    {
        max11 = max8;
    }
    double max12 = max5;
    if( (max12 < max6) )
    {
        max12 = max6;
    }
    if( (max12 < max7) )
    {
        max12 = max7;
    }
    if( (max12 < max8) )
    {
        max12 = max8;
    }
    if( (max12 < max1) )
    {
        max12 = max1;
    }
    if( (max11 < max12) )
    {
        max11 = max12;
    }
    if( (max11 < max14) )
    {
        max11 = max14;
    }
    if( (max11 < max1) )
    {
        max11 = max1;
    }
    if( (max11 < max2) )
    {
        max11 = max2;
    }
    double max16 = max8;
    if( (max16 < fabs(p5_1_p0_1)) )
    {
        max16 = fabs(p5_1_p0_1);
    }
    if( (max16 < fabs(p5_2_p0_2)) )
    {
        max16 = fabs(p5_2_p0_2);
    }
    if( (max16 < fabs(p5_3_p0_3)) )
    {
        max16 = fabs(p5_3_p0_3);
    }
    if( (max16 < fabs(p5_7_p0_7)) )
    {
        max16 = fabs(p5_7_p0_7);
    }
    if( (max16 < fabs(p5_4_p0_4)) )
    {
        max16 = fabs(p5_4_p0_4);
    }
    if( (max16 < fabs(p5_5_p0_5)) )
    {
        max16 = fabs(p5_5_p0_5);
    }
    if( (max16 < fabs(p5_6_p0_6)) )
    {
        max16 = fabs(p5_6_p0_6);
    }
    if( (max16 < fabs(p5_0_p0_0)) )
    {
        max16 = fabs(p5_0_p0_0);
    }
    lower_bound_1 = max4;
    upper_bound_1 = max4;
    if( (max9 < lower_bound_1) )
    {
        lower_bound_1 = max9;
    }
    else
    {
        if( (max9 > upper_bound_1) )
        {
            upper_bound_1 = max9;
        }
    }
    if( (max10 < lower_bound_1) )
    {
        lower_bound_1 = max10;
    }
    else
    {
        if( (max10 > upper_bound_1) )
        {
            upper_bound_1 = max10;
        }
    }
    if( (max11 < lower_bound_1) )
    {
        lower_bound_1 = max11;
    }
    else
    {
        if( (max11 > upper_bound_1) )
        {
            upper_bound_1 = max11;
        }
    }
    if( (max12 < lower_bound_1) )
    {
        lower_bound_1 = max12;
    }
    else
    {
        if( (max12 > upper_bound_1) )
        {
            upper_bound_1 = max12;
        }
    }
    if( (max13 < lower_bound_1) )
    {
        lower_bound_1 = max13;
    }
    if( (max14 < lower_bound_1) )
    {
        lower_bound_1 = max14;
    }
    if( (max15 < lower_bound_1) )
    {
        lower_bound_1 = max15;
    }
    if( (max16 < lower_bound_1) )
    {
        lower_bound_1 = max16;
    }
    else
    {
        if( (max16 > upper_bound_1) )
        {
            upper_bound_1 = max16;
        }
    }
    if( (lower_bound_1 < 8.16477434367976849772e-31) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 1.84467440737095475200e+19) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (1.69006173306177122172e-07 * (((((((((max4 * max16) * max12) * max12) * max11) * max14) * max13) * max15) * max10) * max9));
        if( (r > eps) )
        {
            int_tmp_result_Gx4H = 1;
        }
        else
        {
            if( (r < -eps) )
            {
                int_tmp_result_Gx4H = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    return (Delta_sign * int_tmp_result_Gx4H);
}

} // PCK

} // GEO
