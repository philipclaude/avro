//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//

#include "numerics/predicates.h"

namespace GEO {

namespace PCK {

int avro_side1_3d_filter( const double* p0, const double* p1, const double* q0) {
    double p0_0_p1_0 = (p0[0] - p1[0]);
    double p0_1_p1_1 = (p0[1] - p1[1]);
    double p0_2_p1_2 = (p0[2] - p1[2]);
    double r;
    r = (1 * (((p0_0_p1_0 * p0_0_p1_0) + (p0_1_p1_1 * p0_1_p1_1)) + (p0_2_p1_2 * p0_2_p1_2)));
    double p1_0_p0_0 = (p1[0] - p0[0]);
    double q0_0_p0_0 = (q0[0] - p0[0]);
    double p1_1_p0_1 = (p1[1] - p0[1]);
    double q0_1_p0_1 = (q0[1] - p0[1]);
    double p1_2_p0_2 = (p1[2] - p0[2]);
    double q0_2_p0_2 = (q0[2] - p0[2]);
    r = (r - (2 * (((p1_0_p0_0 * q0_0_p0_0) + (p1_1_p0_1 * q0_1_p0_1)) + (p1_2_p0_2 * q0_2_p0_2))));
    int int_tmp_result;
    double eps;
    double max1 = fabs(p0_0_p1_0);
    if( (max1 < fabs(p0_1_p1_1)) )
    {
        max1 = fabs(p0_1_p1_1);
    }
    if( (max1 < fabs(p0_2_p1_2)) )
    {
        max1 = fabs(p0_2_p1_2);
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
    double max2 = fabs(p0_0_p1_0);
    if( (max2 < fabs(p0_1_p1_1)) )
    {
        max2 = fabs(p0_1_p1_1);
    }
    if( (max2 < fabs(p0_2_p1_2)) )
    {
        max2 = fabs(p0_2_p1_2);
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
    double lower_bound_1;
    double upper_bound_1;
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
    if( (lower_bound_1 < 2.23755023300058943229e-147) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 5.59936185544450928309e+101) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (4.44425370757048798480e-15 * (max1 * max2));
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


int avro_side2_3d_filter( const double* p0, const double* p1, const double* p2, const double* q0, const double* q1) {
    double p1_0_p0_0 = (p1[0] - p0[0]);
    double p1_1_p0_1 = (p1[1] - p0[1]);
    double p1_2_p0_2 = (p1[2] - p0[2]);
    double l1;
    l1 = (1 * (((p1_0_p0_0 * p1_0_p0_0) + (p1_1_p0_1 * p1_1_p0_1)) + (p1_2_p0_2 * p1_2_p0_2)));
    double p2_0_p0_0 = (p2[0] - p0[0]);
    double p2_1_p0_1 = (p2[1] - p0[1]);
    double p2_2_p0_2 = (p2[2] - p0[2]);
    double l2;
    l2 = (1 * (((p2_0_p0_0 * p2_0_p0_0) + (p2_1_p0_1 * p2_1_p0_1)) + (p2_2_p0_2 * p2_2_p0_2)));
    double q0_0_p0_0 = (q0[0] - p0[0]);
    double q0_1_p0_1 = (q0[1] - p0[1]);
    double q0_2_p0_2 = (q0[2] - p0[2]);
    double a10;
    a10 = (2 * (((p1_0_p0_0 * q0_0_p0_0) + (p1_1_p0_1 * q0_1_p0_1)) + (p1_2_p0_2 * q0_2_p0_2)));
    double q1_0_p0_0 = (q1[0] - p0[0]);
    double q1_1_p0_1 = (q1[1] - p0[1]);
    double q1_2_p0_2 = (q1[2] - p0[2]);
    double a11;
    a11 = (2 * (((p1_0_p0_0 * q1_0_p0_0) + (p1_1_p0_1 * q1_1_p0_1)) + (p1_2_p0_2 * q1_2_p0_2)));
    double a20;
    a20 = (2 * (((p2_0_p0_0 * q0_0_p0_0) + (p2_1_p0_1 * q0_1_p0_1)) + (p2_2_p0_2 * q0_2_p0_2)));
    double a21;
    a21 = (2 * (((p2_0_p0_0 * q1_0_p0_0) + (p2_1_p0_1 * q1_1_p0_1)) + (p2_2_p0_2 * q1_2_p0_2)));
    double Delta;
    Delta = (a11 - a10);
    double DeltaLambda0;
    DeltaLambda0 = (a11 - l1);
    double DeltaLambda1;
    DeltaLambda1 = (l1 - a10);
    double r;
    r = (((Delta * l2) - (a20 * DeltaLambda0)) - (a21 * DeltaLambda1));
    double eps;
    double max1 = fabs(q0_0_p0_0);
    if( (max1 < fabs(q0_1_p0_1)) )
    {
        max1 = fabs(q0_1_p0_1);
    }
    if( (max1 < fabs(q0_2_p0_2)) )
    {
        max1 = fabs(q0_2_p0_2);
    }
    if( (max1 < fabs(q1_0_p0_0)) )
    {
        max1 = fabs(q1_0_p0_0);
    }
    if( (max1 < fabs(q1_1_p0_1)) )
    {
        max1 = fabs(q1_1_p0_1);
    }
    if( (max1 < fabs(q1_2_p0_2)) )
    {
        max1 = fabs(q1_2_p0_2);
    }
    double max2 = fabs(p1_0_p0_0);
    if( (max2 < fabs(p1_1_p0_1)) )
    {
        max2 = fabs(p1_1_p0_1);
    }
    if( (max2 < fabs(p1_2_p0_2)) )
    {
        max2 = fabs(p1_2_p0_2);
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
    if( (lower_bound_1 < 2.23755023300058943229e-147) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 3.74144419156711063983e+50) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (4.44425370757048798480e-15 * (max2 * max1));
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
    double max4 = max1;
    if( (max4 < fabs(p2_0_p0_0)) )
    {
        max4 = fabs(p2_0_p0_0);
    }
    if( (max4 < fabs(p2_1_p0_1)) )
    {
        max4 = fabs(p2_1_p0_1);
    }
    if( (max4 < fabs(p2_2_p0_2)) )
    {
        max4 = fabs(p2_2_p0_2);
    }
    if( (max3 < max4) )
    {
        max3 = max4;
    }
    double r_sign;
    int int_tmp_result_FFWKCAA;
    lower_bound_1 = max3;
    upper_bound_1 = max3;
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    }
    if( (max4 < lower_bound_1) )
    {
        lower_bound_1 = max4;
    }
    if( (lower_bound_1 < 2.22985945097100191780e-74) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 3.74144419156711063983e+50) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (8.99983341597279045654e-14 * (((max2 * max4) * max4) * max3));
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
    lower_bound_1 = max3;
    upper_bound_1 = max3;
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    }
    if( (max4 < lower_bound_1) )
    {
        lower_bound_1 = max4;
    }
    if( (lower_bound_1 < 2.22985945097100191780e-74) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 3.74144419156711063983e+50) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (8.99983341597279045654e-14 * (((max2 * max4) * max4) * max3));
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
  //  return int_tmp_result_k60Ocge;
  r_sign = int_tmp_result_k60Ocge;
  return r_sign*Delta_sign;
}


int avro_side3_3d_filter( const double* p0, const double* p1, const double* p2, const double* p3, const double* q0, const double* q1, const double* q2) {
    double p1_0_p0_0 = (p1[0] - p0[0]);
    double p1_1_p0_1 = (p1[1] - p0[1]);
    double p1_2_p0_2 = (p1[2] - p0[2]);
    double l1;
    l1 = (1 * (((p1_0_p0_0 * p1_0_p0_0) + (p1_1_p0_1 * p1_1_p0_1)) + (p1_2_p0_2 * p1_2_p0_2)));
    double p2_0_p0_0 = (p2[0] - p0[0]);
    double p2_1_p0_1 = (p2[1] - p0[1]);
    double p2_2_p0_2 = (p2[2] - p0[2]);
    double l2;
    l2 = (1 * (((p2_0_p0_0 * p2_0_p0_0) + (p2_1_p0_1 * p2_1_p0_1)) + (p2_2_p0_2 * p2_2_p0_2)));
    double p3_0_p0_0 = (p3[0] - p0[0]);
    double p3_1_p0_1 = (p3[1] - p0[1]);
    double p3_2_p0_2 = (p3[2] - p0[2]);
    double l3;
    l3 = (1 * (((p3_0_p0_0 * p3_0_p0_0) + (p3_1_p0_1 * p3_1_p0_1)) + (p3_2_p0_2 * p3_2_p0_2)));
    double q0_0_p0_0 = (q0[0] - p0[0]);
    double q0_1_p0_1 = (q0[1] - p0[1]);
    double q0_2_p0_2 = (q0[2] - p0[2]);
    double a10;
    a10 = (2 * (((p1_0_p0_0 * q0_0_p0_0) + (p1_1_p0_1 * q0_1_p0_1)) + (p1_2_p0_2 * q0_2_p0_2)));
    double q1_0_p0_0 = (q1[0] - p0[0]);
    double q1_1_p0_1 = (q1[1] - p0[1]);
    double q1_2_p0_2 = (q1[2] - p0[2]);
    double a11;
    a11 = (2 * (((p1_0_p0_0 * q1_0_p0_0) + (p1_1_p0_1 * q1_1_p0_1)) + (p1_2_p0_2 * q1_2_p0_2)));
    double q2_0_p0_0 = (q2[0] - p0[0]);
    double q2_1_p0_1 = (q2[1] - p0[1]);
    double q2_2_p0_2 = (q2[2] - p0[2]);
    double a12;
    a12 = (2 * (((p1_0_p0_0 * q2_0_p0_0) + (p1_1_p0_1 * q2_1_p0_1)) + (p1_2_p0_2 * q2_2_p0_2)));
    double a20;
    a20 = (2 * (((p2_0_p0_0 * q0_0_p0_0) + (p2_1_p0_1 * q0_1_p0_1)) + (p2_2_p0_2 * q0_2_p0_2)));
    double a21;
    a21 = (2 * (((p2_0_p0_0 * q1_0_p0_0) + (p2_1_p0_1 * q1_1_p0_1)) + (p2_2_p0_2 * q1_2_p0_2)));
    double a22;
    a22 = (2 * (((p2_0_p0_0 * q2_0_p0_0) + (p2_1_p0_1 * q2_1_p0_1)) + (p2_2_p0_2 * q2_2_p0_2)));
    double a30;
    a30 = (2 * (((p3_0_p0_0 * q0_0_p0_0) + (p3_1_p0_1 * q0_1_p0_1)) + (p3_2_p0_2 * q0_2_p0_2)));
    double a31;
    a31 = (2 * (((p3_0_p0_0 * q1_0_p0_0) + (p3_1_p0_1 * q1_1_p0_1)) + (p3_2_p0_2 * q1_2_p0_2)));
    double a32;
    a32 = (2 * (((p3_0_p0_0 * q2_0_p0_0) + (p3_1_p0_1 * q2_1_p0_1)) + (p3_2_p0_2 * q2_2_p0_2)));
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
    double max1 = fabs(q0_0_p0_0);
    if( (max1 < fabs(q0_1_p0_1)) )
    {
        max1 = fabs(q0_1_p0_1);
    }
    if( (max1 < fabs(q0_2_p0_2)) )
    {
        max1 = fabs(q0_2_p0_2);
    }
    if( (max1 < fabs(q1_0_p0_0)) )
    {
        max1 = fabs(q1_0_p0_0);
    }
    if( (max1 < fabs(q1_1_p0_1)) )
    {
        max1 = fabs(q1_1_p0_1);
    }
    if( (max1 < fabs(q1_2_p0_2)) )
    {
        max1 = fabs(q1_2_p0_2);
    }
    double max2 = fabs(p1_2_p0_2);
    if( (max2 < fabs(p1_0_p0_0)) )
    {
        max2 = fabs(p1_0_p0_0);
    }
    if( (max2 < fabs(p1_1_p0_1)) )
    {
        max2 = fabs(p1_1_p0_1);
    }
    double max3 = fabs(p2_0_p0_0);
    if( (max3 < fabs(p2_1_p0_1)) )
    {
        max3 = fabs(p2_1_p0_1);
    }
    if( (max3 < fabs(p2_2_p0_2)) )
    {
        max3 = fabs(p2_2_p0_2);
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
    double lower_bound_1;
    double upper_bound_1;
    int Delta_sign;
    int int_tmp_result;
    lower_bound_1 = max2;
    upper_bound_1 = max2;
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
    if( (lower_bound_1 < 2.22985945097100191780e-74) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 2.59614842926741294957e+33) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (8.99983341597279045654e-14 * (((max2 * max1) * max3) * max4));
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
    double max5 = max3;
    if( (max5 < max1) )
    {
        max5 = max1;
    }
    if( (max5 < max4) )
    {
        max5 = max4;
    }
    double max6 = max1;
    if( (max6 < fabs(p3_1_p0_1)) )
    {
        max6 = fabs(p3_1_p0_1);
    }
    if( (max6 < fabs(p3_0_p0_0)) )
    {
        max6 = fabs(p3_0_p0_0);
    }
    if( (max6 < fabs(p3_2_p0_2)) )
    {
        max6 = fabs(p3_2_p0_2);
    }
    if( (max5 < max6) )
    {
        max5 = max6;
    }
    double max8 = max3;
    if( (max8 < max1) )
    {
        max8 = max1;
    }
    if( (max8 < max4) )
    {
        max8 = max4;
    }
    if( (max5 < max8) )
    {
        max5 = max8;
    }
    double max7 = max2;
    if( (max7 < max3) )
    {
        max7 = max3;
    }
    if( (max7 < max4) )
    {
        max7 = max4;
    }
    if( (max7 < fabs(p3_1_p0_1)) )
    {
        max7 = fabs(p3_1_p0_1);
    }
    if( (max7 < fabs(p3_0_p0_0)) )
    {
        max7 = fabs(p3_0_p0_0);
    }
    if( (max7 < fabs(p3_2_p0_2)) )
    {
        max7 = fabs(p3_2_p0_2);
    }
    lower_bound_1 = max2;
    upper_bound_1 = max2;
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
    if( (max8 < lower_bound_1) )
    {
        lower_bound_1 = max8;
    }
    if( (lower_bound_1 < 4.84416636653081796592e-50) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 2.59614842926741294957e+33) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (1.72198804259438718181e-12 * (((((max2 * max6) * max8) * max8) * max5) * max7));
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


int avro_side4_3d_filter( const double* p0, const double* p1, const double* p2, const double* p3, const double* p4, const double* q0, const double* q1, const double* q2, const double* q3) {
    double p1_0_p0_0 = (p1[0] - p0[0]);
    double p1_1_p0_1 = (p1[1] - p0[1]);
    double p1_2_p0_2 = (p1[2] - p0[2]);
    double l1;
    l1 = (1 * (((p1_0_p0_0 * p1_0_p0_0) + (p1_1_p0_1 * p1_1_p0_1)) + (p1_2_p0_2 * p1_2_p0_2)));
    double p2_0_p0_0 = (p2[0] - p0[0]);
    double p2_1_p0_1 = (p2[1] - p0[1]);
    double p2_2_p0_2 = (p2[2] - p0[2]);
    double l2;
    l2 = (1 * (((p2_0_p0_0 * p2_0_p0_0) + (p2_1_p0_1 * p2_1_p0_1)) + (p2_2_p0_2 * p2_2_p0_2)));
    double p3_0_p0_0 = (p3[0] - p0[0]);
    double p3_1_p0_1 = (p3[1] - p0[1]);
    double p3_2_p0_2 = (p3[2] - p0[2]);
    double l3;
    l3 = (1 * (((p3_0_p0_0 * p3_0_p0_0) + (p3_1_p0_1 * p3_1_p0_1)) + (p3_2_p0_2 * p3_2_p0_2)));
    double p4_0_p0_0 = (p4[0] - p0[0]);
    double p4_1_p0_1 = (p4[1] - p0[1]);
    double p4_2_p0_2 = (p4[2] - p0[2]);
    double l4;
    l4 = (1 * (((p4_0_p0_0 * p4_0_p0_0) + (p4_1_p0_1 * p4_1_p0_1)) + (p4_2_p0_2 * p4_2_p0_2)));
    double q0_0_p0_0 = (q0[0] - p0[0]);
    double q0_1_p0_1 = (q0[1] - p0[1]);
    double q0_2_p0_2 = (q0[2] - p0[2]);
    double a10;
    a10 = (2 * (((p1_0_p0_0 * q0_0_p0_0) + (p1_1_p0_1 * q0_1_p0_1)) + (p1_2_p0_2 * q0_2_p0_2)));
    double q1_0_p0_0 = (q1[0] - p0[0]);
    double q1_1_p0_1 = (q1[1] - p0[1]);
    double q1_2_p0_2 = (q1[2] - p0[2]);
    double a11;
    a11 = (2 * (((p1_0_p0_0 * q1_0_p0_0) + (p1_1_p0_1 * q1_1_p0_1)) + (p1_2_p0_2 * q1_2_p0_2)));
    double q2_0_p0_0 = (q2[0] - p0[0]);
    double q2_1_p0_1 = (q2[1] - p0[1]);
    double q2_2_p0_2 = (q2[2] - p0[2]);
    double a12;
    a12 = (2 * (((p1_0_p0_0 * q2_0_p0_0) + (p1_1_p0_1 * q2_1_p0_1)) + (p1_2_p0_2 * q2_2_p0_2)));
    double q3_0_p0_0 = (q3[0] - p0[0]);
    double q3_1_p0_1 = (q3[1] - p0[1]);
    double q3_2_p0_2 = (q3[2] - p0[2]);
    double a13;
    a13 = (2 * (((p1_0_p0_0 * q3_0_p0_0) + (p1_1_p0_1 * q3_1_p0_1)) + (p1_2_p0_2 * q3_2_p0_2)));
    double a20;
    a20 = (2 * (((p2_0_p0_0 * q0_0_p0_0) + (p2_1_p0_1 * q0_1_p0_1)) + (p2_2_p0_2 * q0_2_p0_2)));
    double a21;
    a21 = (2 * (((p2_0_p0_0 * q1_0_p0_0) + (p2_1_p0_1 * q1_1_p0_1)) + (p2_2_p0_2 * q1_2_p0_2)));
    double a22;
    a22 = (2 * (((p2_0_p0_0 * q2_0_p0_0) + (p2_1_p0_1 * q2_1_p0_1)) + (p2_2_p0_2 * q2_2_p0_2)));
    double a23;
    a23 = (2 * (((p2_0_p0_0 * q3_0_p0_0) + (p2_1_p0_1 * q3_1_p0_1)) + (p2_2_p0_2 * q3_2_p0_2)));
    double a30;
    a30 = (2 * (((p3_0_p0_0 * q0_0_p0_0) + (p3_1_p0_1 * q0_1_p0_1)) + (p3_2_p0_2 * q0_2_p0_2)));
    double a31;
    a31 = (2 * (((p3_0_p0_0 * q1_0_p0_0) + (p3_1_p0_1 * q1_1_p0_1)) + (p3_2_p0_2 * q1_2_p0_2)));
    double a32;
    a32 = (2 * (((p3_0_p0_0 * q2_0_p0_0) + (p3_1_p0_1 * q2_1_p0_1)) + (p3_2_p0_2 * q2_2_p0_2)));
    double a33;
    a33 = (2 * (((p3_0_p0_0 * q3_0_p0_0) + (p3_1_p0_1 * q3_1_p0_1)) + (p3_2_p0_2 * q3_2_p0_2)));
    double a40;
    a40 = (2 * (((p4_0_p0_0 * q0_0_p0_0) + (p4_1_p0_1 * q0_1_p0_1)) + (p4_2_p0_2 * q0_2_p0_2)));
    double a41;
    a41 = (2 * (((p4_0_p0_0 * q1_0_p0_0) + (p4_1_p0_1 * q1_1_p0_1)) + (p4_2_p0_2 * q1_2_p0_2)));
    double a42;
    a42 = (2 * (((p4_0_p0_0 * q2_0_p0_0) + (p4_1_p0_1 * q2_1_p0_1)) + (p4_2_p0_2 * q2_2_p0_2)));
    double a43;
    a43 = (2 * (((p4_0_p0_0 * q3_0_p0_0) + (p4_1_p0_1 * q3_1_p0_1)) + (p4_2_p0_2 * q3_2_p0_2)));
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
    double max1 = fabs(p2_0_p0_0);
    if( (max1 < fabs(p2_2_p0_2)) )
    {
        max1 = fabs(p2_2_p0_2);
    }
    if( (max1 < fabs(p2_1_p0_1)) )
    {
        max1 = fabs(p2_1_p0_1);
    }
    double max2 = fabs(p1_1_p0_1);
    if( (max2 < fabs(p1_0_p0_0)) )
    {
        max2 = fabs(p1_0_p0_0);
    }
    if( (max2 < fabs(p1_2_p0_2)) )
    {
        max2 = fabs(p1_2_p0_2);
    }
    double max3 = fabs(p3_0_p0_0);
    if( (max3 < fabs(p3_1_p0_1)) )
    {
        max3 = fabs(p3_1_p0_1);
    }
    if( (max3 < fabs(p3_2_p0_2)) )
    {
        max3 = fabs(p3_2_p0_2);
    }
    double max4 = fabs(q1_1_p0_1);
    if( (max4 < fabs(q0_0_p0_0)) )
    {
        max4 = fabs(q0_0_p0_0);
    }
    if( (max4 < fabs(q0_1_p0_1)) )
    {
        max4 = fabs(q0_1_p0_1);
    }
    if( (max4 < fabs(q0_2_p0_2)) )
    {
        max4 = fabs(q0_2_p0_2);
    }
    if( (max4 < fabs(q1_0_p0_0)) )
    {
        max4 = fabs(q1_0_p0_0);
    }
    if( (max4 < fabs(q1_2_p0_2)) )
    {
        max4 = fabs(q1_2_p0_2);
    }
    double max5 = fabs(q1_1_p0_1);
    if( (max5 < fabs(q1_0_p0_0)) )
    {
        max5 = fabs(q1_0_p0_0);
    }
    if( (max5 < fabs(q1_2_p0_2)) )
    {
        max5 = fabs(q1_2_p0_2);
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
    double max6 = fabs(q2_0_p0_0);
    if( (max6 < fabs(q2_1_p0_1)) )
    {
        max6 = fabs(q2_1_p0_1);
    }
    if( (max6 < fabs(q2_2_p0_2)) )
    {
        max6 = fabs(q2_2_p0_2);
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
    double lower_bound_1;
    double upper_bound_1;
    int Delta_sign;
    int int_tmp_result;
    lower_bound_1 = max3;
    upper_bound_1 = max3;
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
    if( (lower_bound_1 < 4.82201624962410495789e-50) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 4.83570327845851562508e+24) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (1.76999652622699071193e-12 * (((((max2 * max4) * max1) * max5) * max3) * max6));
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
    double max7 = max3;
    if( (max7 < max1) )
    {
        max7 = max1;
    }
    if( (max7 < max2) )
    {
        max7 = max2;
    }
    if( (max7 < max6) )
    {
        max7 = max6;
    }
    if( (max7 < fabs(p4_2_p0_2)) )
    {
        max7 = fabs(p4_2_p0_2);
    }
    if( (max7 < fabs(p4_0_p0_0)) )
    {
        max7 = fabs(p4_0_p0_0);
    }
    if( (max7 < fabs(p4_1_p0_1)) )
    {
        max7 = fabs(p4_1_p0_1);
    }
    double max8 = max3;
    if( (max8 < max5) )
    {
        max8 = max5;
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
    if( (max8 < max1) )
    {
        max8 = max1;
    }
    if( (max8 < max6) )
    {
        max8 = max6;
    }
    if( (max8 < fabs(p4_2_p0_2)) )
    {
        max8 = fabs(p4_2_p0_2);
    }
    if( (max8 < fabs(p4_0_p0_0)) )
    {
        max8 = fabs(p4_0_p0_0);
    }
    if( (max8 < fabs(p4_1_p0_1)) )
    {
        max8 = fabs(p4_1_p0_1);
    }
    double max9 = max3;
    double max10 = max4;
    if( (max10 < max5) )
    {
        max10 = max5;
    }
    if( (max10 < max1) )
    {
        max10 = max1;
    }
    if( (max10 < max6) )
    {
        max10 = max6;
    }
    if( (max9 < max10) )
    {
        max9 = max10;
    }
    if( (max9 < max4) )
    {
        max9 = max4;
    }
    if( (max9 < max5) )
    {
        max9 = max5;
    }
    if( (max9 < max11) )
    {
        max9 = max11;
    }
    if( (max9 < max1) )
    {
        max9 = max1;
    }
    if( (max9 < max6) )
    {
        max9 = max6;
    }
    double max12 = max4;
    if( (max12 < fabs(p4_2_p0_2)) )
    {
        max12 = fabs(p4_2_p0_2);
    }
    if( (max12 < fabs(p4_0_p0_0)) )
    {
        max12 = fabs(p4_0_p0_0);
    }
    if( (max12 < fabs(p4_1_p0_1)) )
    {
        max12 = fabs(p4_1_p0_1);
    }
    lower_bound_1 = max10;
    upper_bound_1 = max10;
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
    if( (max11 < lower_bound_1) )
    {
        lower_bound_1 = max11;
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
    if( (lower_bound_1 < 7.04247297406605667281e-38) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 4.83570327845851562508e+24) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (3.67741733685231585719e-11 * (((((((max2 * max12) * max10) * max10) * max9) * max11) * max8) * max7));
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

// side5 in 3d doesn't make sense since simplex dimension > ambient dimension

} // PCK

} // GEO
