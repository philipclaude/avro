#ifndef luma_NUMERICS_ORIENT4D_FILTER_H_
#define luma_NUMERICS_ORIENT4D_FILTER_H_

inline double
orient4dfast( const double *p0 , const double *p1 , const double* p2 , const double* p3 , const double* p4 )
{
  double a11;
  a11 = (p1[0] - p0[0]);
  double a12;
  a12 = (p1[1] - p0[1]);
  double a13;
  a13 = (p1[2] - p0[2]);
  double a14;
  a14 = (p1[3] - p0[3]);
  double a21;
  a21 = (p2[0] - p0[0]);
  double a22;
  a22 = (p2[1] - p0[1]);
  double a23;
  a23 = (p2[2] - p0[2]);
  double a24;
  a24 = (p2[3] - p0[3]);
  double a31;
  a31 = (p3[0] - p0[0]);
  double a32;
  a32 = (p3[1] - p0[1]);
  double a33;
  a33 = (p3[2] - p0[2]);
  double a34;
  a34 = (p3[3] - p0[3]);
  double a41;
  a41 = (p4[0] - p0[0]);
  double a42;
  a42 = (p4[1] - p0[1]);
  double a43;
  a43 = (p4[2] - p0[2]);
  double a44;
  a44 = (p4[3] - p0[3]);
  double Delta;
  Delta = ((((((((((((((((((((((((((a11 * a22) * a33) * a44) - (((a11 * a22) * a34) * a43)) - (((a11 * a23) * a32) * a44)) + (((a11 * a23) * a34) * a42)) + (((a11 * a24) * a32) * a43)) - (((a11 * a24) * a33) * a42)) - (((a12 * a21) * a33) * a44)) + (((a12 * a21) * a34) * a43)) + (((a12 * a23) * a31) * a44)) - (((a12 * a23) * a34) * a41)) - (((a12 * a24) * a31) * a43)) + (((a12 * a24) * a33) * a41)) + (((a13 * a21) * a32) * a44)) - (((a13 * a21) * a34) * a42)) - (((a13 * a22) * a31) * a44)) + (((a13 * a22) * a34) * a41)) + (((a13 * a24) * a31) * a42)) - (((a13 * a24) * a32) * a41)) - (((a14 * a21) * a32) * a43)) + (((a14 * a21) * a33) * a42)) + (((a14 * a22) * a31) * a43)) - (((a14 * a22) * a33) * a41)) - (((a14 * a23) * a31) * a42)) + (((a14 * a23) * a32) * a41));
  return Delta;
}

inline int
orient4d_filter( const double *p0 , const double *p1 , const double* p2 , const double *p3 , const double* p4 )
{
    double a11;
    a11 = (p1[0] - p0[0]);
    double a12;
    a12 = (p1[1] - p0[1]);
    double a13;
    a13 = (p1[2] - p0[2]);
    double a14;
    a14 = (p1[3] - p0[3]);
    double a21;
    a21 = (p2[0] - p0[0]);
    double a22;
    a22 = (p2[1] - p0[1]);
    double a23;
    a23 = (p2[2] - p0[2]);
    double a24;
    a24 = (p2[3] - p0[3]);
    double a31;
    a31 = (p3[0] - p0[0]);
    double a32;
    a32 = (p3[1] - p0[1]);
    double a33;
    a33 = (p3[2] - p0[2]);
    double a34;
    a34 = (p3[3] - p0[3]);
    double a41;
    a41 = (p4[0] - p0[0]);
    double a42;
    a42 = (p4[1] - p0[1]);
    double a43;
    a43 = (p4[2] - p0[2]);
    double a44;
    a44 = (p4[3] - p0[3]);
    double Delta;
    Delta = ((((((((((((((((((((((((((a11 * a22) * a33) * a44) - (((a11 * a22) * a34) * a43)) - (((a11 * a23) * a32) * a44)) + (((a11 * a23) * a34) * a42)) + (((a11 * a24) * a32) * a43)) - (((a11 * a24) * a33) * a42)) - (((a12 * a21) * a33) * a44)) + (((a12 * a21) * a34) * a43)) + (((a12 * a23) * a31) * a44)) - (((a12 * a23) * a34) * a41)) - (((a12 * a24) * a31) * a43)) + (((a12 * a24) * a33) * a41)) + (((a13 * a21) * a32) * a44)) - (((a13 * a21) * a34) * a42)) - (((a13 * a22) * a31) * a44)) + (((a13 * a22) * a34) * a41)) + (((a13 * a24) * a31) * a42)) - (((a13 * a24) * a32) * a41)) - (((a14 * a21) * a32) * a43)) + (((a14 * a21) * a33) * a42)) + (((a14 * a22) * a31) * a43)) - (((a14 * a22) * a33) * a41)) - (((a14 * a23) * a31) * a42)) + (((a14 * a23) * a32) * a41));
    int int_tmp_result;
    double eps;
    double max1 = fabs(a11);
    if( (max1 < fabs(a12)) )
    {
        max1 = fabs(a12);
    }
    if( (max1 < fabs(a13)) )
    {
        max1 = fabs(a13);
    }
    if( (max1 < fabs(a14)) )
    {
        max1 = fabs(a14);
    }
    double max2 = fabs(a21);
    if( (max2 < fabs(a22)) )
    {
        max2 = fabs(a22);
    }
    if( (max2 < fabs(a23)) )
    {
        max2 = fabs(a23);
    }
    if( (max2 < fabs(a24)) )
    {
        max2 = fabs(a24);
    }
    double max3 = fabs(a31);
    if( (max3 < fabs(a32)) )
    {
        max3 = fabs(a32);
    }
    if( (max3 < fabs(a33)) )
    {
        max3 = fabs(a33);
    }
    if( (max3 < fabs(a34)) )
    {
        max3 = fabs(a34);
    }
    double max4 = fabs(a41);
    if( (max4 < fabs(a42)) )
    {
        max4 = fabs(a42);
    }
    if( (max4 < fabs(a43)) )
    {
        max4 = fabs(a43);
    }
    if( (max4 < fabs(a44)) )
    {
        max4 = fabs(a44);
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
    if( (lower_bound_1 < 2.55895084730087935675e-74) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 7.23700557733225900010e+75) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (5.18915514634287573143e-14 * (((max1 * max2) * max3) * max4));
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
    return int_tmp_result;
}

#endif
