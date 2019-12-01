// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXD_DETERMINANT_H
#define MATRIXD_DETERMINANT_H

#include "tools/SANSnumerics.h"
#include "tools/SANSException.h"
#include "MatrixD_Type.h"

namespace numpack
{
namespace DLA
{

namespace MatrixDdet
{
  template< class T >
  inline T
  det2x2(const MatrixDView<T>& M)
  {
    return M(0,0)*M(1,1) - M(0,1)*M(1,0);
  }

  template< class T >
  inline T
  det3x3(const MatrixDView<T>& M)
  {
    return M(0,0)*( M(1,1)*M(2,2) - M(1,2)*M(2,1) ) -
           M(0,1)*( M(1,0)*M(2,2) - M(1,2)*M(2,0) ) +
           M(0,2)*( M(1,0)*M(2,1) - M(1,1)*M(2,0) );
  }
}

//=============================================================================
  //Function to forward the call to the appropriate determinant calculation
  template< class T >
  inline T
  Det(const MatrixDView<T>& Matrix)
  {

    if ( Matrix.m() == 2 && Matrix.n() == 2 )
      return MatrixDdet::det2x2< T >( Matrix );
    if ( Matrix.m() == 3 && Matrix.n() == 3 )
      return MatrixDdet::det3x3< T >( Matrix );

    SANS_DEVELOPER_EXCEPTION(
            "Tried to compute determinant of %dx%d matrix.\nDeterminant only avilable for 2x2 and 3x3 matricies. ", Matrix.m(), Matrix.n());

    return T();
  }

  template<class T>
  inline T
  Det(const MatrixSymD<T>& M)
  {
    if (M.m()==2)
      return M(0,0)*M(1,1) - M(1,0)*M(1,0);
    else if (M.m()==3)
    {
      return M(0,0)*( M(1,1)*M(2,2) - M(1,2)*M(2,1) ) -
             M(0,1)*( M(1,0)*M(2,2) - M(1,2)*M(2,0) ) +
             M(0,2)*( M(1,0)*M(2,1) - M(1,1)*M(2,0) );
    }
    else if (M.m()==4)
    {
      T X1_1=M(0,0),X1_2=M(0,1),X1_3=M(0,2),X1_4=M(0,3);
      T X2_1=M(1,0),X2_2=M(1,1),X2_3=M(1,2),X2_4=M(1,3);
      T X3_1=M(2,0),X3_2=M(2,1),X3_3=M(2,2),X3_4=M(2,3);
      T X4_1=M(3,0),X4_2=M(3,1),X4_3=M(3,2),X4_4=M(3,3);
      return X1_1*X2_2*X3_3*X4_4 - X1_1*X2_2*X3_4*X4_3 - X1_1*X2_3*X3_2*X4_4
          + X1_1*X2_3*X3_4*X4_2 + X1_1*X2_4*X3_2*X4_3 - X1_1*X2_4*X3_3*X4_2
          - X1_2*X2_1*X3_3*X4_4 + X1_2*X2_1*X3_4*X4_3 + X1_2*X2_3*X3_1*X4_4
          - X1_2*X2_3*X3_4*X4_1 - X1_2*X2_4*X3_1*X4_3 + X1_2*X2_4*X3_3*X4_1
          + X1_3*X2_1*X3_2*X4_4 - X1_3*X2_1*X3_4*X4_2 - X1_3*X2_2*X3_1*X4_4
          + X1_3*X2_2*X3_4*X4_1 + X1_3*X2_4*X3_1*X4_2 - X1_3*X2_4*X3_2*X4_1
          - X1_4*X2_1*X3_2*X4_3 + X1_4*X2_1*X3_3*X4_2 + X1_4*X2_2*X3_1*X4_3
          - X1_4*X2_2*X3_3*X4_1 - X1_4*X2_3*X3_1*X4_2 + X1_4*X2_3*X3_2*X4_1;
    }
    else
      SANS_DEVELOPER_EXCEPTION("unsupported matrix size");
    return -1;
  }

} //namespace DLA
} //namespace numpack



#endif //MATRIXD_DETERMINANT_H
