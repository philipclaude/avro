// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXS_DETERMINANT_H
#define MATRIXS_DETERMINANT_H

#include "tools/SANSException.h"
#include "MatrixS_Type.h"

namespace SANS
{
namespace DLA
{
  //Determinant calculations for 1x1, 2x2 and 3x3 matrices. Anything else will create a compiler error.
  template<class T>
  inline T
  Det(const MatrixS<1, 1, T>& M)
  {
    return M(0,0);
  }

  template<class T>
  inline T
  Det(const MatrixSymS<1, T>& M)
  {
    return M(0,0);
  }

  template<class T>
  inline T
  Det(const MatrixS<2, 2, T>& M)
  {
    return M(0,0)*M(1,1) - M(0,1)*M(1,0);
  }

  template<class T>
  inline T
  Det(const MatrixSymS<2, T>& M)
  {
    return M(0,0)*M(1,1) - M(1,0)*M(1,0);
  }

  template<class T>
  inline T
  Det(const MatrixS<3, 3, T>& M)
  {
    return M(0,0)*( M(1,1)*M(2,2) - M(1,2)*M(2,1) ) -
           M(0,1)*( M(1,0)*M(2,2) - M(1,2)*M(2,0) ) +
           M(0,2)*( M(1,0)*M(2,1) - M(1,1)*M(2,0) );
  }

  //This can simplified to take advantage of the symmetry
  template<class T>
  inline T
  Det(const MatrixSymS<3, T>& M)
  {
    return M(0,0)*( M(1,1)*M(2,2) - M(1,2)*M(2,1) ) -
           M(0,1)*( M(1,0)*M(2,2) - M(1,2)*M(2,0) ) +
           M(0,2)*( M(1,0)*M(2,1) - M(1,1)*M(2,0) );
  }

  template<class T>
  inline T
  Det(const MatrixS<4, 4, T>& X)
  {
    T X1_1=X(0,0),X1_2=X(0,1),X1_3=X(0,2),X1_4=X(0,3);
    T X2_1=X(1,0),X2_2=X(1,1),X2_3=X(1,2),X2_4=X(1,3);
    T X3_1=X(2,0),X3_2=X(2,1),X3_3=X(2,2),X3_4=X(2,3);
    T X4_1=X(3,0),X4_2=X(3,1),X4_3=X(3,2),X4_4=X(3,3);
    return X1_1*X2_2*X3_3*X4_4 - X1_1*X2_2*X3_4*X4_3 - X1_1*X2_3*X3_2*X4_4
        + X1_1*X2_3*X3_4*X4_2 + X1_1*X2_4*X3_2*X4_3 - X1_1*X2_4*X3_3*X4_2
        - X1_2*X2_1*X3_3*X4_4 + X1_2*X2_1*X3_4*X4_3 + X1_2*X2_3*X3_1*X4_4
        - X1_2*X2_3*X3_4*X4_1 - X1_2*X2_4*X3_1*X4_3 + X1_2*X2_4*X3_3*X4_1
        + X1_3*X2_1*X3_2*X4_4 - X1_3*X2_1*X3_4*X4_2 - X1_3*X2_2*X3_1*X4_4
        + X1_3*X2_2*X3_4*X4_1 + X1_3*X2_4*X3_1*X4_2 - X1_3*X2_4*X3_2*X4_1
        - X1_4*X2_1*X3_2*X4_3 + X1_4*X2_1*X3_3*X4_2 + X1_4*X2_2*X3_1*X4_3
        - X1_4*X2_2*X3_3*X4_1 - X1_4*X2_3*X3_1*X4_2 + X1_4*X2_3*X3_2*X4_1;
  }

  template<class T>
  inline T
  Det(const MatrixSymS<4, T>& X)
  {
    T X1_1=X(0,0),X1_2=X(0,1),X1_3=X(0,2),X1_4=X(0,3);
    T X2_1=X(1,0),X2_2=X(1,1),X2_3=X(1,2),X2_4=X(1,3);
    T X3_1=X(2,0),X3_2=X(2,1),X3_3=X(2,2),X3_4=X(2,3);
    T X4_1=X(3,0),X4_2=X(3,1),X4_3=X(3,2),X4_4=X(3,3);
    return X1_1*X2_2*X3_3*X4_4 - X1_1*X2_2*X3_4*X4_3 - X1_1*X2_3*X3_2*X4_4
        + X1_1*X2_3*X3_4*X4_2 + X1_1*X2_4*X3_2*X4_3 - X1_1*X2_4*X3_3*X4_2
        - X1_2*X2_1*X3_3*X4_4 + X1_2*X2_1*X3_4*X4_3 + X1_2*X2_3*X3_1*X4_4
        - X1_2*X2_3*X3_4*X4_1 - X1_2*X2_4*X3_1*X4_3 + X1_2*X2_4*X3_3*X4_1
        + X1_3*X2_1*X3_2*X4_4 - X1_3*X2_1*X3_4*X4_2 - X1_3*X2_2*X3_1*X4_4
        + X1_3*X2_2*X3_4*X4_1 + X1_3*X2_4*X3_1*X4_2 - X1_3*X2_4*X3_2*X4_1
        - X1_4*X2_1*X3_2*X4_3 + X1_4*X2_1*X3_3*X4_2 + X1_4*X2_2*X3_1*X4_3
        - X1_4*X2_2*X3_3*X4_1 - X1_4*X2_3*X3_1*X4_2 + X1_4*X2_3*X3_2*X4_1;
  }


} //namespace DLA
} //namespace SANS



#endif //MATRIXS_DETERMINANT_H
