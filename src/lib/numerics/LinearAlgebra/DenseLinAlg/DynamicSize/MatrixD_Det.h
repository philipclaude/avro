// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXD_DETERMINANT_H
#define MATRIXD_DETERMINANT_H

#include "tools/SANSnumerics.h"
#include "tools/SANSException.h"
#include "MatrixD_Type.h"

namespace SANS
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

} //namespace DLA
} //namespace SANS



#endif //MATRIXD_DETERMINANT_H
