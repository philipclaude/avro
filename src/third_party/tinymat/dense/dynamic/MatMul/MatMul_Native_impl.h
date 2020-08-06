// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

// No include blocker here. This file should only be included once in a cpp file for explicit instantiation

#if !defined(MATRIXD_MATMUL_NATIVE_INSTANTIATE) && !defined(SANS_HEADERCOMPILE_CHECK)
#error "This file should only be included in a cpp file for explicit instantiations"
#endif

#include <type_traits> //is_same

#include "tinymat/dense/static/MatrixS.h"
#include "tinymat/dense/static/VectorS.h"
#include "tinymat/dense/dynamic/MatrixD.h"
#include "MatMul_Native.h"
#include "tools/SANSException.h"

namespace tinymat 
{
namespace DLA
{

template<class TL, class TR, class T>
void MatMul_Native<TL, TR, T>::value(const MatrixDView<TL>& ML, const MatrixDView<TR>& MR, const Real sgn, MatrixDView<T>& res )
{
  res = 0;
  plus(ML, MR, sgn, res);
}

template<class TL, class TR, class T>
void MatMul_Native<TL, TR, T>::plus(const MatrixDView<TL>& ML, const MatrixDView<TR>& MR, const Real sgn, MatrixDView<T>& res )
{
  SANS_ASSERT(res.m() == ML.m());
  SANS_ASSERT(res.n() == MR.n());
  SANS_ASSERT(ML.n() == MR.m());

  if (std::is_same<T,TL>::value) SANS_ASSERT( res.ID() != ML.ID());
  if (std::is_same<T,TR>::value) SANS_ASSERT( res.ID() != MR.ID());

  const int mL = ML.m();
  const int nL = ML.n();

  //const int mR = MR.m(); because nL == mR
  const int nR = MR.n();

  const int Rstride = MR.stride();

  const int nRoundCacheL = (nL / CacheItemsL)*CacheItemsL;
  const int nRoundCacheR = (nR / CacheItemsR)*CacheItemsR;

  int i1, j1, k1, j2, k2;
  T *__restrict rC;

  for (i1 = 0; i1 < mL; ++i1)
  {
    T *C = &res(i1,0);
    const TL *L1 = &ML(i1,0);
    const TR *R1 = &MR(0,0);

//---------------
    for (j1 = 0; j1 < nRoundCacheL; j1 += CacheItemsL, R1 += CacheItemsL*Rstride)
    {
      const TL *__restrict rL1 = &L1[j1];
      for (k1 = 0; k1 < nRoundCacheR; k1 += CacheItemsR)
      {
        const TR *__restrict rR1 = &R1[k1];
        rC = &C[k1];
        for (j2 = 0; j2 < CacheItemsL; ++j2, rR1 += Rstride)
        {
          const TL& Lv = rL1[j2];
          for (k2 = 0; k2 < CacheItemsR; ++k2)
            rC[k2] += sgn*(Lv*rR1[k2]);
        }
      }

      if ( nR % CacheItemsR )
      {
        k1 = nRoundCacheR;
        const TR *__restrict rR1 = &R1[k1];
        rC = &C[k1];
        for (j2 = 0; j2 < CacheItemsL; ++j2, rR1 += Rstride)
        {
          const TL& Lv = rL1[j2];
          for (k2 = 0; k2 < nR % CacheItemsR; ++k2)
            rC[k2] += sgn*(Lv*rR1[k2]);
        }
      }
    }
//---------------

    if ( nL % CacheItemsL )
    {
      for (j1 = nRoundCacheL; j1 < nL; ++j1, R1 += Rstride)
      {
        const TL& Lv = L1[j1];
        for (k1 = 0; k1 < nRoundCacheR; k1 += CacheItemsR)
        {
          const TR *__restrict rR1 = &R1[k1];
          rC = &C[k1];
          for (k2 = 0; k2 < CacheItemsR; ++k2)
            rC[k2] += sgn*(Lv*rR1[k2]);
        }

        if ( nR % CacheItemsR )
        {
          k1 = nRoundCacheR;
          const TR *__restrict rR1 = &R1[k1];
          rC = &C[k1];
          for (k2 = 0; k2 < nR % CacheItemsR; ++k2)
            rC[k2] += sgn*(Lv*rR1[k2]);
        }
      }
    }
  }
}

} //namespace DLA
} //namespace tinymat 
