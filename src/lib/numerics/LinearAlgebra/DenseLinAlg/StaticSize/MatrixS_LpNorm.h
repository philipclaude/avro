// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXS_LPNORM_H
#define MATRIXS_LPNORM_H

#include <cmath>  // pow

#include "tools/SANSnumerics.h"

#include "VectorS.h"
#include "Eigen.h"
#include "LinearAlgebra/DenseLinAlg/tools/dot.h"

namespace SANS
{
namespace DLA
{

//=============================================================================
template< int M, class T >
class LpNorm
{

public:
  template<class MatrixType>
  LpNorm(const MatrixType& A)
  {
    //Compute the eigen system with normalized eigen vectors
    EigenSystem( A, L, E );

    //Store the sqrt of the eigen values
    for (int d = 0; d < M; d++ )
      L[d] = sqrt(L[d]);
  }

  T Lp2(const VectorS<M,T>& dX, const int p)
  {
    //Compute the Lp norm square
    return pow(PowP(dX,p), 2./Real(p));
  }

  DLA::VectorS<M, T> dLp2dX(const VectorS<M,T>& dX, const int p)
  {
    //Compute the derivative of the Lp norm wit respect to dX
    DLA::VectorS<M, T> Ej;
    DLA::VectorS<M, T> dLpdX(0);
    const int pm = p-1;
    for ( int i = 0; i < M; i++)
      for ( int j = 0; j < M; j++)
      {
        Ej = E.col(j);
        dLpdX[i] += L[j]*Ej[i]*pow(L[j]*dot(Ej,dX), pm);
      }

    return 2*dLpdX*pow(PowP(dX,p), 2./Real(p) - 1);
  }

  T Lp(const VectorS<M,T>& dX, const int p)
  {
    //Compute the Lp norm
    return pow(PowP(dX,p), 1./Real(p));
  }

  DLA::VectorS<M, T> dLpdX(const VectorS<M,T>& dX, const int p)
  {
    //Compute the derivative of the Lp norm wit respect to dX
    DLA::VectorS<M, T> Ej;
    DLA::VectorS<M, T> dLpdX(0);
    const int pm = p-1;
    for ( int i = 0; i < M; i++)
      for ( int j = 0; j < M; j++)
      {
        Ej = E.col(j);
        dLpdX[i] += L[j]*Ej[i]*pow(L[j]*dot(Ej,dX), pm);
      }

    return dLpdX*pow(PowP(dX,p), 1./Real(p) - 1);
  }

protected:
  T PowP(const VectorS<M,T>& dX, const int p)
  {
    //Computes the pth power
    T powp = 0;
    for ( int i = 0; i < M; i++)
      powp += pow(L[i]*dot(E.col(i),dX), p);

    return powp;
  }

  DLA::MatrixS<M, M, T> E;
  DLA::VectorS<M, T> L;
};

} //namespace DLA
} //namespace SANS



#endif //MATRIXS_LPNORM_H
