// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXD_SUB_H
#define MATRIXD_SUB_H

#include "tools/SANSnumerics.h"
#include "tools/SANSException.h"
#include "MatrixD_Type.h"

#include <vector>

namespace numpack 
{
namespace DLA
{

//=============================================================================
// Extracts a sub-matrix
template< class Ta, class Ts, class S >
void
subMatrixValue(const MatrixDView<Ta>& A, const std::vector<int>& rows, const std::vector<int>& cols,
               const S& sgn, MatrixDView<Ts>& Asub)
{
  SANS_ASSERT(Asub.m() == (int)rows.size());
  SANS_ASSERT(Asub.n() == (int)cols.size());

  for (std::size_t i = 0; i < rows.size(); i++)
    for (std::size_t j = 0; j < cols.size(); j++)
      Asub(i,j) = sgn*A(rows[i], cols[j]);
}

template< class Ta, class Ts, class S >
void
subMatrixPlus(const MatrixDView<Ta>& A, const std::vector<int>& rows, const std::vector<int>& cols,
              const S& sgn, MatrixDView<Ts>& Asub)
{
  SANS_ASSERT(Asub.m() == (int)rows.size());
  SANS_ASSERT(Asub.n() == (int)cols.size());

  for (std::size_t i = 0; i < rows.size(); i++)
    for (std::size_t j = 0; j < cols.size(); j++)
      Asub(i,j) += sgn*A(rows[i], cols[j]);
}

//=============================================================================
// Extracts a sub-vector
template< class Tv, class Ts, class S >
void
subVectorValue(const VectorDView<Tv>& V, const std::vector<int>& rows,
               const S& sgn, VectorDView<Ts>& Vsub)
{
  SANS_ASSERT(Vsub.m() == (int)rows.size());

  for (std::size_t i = 0; i < rows.size(); i++)
    Vsub[i] = sgn*V[rows[i]];
}

template< class Tv, class Ts, class S >
void
subVectorPlus(const VectorDView<Tv>& V, const std::vector<int>& rows,
              const S& sgn, VectorDView<Ts>& Vsub)
{
  SANS_ASSERT(Vsub.m() == (int)rows.size());

  for (std::size_t i = 0; i < rows.size(); i++)
    Vsub[i] += sgn*V[rows[i]];
}

} //namespace DLA
} //namespace numpack 



#endif //MATRIXD_DETERMINANT_H
