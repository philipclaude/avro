// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

//No include block because this file needs to be included twice

#ifdef TM_SPEC

namespace numpack 
{
namespace SLA
{

//=============================================================================
template<class TM>
Real
SparseMatrix_CRS< TM_SPEC >::operator=(const Real& v)
{
  for (int i = 0; i < nnz_; i++)
    (*this)[i] = v;

  return v;
}

//=============================================================================
template<class TM>
template<class TV1, class TV2>
void
SparseMatrix_CRS< TM_SPEC >::mulVec_value(const SparseVector<TV1>& x, const Real sgn, SparseVector<TV2>& b) const
{
  b = 0;
  mulVec_plus(x, sgn, b);
}

//=============================================================================
template<class TM>
template<class TV1, class TV2>
void
SparseMatrix_CRS< TM_SPEC >::mulVec_plus(const SparseVector<TV1>& x, const Real sgn, SparseVector<TV2>& b) const
{

  //If the matrix is empty, it is reated as a zero matrix and hence does nothing with the multiplication
  if (nnz_ == 0)
    return;

  SANS_ASSERT( x.m() == n_ );
  SANS_ASSERT( b.m() == m_ );

  for (int i = 0; i < m_; i++)
    for (int k = row_ptr_[i]; k < row_ptr_[i+1]; k++ )
      b[i] += sgn*((*this)[k]*x[col_ind_[k]]);
}

//=============================================================================
template <class TM>
bool
SparseMatrix_CRS< TM_SPEC >::isNonZero( const int i, const int j ) const
{
  SANS_ASSERT( i >= 0 && i < m_ );
  SANS_ASSERT( j >= 0 && j < n_ );
  if ( nnz_ == 0 ) return false; //Always false if there are no nonzero elements
  int k = 0;
  const int row = row_ptr_[i];
  const int rownnz = rowNonZero(i);
  while ( k < rownnz && j != col_ind_[row+k] ) k++;
  return k < rownnz; //Did we find the column or not?
}

//=============================================================================
template <class TM>
int
SparseMatrix_CRS< TM_SPEC >::sparseIndex( const int i, const int j ) const
{
  SANS_ASSERT_MSG( i >= 0 && i < m_, " i = %d, m_ = %d", i, m_ );
  SANS_ASSERT_MSG( j >= 0 && j < n_, " j = %d, n_ = %d", j, n_ );
  if ( nnz_ == 0 ) return -1; //Always -1 if there are no nonzero elements
  int k = 0;
  const int row = row_ptr_[i];
  const int rownnz = rowNonZero(i);
  while ( k < rownnz && j != col_ind_[row+k] ) k++;
  return k < rownnz ? row+k : -1; //Did we find the column or not?
}


} //namespace SLA
} //namespace numpack 

#endif
