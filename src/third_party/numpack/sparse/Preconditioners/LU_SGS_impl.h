// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

// No include blocker here. This file should only be included once in a cpp file for explicit instantiation

#if !defined(LU_SGS_INSTANTIATE) && !defined(SANS_HEADERCOMPILE_CHECK)
#error "This file should only be included in a cpp file for explicit instantiations"
#endif

#include "LU_SGS.h"
#include "numpack/dense/InverseLU.h"

namespace numpack 
{
namespace SLA
{

//-----------------------------------------------------------------------------
template< class Matrix_type >
void
LU_SGS< Matrix_type >::init()
{
  SystemNonZeroPattern nz(f_.matrixSize());
  if (transpose_)
    f_.jacobianTranspose(nz);
  else
    f_.jacobian(nz);

  A_ = new Matrix_type(nz);

  //Reallocate the inverse diagonal
  Dinv_.resize(*A_);
}

//-----------------------------------------------------------------------------
template< class Matrix_type >
void
LU_SGS< Matrix_type >::factorize()
{
  // update the matrix
  *A_ = 0;
  if (transpose_)
    f_.jacobianTranspose(*A_);
  else
    f_.jacobian(*A_);

  const int nRow = A_->m();

  // compute the inverse diagonal
  for ( int i = 0; i < nRow; i++)
    Dinv_[i] = DLA::InverseLU::Inverse(A_->diag(i));
}

//-----------------------------------------------------------------------------
template< class Matrix_type >
LinearSolveStatus
LU_SGS< Matrix_type >::backsolve( const SparseVectorView_type& b, SparseVectorView_type& x ) const
{
  const Matrix_type& A = *A_;
  const int nRow = A.m();

  const int* row_ptr = A.get_row_ptr();
  const int* col_ind = A.get_col_ind();

  //Perform the forward sweep solving (L + D) x' = b, i.e. compute x' = D^-1 (b - Lx')
  for ( int i = 0; i < nRow; i++)
  {
    TV sum( b[i] );

    int k = row_ptr[i];
    //Loop until we hit the diagonal
    while ( col_ind[k] < i )
    {
      sum -= A[k]*x[col_ind[k]];
      k++;
    }

    x[i] = Dinv_[i]*sum;
  }

  //Back substitute by solving (I + D^-1 U) x = x', i.e. compute x = x' - D^-1 U x
  for ( int i = nRow-2; i >= 0; i--)
  {

    int k = row_ptr[i+1]-1;

    if ( col_ind[k] <= i ) continue;

    TV sum( -A[k]*x[col_ind[k]] );
    k--;

    //Loop until we hit the diagonal
    while ( col_ind[k] > i )
    {
      sum -= A[k]*x[col_ind[k]];
      k--;
    }

    x[i] += Dinv_[i]*sum;
  }

  return LinearSolveStatus(true);
}


//-----------------------------------------------------------------------------
template< class Matrix_type >
void
LU_SGS< DLA::MatrixD<Matrix_type> >::init()
{
  SystemNonZeroPattern nz(f_.matrixSize());
  f_.jacobian(nz);

  A_ = new DLA::MatrixD<Matrix_type>(nz);

  const DLA::MatrixD<Matrix_type>& A = *A_;

  //Reallocate the inverse diagonal
  Dinv_.resize(A.m());

  for (int d = 0; d < A.m(); d++)
  {
    const Matrix_type& Add = A(d,d);

    //Reallocate the inverse diagonal
    Dinv_[d].resize(Add);
  }
}

//-----------------------------------------------------------------------------
template< class Matrix_type >
void
LU_SGS< DLA::MatrixD<Matrix_type> >::factorize()
{
  // update the matrix
  *A_ = 0;
  f_.jacobian(*A_);

  const DLA::MatrixD<Matrix_type>& A = *A_;

  for (int d = 0; d < A.m(); d++)
  {
    const Matrix_type& Add = A(d,d);
    const int nRow = Add.m();

    //Compute the inverse diagonal
    for ( int i = 0; i < nRow; i++)
      Dinv_[d][i] = DLA::InverseLU::Inverse(Add.diag(i));
  }
}



template< class Matrix_type >
LinearSolveStatus
LU_SGS< DLA::MatrixD<Matrix_type> >::backsolve( const SparseVectorView_type& b, SparseVectorView_type& x ) const
{
  const DLA::MatrixD<Matrix_type>& A = *A_;

  //Perform the forward sweep solving (L + D) x' = b, i.e. compute x' = D^-1 (b - Lx')
  for (int d = 0; d < A.m(); d++)
  {
    const Matrix_type& Add = A(d,d);

    const int nRow = Add.m();

    const int* row_ptr = Add.get_row_ptr();
    const int* col_ind = Add.get_col_ind();

    for ( int i = 0; i < nRow; i++)
    {
      TV sum( b[d][i] );

      //Loop over the lower sparse matricies
      for ( int dL = 0; dL < d; dL++ )
      {
        const Matrix_type& AdL = A(d,dL);
        const int* row_ptrL = AdL.get_row_ptr();
        const int* col_indL = AdL.get_col_ind();

        for (int k = row_ptrL[i]; k < row_ptrL[i+1]; k++)
          sum -= AdL[k]*x[dL][col_indL[k]];
      }

      //Account for the lower part of the diagonal sparse block
      int k = row_ptr[i];
      //Loop until we hit the diagonal
      while ( col_ind[k] < i )
      {
        sum -= Add[k]*x[d][col_ind[k]];
        k++;
      }

      x[d][i] = Dinv_[d][i]*sum;
    }
  }


  //Back substitute by solving (I + D^-1 U) x = x', i.e. compute x = x' - D^-1 U x
  for (int d = A.m()-1; d >= 0; d--)
  {
    const Matrix_type& Add = A(d,d);

    const int nRow = Add.m();

    const int* row_ptr = Add.get_row_ptr();
    const int* col_ind = Add.get_col_ind();

    for ( int i = nRow-1; i >= 0; i--)
    {

      //Loop over the upper sparse matricies
      for ( int dU = A.n()-1; dU > d; dU-- )
      {
        const Matrix_type& AdU = A(d,dU);
        const int* row_ptrU = AdU.get_row_ptr();
        const int* col_indU = AdU.get_col_ind();

        int k = row_ptrU[i+1]-1;

        if ( k >= 0 )
        {
          TV sum( -AdU[k]*x[dU][col_indU[k]] );
          k--;

          while ( k >= row_ptrU[i] )
          {
            sum -= AdU[k]*x[dU][col_indU[k]];
            k--;
          }

          x[d][i] += Dinv_[d][i]*sum;
        }
      }

      //Account for the upper part of the diagonal sparse block
      int k = row_ptr[i+1]-1;
      if ( col_ind[k] > i )
      {
        TV sum( -Add[k]*x[d][col_ind[k]] );
        k--;

        //Loop until we hit the diagonal
        while ( col_ind[k] > i )
        {
          sum -= Add[k]*x[d][col_ind[k]];
          k--;
        }

        x[d][i] += Dinv_[d][i]*sum;
      }

    }
  }

  return LinearSolveStatus(true);
}

} //namespace SLA
} //namespace numpack 
