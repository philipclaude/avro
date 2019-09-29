// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

// No include blocker here. This file should only be included once in a cpp file for explicit instantiation

#if !defined(ILU0_INSTANTIATE) && !defined(SANS_HEADERCOMPILE_CHECK)
#error "This file should only be included in a cpp file for explicit instantiations"
#endif

#include "ILU0.h"
#include "numpack/dense/InverseLU.h"

namespace numpack 
{
namespace SLA
{

//-----------------------------------------------------------------------------
template< class Matrix_type >
void
ILU0< Matrix_type >::init()
{
  SystemNonZeroPattern nz(f_.matrixSize());
  if (transpose_)
    f_.jacobianTranspose(nz);
  else
    f_.jacobian(nz);

  A_ = new Matrix_type(nz);
}

//-----------------------------------------------------------------------------
template< class Matrix_type >
void
ILU0< Matrix_type >::factorize()
{
  // update the matrix
  *A_ = 0;
  if (transpose_)
    f_.jacobianTranspose(*A_);
  else
    f_.jacobian(*A_);

  const Matrix_type& A = *A_;

  const int nRow = A.m();

  //Reallocate a copy of A for the LU decomposition
  LU_.copy_from(A);

  const int* row_ptr = LU_.get_row_ptr();
  const int* col_ind = LU_.get_col_ind();

  // loop over all elements
  for (int i = 0; i < nRow; i++)
  {
    /*---------------*/
    /* Factorize row */
    /*---------------*/
    int diag_i = LU_.diagIndex(i)-row_ptr[i];
    for (int j = 0; j < diag_i; j++)
    {
      /*--------------------------*/
      /* Set Lij = Lij * inv(Ujj) */
      /*--------------------------*/

      // Note here that j is the row local non-zero indexing, not the full matrix indexing
      TM& Lij = LU_.sparseRow(i,j);
      int jcol = col_ind[row_ptr[i]+j];

      // Set Lij = Lij*inv(Ujj) since LU_ stores inv(Ujj) on the diagonal
      TM Temp = Lij*LU_.diag(jcol);
      Lij = Temp;

      /*------------------------*/
      /* Set Lik -= Lij * Ujk   */
      /*------------------------*/
      for (int k = row_ptr[i]+j+1; k < row_ptr[i+1]; k++)
      {
        // Get column index for Ujk
        int jk = LU_.sparseIndex(jcol, col_ind[k]);
        if (jk == -1) continue;

        TM& Lik = LU_[k];
        TM& Ujk = LU_[jk];

        Lik -= Lij * Ujk;
      }
    }

    /*-----------------------*/
    /* Invert Block diagonal */
    /*-----------------------*/
    TM& Uii = LU_.diag(i);

    Uii = DLA::InverseLU::Inverse(Uii);
  }
}



template< class Matrix_type >
LinearSolveStatus
ILU0< Matrix_type >::backsolve( const SparseVectorView_type& b, SparseVectorView_type& x ) const
{
  const int nRow = LU_.m();

  const int* row_ptr = LU_.get_row_ptr();
  const int* col_ind = LU_.get_col_ind();

  //Perform the forward sweep solving L x' = b, i.e. compute x' = b - L x'
  for ( int i = 0; i < nRow; i++)
  {
    TV sum( b[i] );

    int k = row_ptr[i];
    //Loop until we hit the diagonal
    while ( col_ind[k] < i )
    {
      sum -= LU_[k]*x[col_ind[k]];
      k++;
    }

    x[i] = sum;
  }

  //Back substitute by solving (D + U) x = x', i.e. compute x = D^-1 (x' - U x)
  for ( int i = nRow-1; i >= 0; i--)
  {
    int k = row_ptr[i+1]-1;

    TV sum( x[i] );

    //Loop until we hit the diagonal
    while ( col_ind[k] > i )
    {
      sum -= LU_[k]*x[col_ind[k]];
      k--;
    }

    //Note Uii has already been set to inv(Uii)
    x[i] = LU_.diag(i)*sum;
  }

  return LinearSolveStatus(true);
}



//-----------------------------------------------------------------------------
template< class Matrix_type >
void
ILU0< DLA::MatrixD<Matrix_type> >::init()
{
  SANS_DEVELOPER_EXCEPTION("Not implemented...");
}

//-----------------------------------------------------------------------------
template< class Matrix_type >
void
ILU0< DLA::MatrixD<Matrix_type> >::factorize()
{
  SANS_DEVELOPER_EXCEPTION("Not implemented...");
  /*
  Base_type::factorize(A);
  Dinv_.resize(A.m());

  for (int d = 0; d < A.m(); d++)
  {
    const Matrix_type& Add = A(d,d);
    const int nRow = Add.m();

    //Reallocate the inverse diagonal
    Dinv_[d].resize(Add);

    for ( int i = 0; i < nRow; i++)
      Dinv_[d][i] = DLA::InverseLU::Inverse(Add.diag(i));
  }
  */
}

//-----------------------------------------------------------------------------
template< class Matrix_type >
LinearSolveStatus
ILU0< DLA::MatrixD<Matrix_type> >::backsolve( const SparseVectorView_type& b, SparseVectorView_type& x ) const
{
  SANS_DEVELOPER_EXCEPTION("Not implemented...");
  /*
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
      TV sum( b(d,0)[i] );

      //Loop over the lower sparse matricies
      for ( int dL = 0; dL < d; dL++ )
      {
        const Matrix_type& AdL = A(d,dL);
        const int* row_ptrL = AdL.get_row_ptr();
        const int* col_indL = AdL.get_col_ind();

        for (int k = row_ptrL[i]; k < row_ptrL[i+1]; k++)
          sum -= AdL[k]*x(dL,0)[col_indL[k]];
      }

      //Account for the lower part of the diagonal sparse block
      int k = row_ptr[i];
      //Loop until we hit the diagonal
      while ( col_ind[k] < i )
      {
        sum -= Add[k]*x(d,0)[col_ind[k]];
        k++;
      }

      x(d,0)[i] = Dinv_[d][i]*sum;
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
          TV sum( -AdU[k]*x(dU,0)[col_indU[k]] );
          k--;

          while ( k >= row_ptrU[i] )
          {
            sum -= AdU[k]*x(dU,0)[col_indU[k]];
            k--;
          }

          x(d,0)[i] += Dinv_[d][i]*sum;
        }
      }

      //Account for the upper part of the diagonal sparse block
      int k = row_ptr[i+1]-1;
      if ( col_ind[k] > i )
      {
        TV sum( -Add[k]*x(d,0)[col_ind[k]] );
        k--;

        //Loop until we hit the diagonal
        while ( col_ind[k] > i )
        {
          sum -= Add[k]*x(d,0)[col_ind[k]];
          k--;
        }

        x(d,0)[i] += Dinv_[d][i]*sum;
      }

    }
  }
  */

  return LinearSolveStatus(true);
}

} //namespace SLA
} //namespace numpack 
