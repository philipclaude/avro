// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef SPARSEMATRIX_CRS_TRANSPOSE_H
#define SPARSEMATRIX_CRS_TRANSPOSE_H

#include "tools/SANSnumerics.h"
#include "SparseMatrix_CRS.h"
#include "SparseLinAlg_Type.h"
#include "numpack/Transpose.h"
#include "numpack/MatrixScatterAdd.h"

namespace numpack 
{
namespace SLA
{

//=============================================================================
//Represents the transpose of a matrix
template<class TM_>
class SparseMatrix_CRS_Transpose : public SparseLinAlgType< SparseMatrix_CRS_Transpose<TM_>, true >,
                                   public MatrixScatterAdd< typename TransposeTraits<TM_>::type >
{
protected:
  typedef TM_ Ttype_;
  typedef MatrixScatterAdd< typename TransposeTraits<TM_>::type > BaseScatterAdd;

public:
  typedef typename TransposeTraits<TM_>::type TM;

  //Transpose, so m = M.n
  //Transpose, so n = M.m
  // cppcheck-suppress noExplicitConstructor
  SparseMatrix_CRS_Transpose( SparseMatrix_CRS<TM_>& M ) : BaseScatterAdd(LA::eSparseMatrix_CRS_Transpose, M.n(), M.m()), M_(M) {}

  //Adds a matrix of values to the transposed sparse matrix
  void scatterAdd( const DLA::MatrixDView<TM>& M, const int Map[], const int nMap );
  void scatterAdd( const DLA::MatrixDView<TM>& M, const int rowMap[], const int nrow, const int colMap[], const int ncol );

  using BaseScatterAdd::m;
  using BaseScatterAdd::n;

private:
  SparseMatrix_CRS<Ttype_>& M_; //The transposed matrix
  using BaseScatterAdd::m_;
  using BaseScatterAdd::n_;
};

//=============================================================================
//Represents the transpose of a matrix
template<class TM_>
class SparseMatrix_CRS_Transpose< DLA::MatrixD<TM_> > : public SparseLinAlgType< SparseMatrix_CRS_Transpose< DLA::MatrixD<TM_> >, true >
{
protected:
  typedef DLA::MatrixD<TM_> Ttype_;

public:
  typedef typename TransposeTraits<TM_>::type TM;
  typedef DLA::MatrixD<TM> Ttype;

  // cppcheck-suppress noExplicitConstructor
  SparseMatrix_CRS_Transpose( SparseMatrix_CRS<Ttype_>& M ) : M_(M) {}

  //Adds a matrix of values to the transposed sparse matrix
  void scatterAdd( const DLA::MatrixDView<TM>& M, const int row, const int col );

  int m() const { return M_.n(); } //Transpose, so m = M_.n
  int n() const { return M_.m(); } //Transpose, so n = M_.m

private:
  SparseMatrix_CRS<Ttype_>& M_; //The transposed matrix
};


//=============================================================================
template <class TM_>
void
SparseMatrix_CRS_Transpose<TM_>::scatterAdd( const DLA::MatrixDView<TM>& M, const int Map[], const int nMap )
{
  int *row_ptr = M_.get_row_ptr();
  int *col_ind = M_.get_col_ind();

  //Add M transposed to the sparse matrix such that it is transposed
  for (int i = 0; i < nMap; i++)
  {
    int row = Map[i];
    const int col = row_ptr[row];
    for (int j = 0; j < nMap; j++)
    {
      int k = 0;
      while ( Map[j] != col_ind[col+k] ) k++;
      M_.sparseRow(row, k) += Transpose(M(j,i));
    }
  }
}


template <class TM_>
void
SparseMatrix_CRS_Transpose<TM_>::scatterAdd( const DLA::MatrixDView<TM>& M,
                                              const int rowMap[], const int nrow,
                                              const int colMap[], const int ncol )
{
  int *row_ptr = M_.get_row_ptr();
  int *col_ind = M_.get_col_ind();

  //Add M transposed to the sparse matrix such that it is transposed
  for (int i = 0; i < ncol; i++)
  {
    int row = colMap[i];
    const int col = row_ptr[row];
    for (int j = 0; j < nrow; j++)
    {
      int k = 0;
      while ( rowMap[j] != col_ind[col+k] ) k++;
      M_.sparseRow(row, k) += Transpose(M(j,i));
    }
  }
}


//=============================================================================
// Specialization for SparseMatrix_CRS< MatrixD<TM> >
template <class TM_>
void
SparseMatrix_CRS_Transpose<DLA::MatrixD<TM_>>::scatterAdd( const DLA::MatrixDView<TM>& M, const int row, const int col )
{
  int *row_ptr = M_.get_row_ptr();
  int *col_ind = M_.get_col_ind();

  //Add M transposed to the sparse matrix such that it is transposed
  int k = 0;
  const int col_start = row_ptr[col];
  while ( row != col_ind[col_start+k] ) k++;
  M_.sparseRow(col, k) += Transpose(M);
}

} //namespace SLA

} //namespace numpack 



#endif //SPARSEMATRIX_CRS_TRANSPOSE_H
