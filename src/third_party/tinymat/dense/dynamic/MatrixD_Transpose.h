// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXD_TRANSPOSE_H
#define MATRIXD_TRANSPOSE_H

#include "tools/SANSnumerics.h"
#include "tools/SANSException.h"
#include "MatrixD_Type.h"
#include "tinymat/Transpose.h"
#include "tinymat/MatrixScatterAdd.h"

namespace tinymat 
{
namespace DLA
{

//Represents the transposed view of a matrix
template<class T>
class MatrixDTranspose : public MatrixDType< MatrixDTranspose<T>, true >,
                         public MatrixScatterAdd< typename TransposeTraits<T>::type >
{
  typedef MatrixScatterAdd< typename TransposeTraits<T>::type > BaseScatterAdd;
public:
  typedef typename TransposeTraits<T>::type node_type;

  //Transpose, so m = M.n()
  //Transpose, so n = M.m()
  explicit MatrixDTranspose( MatrixDView<T>& M ) : BaseScatterAdd(LA::eMatrixDTranspose, M.n(), M.m()), M_(M) {}
  explicit MatrixDTranspose( const MatrixDView<T>& M ) : BaseScatterAdd(LA::eMatrixDTranspose, M.n(), M.m()), M_(const_cast<MatrixDView<T>&>(M)) {}

  //Operator to access the matrix values
  inline       typename TransposeViewTraits<T>::type operator()(const int i, const int j)       { return Transpose(M_(j,i)); }
  inline const typename TransposeViewTraits<T>::type operator()(const int i, const int j) const { return Transpose(M_(j,i)); }

  // Lazy expression operations
  template<class Tr>
  inline void value(const Real sgn, MatrixDView<Tr>& res) const
  {
    const int m = this->m();
    const int n = this->n();

    SANS_ASSERT( M_.ID() != res.ID() );
    SANS_ASSERT( m == res.m() );
    SANS_ASSERT( n == res.n() );

    for (int i = 0; i < m; ++i)
      for (int j = 0; j < n; ++j)
        res(i,j) = sgn*Transpose(M_(j,i));
  }

  template<class Tr>
  inline void plus(const Real sgn, MatrixDView<Tr>& res) const
  {
    const int m = this->m();
    const int n = this->n();

    SANS_ASSERT( M_.ID() != res.ID() );
    SANS_ASSERT( m == res.m() );
    SANS_ASSERT( n == res.n() );

    for (int i = 0; i < m; ++i)
      for (int j = 0; j < n; ++j)
        res(i,j) += sgn*Transpose(M_(j,i));
  }

  int size() const { return M_.size(); }
  using BaseScatterAdd::m;
  using BaseScatterAdd::n;
  int stride() const { return M_.stride(); } //Memory stride does not change

  //A unique ID for the matrix so the matrix can be uniquely identified
  uintptr_t ID() const { return M_.ID(); }

  // Sparse matrix add to the matrix
  void scatterAdd( const MatrixDView< typename TransposeTraits<T>::type >& M, const int Map[], int nMap )
  {
    for (int i = 0; i < nMap; i++)
    {
      int iGlobal = Map[i];
      for (int j = 0; j < nMap; j++)
      {
        int jGlobal = Map[j];
        M_(jGlobal,iGlobal) += Transpose(M(i,j));
      }
    }
  }

  void scatterAdd( const MatrixDView< typename TransposeTraits<T>::type >& M, const int rowMap[], int nRow, const int colMap[], int nCol )
  {
    for (int i = 0; i < nRow; i++)
    {
      int iGlobal = rowMap[i];
      for (int j = 0; j < nCol; j++)
      {
        int jGlobal = colMap[j];
        M_(jGlobal,iGlobal) += Transpose(M(i,j));
      }
    }
  }

private:
  MatrixDView<T>& M_; //The transposed matrix
  using BaseScatterAdd::m_;
  using BaseScatterAdd::n_;
};

} //namespace DLA
} //namespace tinymat 



#endif //MATRIXD_TRANSPOSE_H
