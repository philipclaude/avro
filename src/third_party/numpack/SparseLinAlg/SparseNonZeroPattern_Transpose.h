// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef SPARSENONZEROPATTERN_TRANSPOSE_H
#define SPARSENONZEROPATTERN_TRANSPOSE_H

#include "tools/SANSnumerics.h" //Real

#include "numpack/Transpose.h"
#include "numpack/MatrixScatterAdd.h"

namespace numpack 
{

namespace SLA
{

//Forward declare
template<class TM>
class SparseNonZeroPattern;

//Represents the transpose of a non-zero pattern
template<class TM_>
class SparseNonZeroPattern_Transpose : public MatrixScatterAdd< typename TransposeTraits<TM_>::type >
{
  typedef MatrixScatterAdd< typename TransposeTraits<TM_>::type > BaseScatterAdd;

public:
  friend class SparseNonZeroPattern<TM_>;
  typedef typename TransposeTraits<TM_>::type TM;
  typedef TM Ttype;

  //Transpose, so m = M.n
  //Transpose, so n = M.m
  explicit SparseNonZeroPattern_Transpose( SparseNonZeroPattern<TM_>& M ) :
      BaseScatterAdd(LA::eSparseNonZeroPattern_Transpose, M.n(), M.m()), M_(M) {}

  SparseNonZeroPattern_Transpose( const SparseNonZeroPattern_Transpose<TM_>& nz ) : BaseScatterAdd(nz), M_(nz.M_) {}
  SparseNonZeroPattern_Transpose& operator=( const SparseNonZeroPattern_Transpose<TM_>& nz ) = delete;

  //Functions to add all the non-zero elements
  void scatterAdd( const DLA::MatrixDView< TM >&, const int Map[], const int nMap )
  {
    scatterAdd(Map,nMap);
  }
  void scatterAdd( const DLA::MatrixDView< TM >&, const int rowMap[], const int nRow, const int colMap[], const int nCol )
  {
    scatterAdd(rowMap,nRow,colMap,nCol);
  }

  void scatterAdd( const int Map[], const int nMap );
  void scatterAdd( const int rowMap[], const int nRow, const int colMap[], const int nCol );

  using BaseScatterAdd::m;
  using BaseScatterAdd::n;

protected:
  SparseNonZeroPattern<TM_>& M_; //The transposed non-zero pattern
  using BaseScatterAdd::m_;
  using BaseScatterAdd::n_;
};


//This is needed for MatrixD_Transpose
template<class TM>
const SparseNonZeroPattern_Transpose<TM>&
operator*(const Real, const SparseNonZeroPattern_Transpose<TM>& PatternT)
{
  return PatternT;
}

//Fill the non-zero pattern
template< class TM >
void
SparseNonZeroPattern_Transpose<TM>::scatterAdd( const int Map[], const int nMap )
{
  for (int i = 0; i < nMap; i++)
  {
    int iGlobal = Map[i];
    for (int j = 0; j < nMap; j++)
    {
      int jGlobal = Map[j];
      M_.add(jGlobal,iGlobal); // transposed i, j
    }
  }
}

//Fill the non-zero pattern
template< class TM >
void
SparseNonZeroPattern_Transpose<TM>::scatterAdd( const int rowMap[], const int nRow, const int colMap[], const int nCol )
{
  for (int i = 0; i < nRow; i++)
  {
    int iGlobal = rowMap[i];
    for (int j = 0; j < nCol; j++)
    {
      int jGlobal = colMap[j];
      M_.add(jGlobal,iGlobal); // transposed i, j
    }
  }
}


} //namespace SLA
} //namespace numpack 


#endif //SPARSENONZEROPATTERN_TRANSPOSE_H
