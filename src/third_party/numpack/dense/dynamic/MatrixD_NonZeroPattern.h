// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXD_NONZEROPATTERN_H
#define MATRIXD_NONZEROPATTERN_H

//===========================================================================//
// These are dummy non-zero pattern classes to keep the dense linear algebra
// consistent with sparse linear algebra

#include "tools/SANSnumerics.h"
#include "tools/SANSException.h"
#include "tools/SANSTraitsInitListAssign.h"
#include "MatrixD_Type.h"
#include "numpack/Transpose.h"
#include "numpack/MatrixScatterAdd.h"


namespace numpack 
{
#ifdef __INTEL_COMPILER
namespace DLA
{

//Forward declaration
template< class T >
class DenseNonZeroPattern;
}

//Create a specialization so to allow for the syntax
//   DLA::MatrixD< DenseNonZeroPattern<Real> >
//      pA = {{ {3,3}, {3,2} },
//            { {2,3}, {2,2} }};
//
// This is completely unnessary if the intel compiler could use templated initializer_list functions....
//
template<class T>
struct initializer_list_assign< DLA::DenseNonZeroPattern< T > >
{
  initializer_list_assign(DLA::DenseNonZeroPattern< T >& val, const std::initializer_list<int>& s) { val = s; }
};
#endif

namespace DLA
{
// A class to represent the size of a MatrixD similar to sparse
class DenseMatrixSize
{
public:
  DenseMatrixSize(const int m, const int n) : m_(m), n_(n) {}
  DenseMatrixSize() : m_(0), n_(0) {}
  // cppcheck-suppress noExplicitConstructor
  DenseMatrixSize( const std::initializer_list<int>& s )
  {
    SANS_ASSERT(s.size() == 2);
    m_ = *s.begin();
    n_ = *(s.begin()+1);
  }

  DenseMatrixSize& operator=( const std::initializer_list<int>& s )
  {
    SANS_ASSERT(s.size() == 2);
    m_ = *s.begin();
    n_ = *(s.begin()+1);
    return *this;
  }

  int m() const {return m_;}
  int n() const {return n_;}
  void resize(const int m, const int n)
  {
    m_ = m;
    n_ = n;
  }

protected:
  int m_;
  int n_;
};

// A class to represent the size of a VectorD
class DenseVectorSize
{
public:
  explicit DenseVectorSize(const int m) : m_(m) {}
  DenseVectorSize() : m_(0) {}
  // cppcheck-suppress noExplicitConstructor
  DenseVectorSize( const std::initializer_list<int>& s )
  {
    SANS_ASSERT(s.size() == 1);
    m_ = *s.begin();
  }
  void operator=(const int m) { m_ = m; };
  int m() const {return m_;}
  void resize(const int m)
  {
    m_ = m;
  }

protected:
  int m_;
};

} //namespace DLA


#ifdef __INTEL_COMPILER
//
// This is completely unnessary if the intel compiler could use templated initializer_list functions....
//
struct initializer_list_assign< DLA::DenseMatrixSize >
{
  initializer_list_assign(DLA::DenseMatrixSize& val, const std::initializer_list<int>& s) { val = s; }
};
#endif

namespace DLA
{

template< class T >
class DenseNonZeroPattern : public MatrixScatterAdd< T >
{
public:
  DenseNonZeroPattern(const unsigned int m, const unsigned int n) : MatrixScatterAdd<T>(LA::eDenseNonZeroPattern, m, n) {}

  DenseNonZeroPattern(const DenseNonZeroPattern& nz) : DenseNonZeroPattern<T>() { operator=(nz); }
  // cppcheck-suppress noExplicitConstructor
  DenseNonZeroPattern( const std::initializer_list<int>& s ) : DenseNonZeroPattern<T>() { operator=(s); }
  // cppcheck-suppress noExplicitConstructor
  DenseNonZeroPattern( const DenseMatrixSize& spsize ) : DenseNonZeroPattern<T>() { operator=(spsize); }

  friend class DLA::MatrixD< DenseNonZeroPattern<T> >;
#ifdef __INTEL_COMPILER
  friend class initializer_list_assign< DenseNonZeroPattern<T> >;
#endif
protected:
  DenseNonZeroPattern() : MatrixScatterAdd<T>(LA::eDenseNonZeroPattern, 0, 0) {}
  DenseNonZeroPattern& operator=( const std::initializer_list<int>& s )
  {
    SANS_ASSERT(s.size() == 2);
    m_ = *s.begin();
    n_ = *(s.begin()+1);
    SANS_ASSERT(m_ >= 0);
    SANS_ASSERT(n_ >= 0);

    return *this;
  }

  DenseNonZeroPattern& operator=( const DenseMatrixSize& spsize )
  {
    // construct from spsize object
    m_ = spsize.m();
    n_ = spsize.n();
    SANS_ASSERT(m_ >= 0);
    SANS_ASSERT(n_ >= 0);

    return *this;
  }

public:

  //Allows the non-zero pattern to be copied from another DenseNonZeroPattern
  DenseNonZeroPattern<T>& operator=(const DenseNonZeroPattern<T>& Pattern)
  {
    m_ = Pattern.m();
    n_ = Pattern.n();

    return *this;
  }

  //Adds a non-zero index to the matrix
  void add(const unsigned int row, const unsigned int col)
  {
    SANS_ASSERT_MSG( row < (unsigned int)m_ && col < (unsigned int)n_, "with m_=%d, n_=%d, row=%d, col=%d", m_, n_, row, col);
  }

  //Computes the number of non-zeros
  int nnz() const
  {
    return m_*n_;
  }

  //Number of non-zero columns_ in a row
  std::size_t rowSize(const unsigned int row) const
  {
    return n_;
  }

  //Functions to add all the non-zero elements
  void scatterAdd( const DLA::MatrixDView< T >&, const int Map[], const int nMap )
  {
    scatterAdd(Map,nMap);
  }
  void scatterAdd( const DLA::MatrixDView< T >&, const int rowMap[], const int nRow, const int colMap[], const int nCol )
  {
    scatterAdd(rowMap,nRow,colMap,nCol);
  }

  void scatterAdd( const int Map[], const int nMap );
  void scatterAdd( const int rowMap[], const int nRow, const int colMap[], const int nCol );

  T operator()(const int i, const int j) const
  {
    SANS_ASSERT( i >= 0 && i < m_ );
    SANS_ASSERT( j >= 0 && j < n_ );
    return T(1);
  }

  using MatrixScatterAdd< T >::m;
  using MatrixScatterAdd< T >::n;

protected:
  using MatrixScatterAdd< T >::m_; //Number of rows
  using MatrixScatterAdd< T >::n_; //Number of columns_
};

//Fill the non-zero pattern
template< class T >
void
DenseNonZeroPattern<T>::scatterAdd( const int Map[], const int nMap )
{
  for (int i = 0; i < nMap; i++)
  {
    int iGlobal = Map[i];
    for (int j = 0; j < nMap; j++)
    {
      int jGlobal = Map[j];
      add(iGlobal,jGlobal);
    }
  }
}

//Fill the non-zero pattern
template< class T >
void
DenseNonZeroPattern<T>::scatterAdd( const int rowMap[], const int nRow, const int colMap[], const int nCol )
{
  for (int i = 0; i < nRow; i++)
  {
    int iGlobal = rowMap[i];
    for (int j = 0; j < nCol; j++)
    {
      int jGlobal = colMap[j];
      add(iGlobal,jGlobal);
    }
  }
}

//=============================================================================
//This is needed for MatrixD_Transpose
template<class T>
const DenseNonZeroPattern<T>&
operator*(const Real, const DenseNonZeroPattern<T>& nz)
{
  return nz;
}

} //namespace DLA
} //namespace numpack 


#endif //MATRIXD_NONZEROPATTERN_H
