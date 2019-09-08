// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef SPARSENONZEROPATTERN_H
#define SPARSENONZEROPATTERN_H

#include "tools/SANSException.h"
#include "tools/SANSTraitsInitListAssign.h"
#include "LinearAlgebra/DenseLinAlg/DynamicSize/MatrixD_Type.h"
#include "SparseNonZeroPattern_Transpose.h"

#include "LinearAlgebra/MatrixScatterAdd.h"

#include <set>

namespace SANS
{
#ifdef __INTEL_COMPILER
namespace SLA
{

//Forward declaration
template< class TM >
class SparseNonZeroPattern;
}

//Create a specialization so to allow for the syntax
//   DLA::MatrixD< SparseNonZeroPattern<Real> >
//      pA = {{ {3,3}, {3,2} },
//            { {2,3}, {2,2} }};
//
// This is completely unnessary if the intel compiler could use templated initializer_list functions....
//
template<class TM>
struct initializer_list_assign< SLA::SparseNonZeroPattern< TM > >
{
  initializer_list_assign(SLA::SparseNonZeroPattern< TM >& val, const std::initializer_list<int>& s) { val = s; }
};
#endif

namespace SLA
{
// A class to represent the size of a sparse matrix
class SparseMatrixSize
{
public:
  SparseMatrixSize(const int m, const int n) : m_(m), n_(n) {}
  SparseMatrixSize() : m_(0), n_(0) {}
  // cppcheck-suppress noExplicitConstructor
  SparseMatrixSize( const std::initializer_list<int>& s )
  {
    SANS_ASSERT(s.size() == 2);
    m_ = *s.begin();
    n_ = *(s.begin()+1);
  }
  int m() const {return m_;}
  int n() const {return n_;}
  void resize(const int m, const int n)
  {
    m_ = m;
    n_ = n;
  }

  SparseMatrixSize& operator= ( const std::initializer_list<int>& s )
  {
    SANS_ASSERT(s.size() == 2);
    m_ = *s.begin();
    n_ = *(s.begin()+1);
    return *this;
  }

protected:
  int m_;
  int n_;
};

} //namespace SLA


#ifdef __INTEL_COMPILER
//
// This is completely unnessary if the intel compiler could use templated initializer_list functions....
//
struct initializer_list_assign< SLA::SparseMatrixSize >
{
  initializer_list_assign(SLA::SparseMatrixSize& val, const std::initializer_list<int>& s) { val = s; }
};
#endif

namespace SLA
{

template< class TM >
class SparseNonZeroPattern : public MatrixScatterAdd< TM >
{
public:
  SparseNonZeroPattern(const unsigned int m, const unsigned int n) : MatrixScatterAdd< TM >(LA::eSparseNonZeroPattern, m, n),
    columns_(new std::set<unsigned int>[m]) {}

  SparseNonZeroPattern(const SparseNonZeroPattern& Pattern) = delete;
  explicit SparseNonZeroPattern( const std::initializer_list<int>& s ) : SparseNonZeroPattern() { operator=(s); }
  explicit SparseNonZeroPattern( const SparseMatrixSize& spsize ) : SparseNonZeroPattern() { operator=(spsize); }
  explicit SparseNonZeroPattern( const SparseNonZeroPattern_Transpose<TM>& PatternT ) : MatrixScatterAdd< TM >(LA::eSparseNonZeroPattern, 0, 0),
    columns_(NULL)
  {
    operator=(PatternT);
  }

  friend class DLA::MatrixD< SparseNonZeroPattern<TM> >;
#ifdef __INTEL_COMPILER
  friend class initializer_list_assign< SparseNonZeroPattern<TM> >;
#endif
protected:
  SparseNonZeroPattern() : MatrixScatterAdd< TM >(LA::eSparseNonZeroPattern, 0, 0), columns_(NULL) {}
  SparseNonZeroPattern& operator=( const std::initializer_list<int>& s )
  {
    deallocate();
    SANS_ASSERT(s.size() == 2);
    m_ = *s.begin();
    n_ = *(s.begin()+1);
    SANS_ASSERT(m_ >= 0);
    SANS_ASSERT(n_ >= 0);

    columns_ = new std::set<unsigned int>[m_];
    return *this;
  }

  SparseNonZeroPattern& operator=( const SparseMatrixSize& spsize )
  { // construct from spsize object
    deallocate();
    m_ = spsize.m();
    n_ = spsize.n();
    SANS_ASSERT(m_ >= 0);
    SANS_ASSERT(n_ >= 0);

    columns_ = new std::set<unsigned int>[m_];
    return *this;
  }

public:

  ~SparseNonZeroPattern() { deallocate(); }

  //Allows the non-zero pattern to be assigned a transposed pattern
  void operator=(const SparseNonZeroPattern_Transpose<TM>& PatternT);

  //Allows the non-zero pattern to be copied from another SparseNonZeroPattern
  void operator=(const SparseNonZeroPattern<TM>& Pattern);

  //Adds a non-zero index to the matrix
  void add(const unsigned int row, const unsigned int col)
  {
    SANS_ASSERT_MSG( row < (unsigned int)m_ && col < (unsigned int)n_, "with m_=%d, n_=%d, row=%d, col=%d", m_, n_, row, col);
    columns_[row].insert(col);
  }

  void deallocate()
  {
    m_ = 0;
    n_ = 0;
    delete [] columns_; columns_ = NULL;
  }

  //Computes the number of non-zeros
  int nnz() const
  {
    int sum = 0;
    for ( int i = 0; i < m_; i++ )
      sum += columns_[i].size();
    return sum;
  }

  //Number of non-zero columns_ in a row
  std::size_t rowSize(const unsigned int row) const
  {
    return columns_[row].size();
  }

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

  // Functions to access a group of non-zero elements
  void scatterGet( DLA::MatrixDView< TM >& mtx, const int Map[], const int nMap ) { scatterGet(mtx,Map,nMap,Map,nMap); }
  void scatterGet( DLA::MatrixDView< TM >& mtx, const int rowMap[], const int nRow, const int colMap[], const int nCol );

  //Dense matrix like accessor
  TM operator()(const int i, const int j) const
  {
    SANS_ASSERT( i >= 0 && i < m_ );
    SANS_ASSERT( j >= 0 && j < n_ );
    if ( columns_[i].find(j) != columns_[i].end() )
      return TM(1);
    else
      return TM(0);
  }

  std::set<unsigned int>::const_iterator rowBegin(const unsigned int row) const { return columns_[row].cbegin(); }
  std::set<unsigned int>::const_iterator rowEnd(const unsigned int row) const { return columns_[row].cend(); }

  using MatrixScatterAdd< TM >::m;
  using MatrixScatterAdd< TM >::n;

protected:
  using MatrixScatterAdd< TM >::m_; //Number of rows
  using MatrixScatterAdd< TM >::n_; //Number of columns_
  std::set<unsigned int>* columns_;
};

//Fill the non-zero pattern transposed from PatternT
template< class TM >
void
SparseNonZeroPattern<TM>::operator=(const SparseNonZeroPattern_Transpose<TM>& PatternT)
{
  deallocate();
  m_ = PatternT.M_.n();
  n_ = PatternT.M_.m();
  columns_ = new std::set<unsigned int>[m_];
  for ( int row = 0; row < PatternT.M_.m(); row++ )
    for (auto col = PatternT.M_.rowBegin(row); col != PatternT.M_.rowEnd(row); col++)
      add(*col, row);
}

//Copy the non-zero pattern from Pattern
template< class TM >
void
SparseNonZeroPattern<TM>::operator=(const SparseNonZeroPattern<TM>& Pattern)
{
  deallocate();
  m_ = Pattern.m();
  n_ = Pattern.n();
  columns_ = new std::set<unsigned int>[m_];
  for ( int row = 0; row < m_; row++ )
    columns_[row] = Pattern.columns_[row];
}

//Fill the non-zero pattern
template< class TM >
void
SparseNonZeroPattern<TM>::scatterAdd( const int Map[], const int nMap )
{
  scatterAdd(Map, nMap, Map, nMap);
}

//Fill the non-zero pattern
template< class TM >
void
SparseNonZeroPattern<TM>::scatterAdd( const int rowMap[], const int nRow, const int colMap[], const int nCol )
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

// Functions to access a group of non-zero elements
template< class TM >
void
SparseNonZeroPattern<TM>::scatterGet( DLA::MatrixDView< TM >& mtx, const int rowMap[], const int nrow, const int colMap[], const int ncol )
{
  //Get a part of the values the sparse matrix and store them into M
  // rowMap/colMap contain dense-matrix-like row/column indexing as in the sparse matrix

  // check matrix and row/column map sizing
  SANS_ASSERT(mtx.m() == nrow);
  SANS_ASSERT(mtx.n() == ncol);

  for (int i = 0; i < nrow; i++) // loop over rows in M
  {
    int row = rowMap[i];
    for (int j = 0; j < ncol; j++) // loop over columns in M
    {
      int col = colMap[j];

      if ( columns_[row].find(col) != columns_[row].end() )
        mtx(i,j) = TM(1);
      else
        mtx(i,j) = TM(0);
    }
  }
}

//-----------------------------------------------------------------------------
template< class TM >
class SparseNonZeroPattern< DLA::MatrixD<TM> > : public MatrixScatterAdd< DLA::MatrixD<TM> >
{
public:
  SparseNonZeroPattern(const unsigned int m, const unsigned int n) : MatrixScatterAdd< DLA::MatrixD<TM> >(LA::eSparseNonZeroPattern, m, n),
    columns_(new std::set<unsigned int>[m]),
    block_m_(new int[m]),
    block_n_(new int[n])
  {
    for ( unsigned int i = 0; i < m; i++ ) block_m_[i] = 0;
    for ( unsigned int i = 0; i < n; i++ ) block_n_[i] = 0;
  }

  SparseNonZeroPattern(const SparseNonZeroPattern& Pattern) = delete;
  explicit SparseNonZeroPattern( const std::initializer_list<int>& s ) : SparseNonZeroPattern() { operator=(s); }
  explicit SparseNonZeroPattern( const SparseMatrixSize& spsize ) : SparseNonZeroPattern() { operator=(spsize); }
  explicit SparseNonZeroPattern( const SparseNonZeroPattern_Transpose<DLA::MatrixD<TM>>& PatternT ) : SparseNonZeroPattern()
  {
    operator=(PatternT);
  }

  friend class DLA::MatrixD< SparseNonZeroPattern< DLA::MatrixD<TM> > >;
#ifdef __INTEL_COMPILER
  friend class initializer_list_assign< SparseNonZeroPattern< DLA::MatrixD<TM> > >;
#endif
protected:
  SparseNonZeroPattern() : MatrixScatterAdd< DLA::MatrixD<TM> >(LA::eSparseNonZeroPattern, 0, 0),
  columns_(nullptr),
  block_m_(nullptr),
  block_n_(nullptr)
  {}
  SparseNonZeroPattern& operator=( const std::initializer_list<int>& s )
  {
    deallocate();
    SANS_ASSERT(s.size() == 2);
    m_ = *s.begin();
    n_ = *(s.begin()+1);
    SANS_ASSERT(m_ >= 0);
    SANS_ASSERT(n_ >= 0);

    columns_ = new std::set<unsigned int>[m_];

    block_m_ = new int[m_];
    block_n_ = new int[n_];
    for ( int i = 0; i < m_; i++ ) block_m_[i] = 0;
    for ( int i = 0; i < n_; i++ ) block_n_[i] = 0;

    return *this;
  }

  SparseNonZeroPattern& operator=( const SparseMatrixSize& spsize )
  { // construct from spsize object
    deallocate();
    m_ = spsize.m();
    n_ = spsize.n();
    SANS_ASSERT(m_ >= 0);
    SANS_ASSERT(n_ >= 0);

    columns_ = new std::set<unsigned int>[m_];

    block_m_ = new int[m_];
    block_n_ = new int[n_];
    for ( int i = 0; i < m_; i++ ) block_m_[i] = 0;
    for ( int i = 0; i < n_; i++ ) block_n_[i] = 0;

    return *this;
  }

public:

  ~SparseNonZeroPattern() { deallocate(); }

  //Allows the non-zero pattern to be assigned a transposed pattern
  void operator=(const SparseNonZeroPattern_Transpose<DLA::MatrixD<TM>>& PatternT);

  //Allows the non-zero pattern to be copied from another SparseNonZeroPattern
  void operator=(const SparseNonZeroPattern<DLA::MatrixD<TM>>& Pattern);

  //Adds a non-zero index to the matrix
  void add(const unsigned int row, const unsigned int col, const unsigned int bm, const unsigned int bn)
  {
    SANS_ASSERT_MSG( row < (unsigned int)m_ && col < (unsigned int)n_, "with m_=%d, n_=%d, row=%d, col=%d", m_, n_, row, col);
    columns_[row].insert(col);
    block_m_[row] = bm;
    block_n_[col] = bn;
  }

  //Number rows and columns in a block for the given sparse row/column
  int block_m(const unsigned int row) const { return block_m_[row]; }
  int block_n(const unsigned int col) const { return block_n_[col]; }

  void deallocate()
  {
    m_ = 0;
    n_ = 0;
    delete [] columns_; columns_ = NULL;
    delete [] block_m_; block_m_ = NULL;
    delete [] block_n_; block_n_ = NULL;
  }

  //Computes the number of non-zeros
  int nnz() const
  {
    int sum = 0;
    for ( int i = 0; i < m_; i++ )
      sum += columns_[i].size();
    return sum;
  }

  //Number of non-zero columns_ in a row
  std::size_t rowSize(const unsigned int row) const
  {
    return columns_[row].size();
  }

  //Functions to add all the non-zero elements
  void scatterAdd( const DLA::MatrixDView< TM >& M, const int row, const int col )
  {
    add(row,col,M.n(),M.m());
  }

  // Functions to access a group of non-zero elements
  void scatterGet( DLA::MatrixDView< TM >& M, const int row, const int col );

  //Dense matrix like accessor
  DLA::MatrixD<TM> operator()(const int i, const int j) const
  {
    SANS_ASSERT( i >= 0 && i < m_ );
    SANS_ASSERT( j >= 0 && j < n_ );
    SANS_ASSERT( block_m_[i] > 0 );
    SANS_ASSERT( block_n_[j] > 0 );
    if ( columns_[i].find(j) != columns_[i].end() )
      return DLA::MatrixD<TM>(block_m_[i],block_n_[j], TM(1));
    else
      return DLA::MatrixD<TM>(block_m_[i],block_n_[j], TM(0));
  }

  std::set<unsigned int>::const_iterator rowBegin(const unsigned int row) const { return columns_[row].cbegin(); }
  std::set<unsigned int>::const_iterator rowEnd(const unsigned int row) const { return columns_[row].cend(); }

  using MatrixScatterAdd< DLA::MatrixD<TM> >::m;
  using MatrixScatterAdd< DLA::MatrixD<TM> >::n;

protected:
  using MatrixScatterAdd< DLA::MatrixD<TM> >::m_; //Number of rows
  using MatrixScatterAdd< DLA::MatrixD<TM> >::n_; //Number of columns_
  std::set<unsigned int>* columns_;
  int* block_m_;
  int* block_n_;
};

//Fill the non-zero pattern transposed from PatternT
template< class TM >
void
SparseNonZeroPattern<DLA::MatrixD<TM>>::operator=(const SparseNonZeroPattern_Transpose<DLA::MatrixD<TM>>& PatternT)
{
  deallocate();
  m_ = PatternT.M_.n();
  n_ = PatternT.M_.m();
  columns_ = new std::set<unsigned int>[m_];
  block_m_ = new int[m_];
  block_n_ = new int[n_];
  for ( int row = 0; row < PatternT.M_.m(); row++ )
    for (auto col = PatternT.M_.rowBegin(row); col != PatternT.M_.rowEnd(row); col++)
      add(*col, row, PatternT.M_.block_n(*col), PatternT.M_.block_m(row));
}

//Copy the non-zero pattern from Pattern
template< class TM >
void
SparseNonZeroPattern<DLA::MatrixD<TM>>::operator=(const SparseNonZeroPattern<DLA::MatrixD<TM>>& Pattern)
{
  deallocate();
  m_ = Pattern.m();
  n_ = Pattern.n();
  columns_ = new std::set<unsigned int>[m_];
  block_m_ = new int[m_];
  block_n_ = new int[n_];
  for ( int row = 0; row < Pattern.m(); row++ )
    for (auto col = Pattern.rowBegin(row); col != Pattern.rowEnd(row); col++)
      add(row, *col, Pattern.block_m(row), Pattern.block_n(*col));
}

// Functions to access a group of non-zero elements
template< class TM >
void
SparseNonZeroPattern<DLA::MatrixD<TM>>::scatterGet( DLA::MatrixDView< TM >& M, const int row, const int col )
{
  //Get a part of the values the sparse matrix and store them into M

  // check matrix row/column sizing
  SANS_ASSERT(M.m() == block_m_[row]);
  SANS_ASSERT(M.n() == block_n_[col]);

  if ( columns_[row].find(col) != columns_[row].end() )
    M = TM(1);
  else
    M = TM(0);
}

} //namespace SLA
} //namespace SANS


#endif //SPARSENONZEROPATTERN_H
