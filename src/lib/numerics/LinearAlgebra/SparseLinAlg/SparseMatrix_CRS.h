// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef SPARSEMATRIX_CRS_H
#define SPARSEMATRIX_CRS_H

#include "tools/SANSnumerics.h"

#include "SparseLinAlg_Type.h"
#include "SparseVector.h"
#include "SparseNonZeroPattern.h"

#include "LinearAlgebra/MatrixScatterAdd.h"

namespace SANS
{
namespace SLA
{

template<class TM>
class SparseMatrix_CRS : public SparseLinAlgType< SparseMatrix_CRS<TM>, true >,
                         public MatrixScatterAdd< TM >
{
public:
  typedef TM Ttype;
  typedef SparseNonZeroPattern<TM> NonZeroPattern;

  friend class DLA::MatrixD< SparseMatrix_CRS<TM> >;
protected:
  //Used in DLA::MatrixD to allocate an array of SparseMatrix_CRS
  SparseMatrix_CRS& operator=(const SparseNonZeroPattern<TM>& Pattern) { resize(Pattern); return *this; }

public:
  SparseMatrix_CRS() : MatrixScatterAdd< TM >(LA::eSparseMatrix_CRS, 0 ,0),
    val_(NULL), col_ind_(NULL), row_ptr_(NULL), nnz_(0), diag_ind_(0) {}

  explicit SparseMatrix_CRS(const SparseNonZeroPattern<TM>& Pattern) : SparseMatrix_CRS() { resize(Pattern); }

  //Copy constructor is intentionally removed to minimize memory usage.
  SparseMatrix_CRS(const SparseMatrix_CRS&) = delete;
  SparseMatrix_CRS& operator=(const SparseMatrix_CRS&) = delete;

  ~SparseMatrix_CRS() { deallocate(); }

  //Release all the memory and reset the pointers
  void deallocate()
  {
    m_ = 0; n_ = 0; nnz_ = 0;
    delete [] val_;      val_ = NULL;
    delete [] col_ind_;  col_ind_ = NULL;
    delete [] row_ptr_;  row_ptr_ = NULL;
    delete [] diag_ind_; diag_ind_ = NULL;
  }

  //Assignment operators
  Real operator=(const Real& v);

  // Reallocates memory and non-zero pattern
  void resize(const SparseNonZeroPattern<TM>& Pattern);

  // Copies A into this
  void copy_from( const SparseMatrix_CRS& A);

  //Access operator where the 1st index is the row and
  //the 2nd index is the count based on the number of non-zero columns in that row
        TM& sparseRow(const unsigned int row, const unsigned int col)       { return val_[row_ptr_[row] + col]; }
  const TM& sparseRow(const unsigned int row, const unsigned int col) const { return val_[row_ptr_[row] + col]; }

  //Returns the (dense) column location for a given row and non-zero column
  int columnIndex(const unsigned int row, const unsigned int col) const { return col_ind_[row_ptr_[row] + col]; }

  //Operator to directly access the non-zero values
        TM& operator[](const unsigned int k)       { return val_[k]; }
  const TM& operator[](const unsigned int k) const { return val_[k]; }

  //Dense matrix like accessor
        TM& slowAccess( const int i, const int j );
  const TM& slowAccess( const int i, const int j ) const;
  const TM& operator()( const int i, const int j ) const;

  //Simple access operator for the diagonal entries
        TM& diag(const unsigned int row)       { return val_[diag_ind_[row]]; }
  const TM& diag(const unsigned int row) const { return val_[diag_ind_[row]]; }

  //Gives the index for the diagonal element
  int diagIndex(unsigned int row) const { return diag_ind_[row]; }

  //Determines if an entry in the sparse matrix is non-zero
  bool isNonZero( const int i, const int j ) const;

  //Returns the row index if i,j is non-zero, else returns -1
  int sparseIndex( const int i, const int j ) const;

  //Adds a matrix of values to the sparse matrix
  void scatterAdd( const DLA::MatrixDView< TM >& M, const int Map[], const int nMap );
  void scatterAdd( const DLA::MatrixDView< TM >& M, const int rowMap[], const int nrow, const int colMap[], const int ncol );

  //Adds a matrix of values to the sparse matrix
  void scatterGet( DLA::MatrixDView< TM >& M, const int Map[], const int nMap ) const { scatterGet(M, Map, nMap, Map, nMap); }
  void scatterGet( DLA::MatrixDView< TM >& M, const int rowMap[], const int nrow, const int colMap[], const int ncol ) const;

  template<class TV1, class TV2>
  void mulVec_value(const SparseVector<TV1>& x, const Real sgn, SparseVector<TV2>& b) const;

  template<class TV1, class TV2>
  void mulVec_plus(const SparseVector<TV1>& x, const Real sgn, SparseVector<TV2>& b) const;

  using MatrixScatterAdd< TM >::m;
  using MatrixScatterAdd< TM >::n;
  int rowNonZero(const int i) const { return row_ptr_[i+1]-row_ptr_[i]; }
  int getNumNonZero()         const { return nnz_; }
  int* get_row_ptr()          const { SANS_ASSERT( row_ptr_ != NULL ); return row_ptr_; }
  int* get_col_ind()          const { SANS_ASSERT( col_ind_ != NULL ); return col_ind_; }
  TM* get_values()            const { SANS_ASSERT( val_ != NULL );     return val_; }

protected:
  using MatrixScatterAdd< TM >::m_;
  using MatrixScatterAdd< TM >::n_;
  TM *val_;      //Non zero values
  int *col_ind_; //Sparse non-zero to dense matrix index map (typical CRS notation)
  int *row_ptr_; //Index into val_ for start of non-zero row (typical CRS notation)
  int nnz_;      //Total number of non-zero elements

  int *diag_ind_; //The sparse row index for the diagonal elements
};

template<class TM>
void
SparseMatrix_CRS<TM>::resize(const SparseNonZeroPattern<TM>& Pattern)
{
  //Release the memory
  deallocate();

  //Initialize indexes for CRS storage
  m_ = Pattern.m();
  n_ = Pattern.n();
  nnz_ = Pattern.nnz();

  //No point in doing anything if there are no non-zero elements
  if (nnz_ == 0) return;

  //nnz_ could be negative if we exceed the maximum int value
  SANS_ASSERT( nnz_ > 0 );

  val_ = new TM[nnz_];
  col_ind_ = new int[nnz_];
  row_ptr_ = new int[m_+1];
  diag_ind_ = new int[m_];

  for (int i = 0; i < nnz_; i++)
    val_[i] = 0;

  row_ptr_[0] = 0;
  for ( int i = 0; i < m_; i++ )
  {
    row_ptr_[i+1] = row_ptr_[i] + Pattern.rowSize(i);
    int j = 0;
    for (auto col = Pattern.rowBegin(i); col != Pattern.rowEnd(i); col++, j++)
    {
      col_ind_[row_ptr_[i] + j] = *col;
      if ( i == (int)*col ) diag_ind_[i] = row_ptr_[i] + j;
    }
  }

  return;
}

//=============================================================================
template <class TM>
void
SparseMatrix_CRS< TM >::copy_from( const SparseMatrix_CRS< TM >& A )
{
  //Release the memory
  deallocate();

  //Initialize indexes for CRS storage
  m_ = A.m_;
  n_ = A.n_;
  nnz_ = A.nnz_;

  //No point in doing anything if there are no non-zero elements
  if (nnz_ == 0) return;

  //nnz_ could be negative if we exceed the maximum int value
  SANS_ASSERT( nnz_ > 0 );

  const int m = m_;
  const int nnz = nnz_;

  // Allocate local arrays and copy over everything from A
  val_ = new TM[nnz];
  col_ind_ = new int[nnz];
  row_ptr_ = new int[m+1];
  diag_ind_ = new int[m];

  for (int i = 0; i < nnz; i++)
  {
    val_[i] = A.val_[i];
    col_ind_[i] = A.col_ind_[i];
  }

  for (int i = 0; i < m; i++)
  {
    diag_ind_[i] = A.diag_ind_[i];
    row_ptr_[i] = A.row_ptr_[i];
  }
  row_ptr_[m_] = A.row_ptr_[m_];
}

// Add a square matrix of values
template <class TM>
void
SparseMatrix_CRS< TM >::scatterAdd( const DLA::MatrixDView< TM >& M, const int Map[], const int nMap )
{
  //Add the values of M to the sparse matrix
  scatterAdd(M, Map, nMap, Map, nMap);
}

// Add a general (rectangular) matrix of values
template <class TM>
void
SparseMatrix_CRS< TM >::scatterAdd( const DLA::MatrixDView< TM >& M,
                                    const int rowMap[], const int nrow,
                                    const int colMap[], const int ncol )
{
  //Add the values of M to the sparse matrix
  // rowMap/colMap contain dense-matrix-like row/column indexing as in the sparse matrix

  // check matrix and row/column map sizing
  SANS_ASSERT(nrow <= M.m());
  SANS_ASSERT(ncol == M.n());

  for (int i = 0; i < nrow; i++) // loop over rows in M
  {
    int row = rowMap[i];
    const int col = row_ptr_[row]; // global CRS index for the first nonzero entry in row
    for (int j = 0; j < ncol; j++) // loop over columns in M
    {
      // k_th nonzero entry in row corresponds to (rowMap[i], colMap[j]) in the sparse matrix
      int k = 0;
      while ( colMap[j] != col_ind_[col+k] ) k++;
      sparseRow(row, k) += M(i,j);
    }
  }
}

// Get a general (rectangular) matrix of values
template <class TM>
void
SparseMatrix_CRS< TM >::scatterGet( DLA::MatrixDView< TM >& M, const int rowMap[], const int nrow, const int colMap[], const int ncol ) const
{
  //Add the values of M to the sparse matrix
  // rowMap/colMap contain dense-matrix-like row/column indexing as in the sparse matrix

  // check matrix and row/column map sizing
  SANS_ASSERT(M.m() == nrow);
  SANS_ASSERT(M.n() == ncol);

  for (int i = 0; i < nrow; i++) // loop over rows in M
  {
    const int row = rowMap[i];
    SANS_ASSERT( row >= 0 && row < m_ );

    const int index_CRS = row_ptr_[row]; // global CRS index for the first nonzero entry in row
    const int rownnz = rowNonZero(row);

    for (int j = 0; j < ncol; j++) // loop over columns in M
    {
      const int col = colMap[j];
      SANS_ASSERT( col >= 0 && col < n_ );

      int k = 0;
      while ( k < rownnz && col != col_ind_[index_CRS+k] ) k++;

      if ( k < rownnz )
        M(i,j) = val_[index_CRS + k];
      else
        M(i,j) = TM(0);
    }
  }
}


template <class TM>
TM&
SparseMatrix_CRS< TM >::slowAccess( const int i, const int j )
{
  SANS_ASSERT( i >= 0 && i < m_ );
  SANS_ASSERT( j >= 0 && j < n_ );
  int k = 0;
  const int row = row_ptr_[i];
  const int rownnz = rowNonZero(i);
  while ( k < rownnz && j != col_ind_[row+k] ) k++;
  SANS_ASSERT( k < rownnz );
  return val_[row + k];
}

template <class TM>
const TM&
SparseMatrix_CRS< TM >::slowAccess( const int i, const int j ) const
{
  return const_cast<SparseMatrix_CRS< TM >*>(this)->slowAccess(i,j);
}

template <class TM>
const TM&
SparseMatrix_CRS< TM >::operator()( const int i, const int j ) const
{
  return const_cast<SparseMatrix_CRS< TM >*>(this)->slowAccess(i,j);
}

} //namespace SLA
} //namespace SANS

#define TM_SPEC TM
#include "SparseMatrix_CRS_impl.h"
#undef TM_SPEC


namespace SANS
{
namespace SLA
{
//=============================================================================
//Specialization for a Block Sparse Matrix
template<class TM>
class SparseMatrix_CRS< DLA::MatrixD<TM> > : public SparseLinAlgType< SparseMatrix_CRS< DLA::MatrixD<TM> >, true >,
                                             public MatrixScatterAdd< DLA::MatrixD<TM> >
{
public:
  typedef DLA::MatrixD<TM> Ttype;
  typedef DLA::MatrixDView<TM> MatrixView_type;
  typedef SparseNonZeroPattern<Ttype> NonZeroPattern;

  friend class DLA::MatrixD< SparseMatrix_CRS< DLA::MatrixD<TM> > >;
protected:
  //Used in DLA::MatrixD to allocate an array of SparseMatrix_CRS
  SparseMatrix_CRS& operator=(const SparseNonZeroPattern< DLA::MatrixD<TM> >& Pattern) { resize(Pattern); return *this; }

public:
  SparseMatrix_CRS() : MatrixScatterAdd< DLA::MatrixD<TM> >(LA::eSparseMatrix_CRS, 0 ,0),
    block_m_( NULL ), block_n_( NULL ), block_i_( NULL ),
    val_(NULL),
    row_ind_(NULL), col_ind_(NULL), row_ptr_(NULL),
    nnz_(0),
    diag_ind_(0)
  {}
  explicit SparseMatrix_CRS(const SparseNonZeroPattern< DLA::MatrixD<TM> >& Pattern) : SparseMatrix_CRS() { resize(Pattern); }

  SparseMatrix_CRS(const SparseMatrix_CRS&) = delete;

  // The destructor is intentionally not virtual. Make it virtual will just slow down the code needlessly
  ~SparseMatrix_CRS() { deallocate(); }

  //Release all the memory and rest the pointers
  void deallocate()
  {
    m_ = 0; n_ = 0; nnz_ = 0;
    delete [] val_;     val_ = NULL;
    delete [] row_ind_; row_ind_ = NULL;
    delete [] col_ind_; col_ind_ = NULL;
    delete [] row_ptr_; row_ptr_ = NULL;
    delete [] block_m_; block_m_ = NULL;
    delete [] block_n_; block_n_ = NULL;
    delete [] block_i_; block_i_ = NULL;
    delete [] diag_ind_; diag_ind_ = NULL;
  }

  //Assignment operators
  Real operator=(const Real& v);

  // Reallocates memory and non-zero pattern
  void resize(const SparseNonZeroPattern< DLA::MatrixD<TM> >& Pattern);

  // Copies A into this
  void copy_from( const SparseMatrix_CRS& A);

  //Access operator where the 1st index is the row and
  //the 2nd index is the count based on the number of non-zero columns in that row
        MatrixView_type sparseRow(unsigned int row, unsigned int col)       { return operator[]( row_ptr_[row] + col ); }
  const MatrixView_type sparseRow(unsigned int row, unsigned int col) const { return operator[]( row_ptr_[row] + col ); }

  //Operator to directly access the non-zero values
  MatrixView_type operator[](unsigned int k)
  {
    return MatrixView_type( val_ + block_i_[k], block_m_[row_ind_[k]], block_n_[col_ind_[k]]);
  }
  const MatrixView_type operator[](unsigned int k) const
  {
    return MatrixView_type( val_ + block_i_[k], block_m_[row_ind_[k]], block_n_[col_ind_[k]]);
  }

  //Dense matrix like accessor
        MatrixView_type slowAccess( const int i, const int j );
  const MatrixView_type slowAccess( const int i, const int j ) const;
  const MatrixView_type operator()( const int i, const int j ) const;

  //Simple access operator for the diagonal entries
        MatrixView_type diag(unsigned int row)       { return operator[](diag_ind_[row]); }
  const MatrixView_type diag(unsigned int row) const { return operator[](diag_ind_[row]); }

  //Gives the index for the diagonal element
  int diagIndex(unsigned int row) const { return diag_ind_[row]; }

  //Adds a matrix of values to sparse matrix
  void scatterAdd( const DLA::MatrixDView< TM >& M, const int row, const int col );

  //Determines if an entry in the sparse matrix is non-zero
  bool isNonZero( const int i, const int j ) const;

  //Returns the row index if i,j is non-zero, else returns -1
  int sparseIndex( const int i, const int j ) const;

  template<class TV1, class TV2>
  void mulVec_value(const SparseVector<TV1>& x, const Real sgn, SparseVector<TV2>& b) const;

  template<class TV1, class TV2>
  void mulVec_plus(const SparseVector<TV1>& x, const Real sgn, SparseVector<TV2>& b) const;

  int block_m( unsigned int row ) const { return block_m_[row]; }
  int block_n( unsigned int col ) const { return block_n_[col]; }
  int block_i( unsigned int k   ) const { return block_i_[k];   }

  using MatrixScatterAdd< DLA::MatrixD<TM> >::m;
  using MatrixScatterAdd< DLA::MatrixD<TM> >::n;
  int rowNonZero(const int i) const { return row_ptr_[i+1]-row_ptr_[i]; } // number of nonzero components in row i
  int getNumNonZero()         const { return nnz_; }
  int* get_row_ptr()          const { SANS_ASSERT( row_ptr_ != NULL ); return row_ptr_; }
  int* get_col_ind()          const { SANS_ASSERT( col_ind_ != NULL ); return col_ind_; }
  TM* get_values()            const { SANS_ASSERT( val_ != NULL );     return val_; }

protected:
  using MatrixScatterAdd< DLA::MatrixD<TM> >::m_;
  using MatrixScatterAdd< DLA::MatrixD<TM> >::n_;
  int *block_m_; //Block rows sizes
  int *block_n_; //Block column sizes
  int *block_i_; //Block starting index for a non-zero block
  TM *val_;      //Non zero values
  int *row_ind_; //Sparse non-zero to dense matrix index map
  int *col_ind_; //Sparse non-zero to dense matrix index map (typical CRS notation)
  int *row_ptr_; //Index into val_ for start of non-zero row (typical CRS notation)
  int nnz_;      //Total number of non-zero elements

  int *diag_ind_; //The sparse row index for the diagonal elements
};

template<class TM>
void
SparseMatrix_CRS< DLA::MatrixD<TM> >::resize(const SparseNonZeroPattern< DLA::MatrixD<TM> >& Pattern)
{
  //Release the memory
  deallocate();

  //Initialize indexes for CRS storage
  m_ = Pattern.m();
  n_ = Pattern.n();
  nnz_ = Pattern.nnz();

  //No point in doing anything if there are no non-zero elements
  if (nnz_ == 0) return;

  //nnz_ could be negative if we exceed the maximum int value
  SANS_ASSERT( nnz_ > 0 );

  col_ind_ = new int[nnz_];
  row_ptr_ = new int[m_+1];
  diag_ind_ = new int[m_];

  //Initialize the block indexes
  block_m_ = new int[m_];
  block_n_ = new int[n_];
  block_i_ = new int[nnz_];
  row_ind_ = new int[nnz_];


  row_ptr_[0] = 0;
  for ( int i = 0; i < m_; i++ )
  {
    row_ptr_[i+1] = row_ptr_[i] + Pattern.rowSize(i);
    int j = 0;
    for (auto col = Pattern.rowBegin(i); col != Pattern.rowEnd(i); col++, j++)
    {
      col_ind_[row_ptr_[i] + j] = *col;
      if ( i == (int)*col ) diag_ind_[i] = row_ptr_[i] + j;
    }

    //Initialize the block size array and the row index array
    block_m_[i] = Pattern.block_m(i);
    for (int k = row_ptr_[i]; k < row_ptr_[i+1]; k++)
      row_ind_[k] = i;
  }

  for ( int i = 0; i < n_; i++ )
    block_n_[i] = Pattern.block_n(i);

  //Create the array for the memory index for each block
  block_i_[0] = 0;
  for (int k = 1; k < nnz_; k++)
    block_i_[k] = block_i_[k-1] + block_m_[row_ind_[k-1]]*block_n_[col_ind_[k-1]];

  //The size of the last block is lacking in the loop above
  int memorysize = block_i_[nnz_-1] + block_m_[row_ind_[nnz_-1]]*block_n_[col_ind_[nnz_-1]];

  //memorysize could be negative if we exceed the maximum int value
  SANS_ASSERT( memorysize > 0 );

  //Allocate the actual memory storage and set it to zero
  val_ = new TM[ memorysize ];

  for ( int i = 0; i < memorysize; i++)
    val_[i] = 0;

  return;
}

template<class TM>
void
SparseMatrix_CRS< DLA::MatrixD<TM> >::copy_from(const SparseMatrix_CRS& A)
{
  //Release the memory
  deallocate();

  //Initialize indexes for CRS storage
  m_ = A.m_;
  n_ = A.n_;
  nnz_ = A.nnz_;

  //No point in doing anything if there are no non-zero elements
  if (nnz_ == 0) return;

  // Allocate local arrays and copy over everything from A
  col_ind_ = new int[nnz_];
  row_ptr_ = new int[m_+1];
  diag_ind_ = new int[m_];

  //Initialize the block indexes
  block_m_ = new int[m_];
  block_n_ = new int[n_];
  block_i_ = new int[nnz_];
  row_ind_ = new int[nnz_];

  for (int i = 0; i < nnz_; i++)
  {
    col_ind_[i] = A.col_ind_[i];
    block_i_[i] = A.block_i_[i];
    row_ind_[i] = A.row_ind_[i];
  }

  for (int i = 0; i < m_; i++)
  {
    diag_ind_[i] = A.diag_ind_[i];
    block_m_[i] = A.block_m_[i];
    row_ptr_[i] = A.row_ptr_[i];
  }
  row_ptr_[m_] = A.row_ptr_[m_];

  for (int i = 0; i < n_; i++)
    block_n_[i] = A.block_n_[i];

  //The size of the last block is lacking in the loop above
  int memorysize = block_i_[nnz_-1] + block_m_[row_ind_[nnz_-1]]*block_n_[col_ind_[nnz_-1]];

  //Allocate the actual memory storage and copy if from A
  val_ = new TM[ memorysize ];

  for ( int i = 0; i < memorysize; i++)
    val_[i] = A.val_[i];
}


template <class TM>
void
SparseMatrix_CRS< DLA::MatrixD<TM> >::scatterAdd( const DLA::MatrixDView< TM >& M, const int row, const int col )
{
  //Add the values of M to the sparse matrix
  int k = 0;
  const int col_start = row_ptr_[row];
  while ( col != col_ind_[col_start+k] ) k++;
  sparseRow(row, k) += M;
}

template <class TM>
typename SparseMatrix_CRS< DLA::MatrixD<TM> >::MatrixView_type
SparseMatrix_CRS< DLA::MatrixD<TM> >::slowAccess( const int i, const int j )
{
  SANS_ASSERT( i >= 0 && i < m_ );
  SANS_ASSERT( j >= 0 && j < n_ );
  int k = 0;
  const int row = row_ptr_[i];
  const int rownnz = rowNonZero(i);
  while ( k < rownnz && j != col_ind_[row+k] ) k++;
  SANS_ASSERT( k < rownnz );
  return operator[](row + k);
}

template <class TM>
const typename SparseMatrix_CRS< DLA::MatrixD<TM> >::MatrixView_type
SparseMatrix_CRS< DLA::MatrixD<TM> >::slowAccess( const int i, const int j ) const
{
  return const_cast<SparseMatrix_CRS< DLA::MatrixD<TM> >*>(this)->slowAccess(i,j);
}

template <class TM>
const typename SparseMatrix_CRS< DLA::MatrixD<TM> >::MatrixView_type
SparseMatrix_CRS< DLA::MatrixD<TM> >::operator()( const int i, const int j ) const
{
  return const_cast<SparseMatrix_CRS< DLA::MatrixD<TM> >*>(this)->slowAccess(i,j);
}

} //namespace SLA
} //namespace SANS

#define TM_SPEC DLA::MatrixD<TM>
#include "SparseMatrix_CRS_impl.h"
#undef TM_SPEC

#endif //SPARSEMATRIX_CRS_H
