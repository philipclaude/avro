// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef SPARSEMATRIX_DIAG_H
#define SPARSEMATRIX_DIAG_H

#include "SparseMatrix_CRS.h"

namespace SANS
{
namespace SLA
{

// SparseMatrix_Diag represents diagonal matrix

template< class TM >
class SparseMatrix_Diag
{
public:
  SparseMatrix_Diag() : m_ (0), val_(NULL) {}

  void resize( const SparseMatrix_CRS<TM>& M )
  {
    deallocate();

    //Assume we are working with a square matrix
    SANS_ASSERT( M.m() == M.n() );

    m_ = M.m();
    val_ = new TM[m_];
  }

  ~SparseMatrix_Diag() { deallocate(); }

  //Operator to directly access the diagonal values
        TM& operator[](unsigned int k)       { return val_[k]; }
  const TM& operator[](unsigned int k) const { return val_[k]; }

  int m() const { return m_; }
  int n() const { return m_; }

protected:

  //Release all the memory and rest the pointers
  void deallocate()
  {
    m_ = 0;
    delete [] val_; val_ = NULL;
  }

protected:
  unsigned int m_; //Number of rows in the diagonal matrix
  TM *val_;        //The diagonal values
};

//=============================================================================
template< class TM >
class SparseMatrix_Diag< DLA::MatrixD<TM> >
{
public:
  typedef TM T;
  typedef DLA::MatrixDView<TM> MatrixView_type;

  SparseMatrix_Diag() : m_ (0), block_m_(NULL), block_n_(NULL), block_i_(NULL), val_(NULL) {}

  void resize( const SparseMatrix_CRS< DLA::MatrixD<TM> >& M )
  {
    deallocate();

    //Assume we are working with a square matrix
    SANS_ASSERT( M.m() == M.n() );

    unsigned int m = m_ = M.m();

    // Nothing to do with an empty matrix (also suppress clang analyzer)
    if ( M.getNumNonZero() == 0 )
    {
      // suppress clang analyzer warning
      SANS_ASSERT( m_ == 0 );
      return;
    }

    block_m_ = new unsigned int[m_];
    block_n_ = new unsigned int[m_];
    block_i_ = new unsigned int[m_];

    // Using a local iterator variable improves performance
    for ( unsigned int i = 0; i < m; i++ )
    {
      block_m_[i] = M.block_m(i);
      block_n_[i] = M.block_n(i);
    }

    //Create the array for the memory index for each block
    block_i_[0] = 0;
    for ( unsigned int k = 1; k < m; k++)
      block_i_[k] = block_i_[k-1] + block_m_[k-1]*block_n_[k-1];

    //Allocate the memory now that we know the total memory size
    val_ = new TM[block_i_[m_-1] + block_m_[m_-1]*block_n_[m_-1]];
  }

  ~SparseMatrix_Diag() { deallocate(); }

  //Operator to directly access the diagonal values
  MatrixView_type operator[](unsigned int k)
  {
    return MatrixView_type( val_ + block_i_[k], block_m_[k], block_n_[k]);
  }
  const MatrixView_type operator[](unsigned int k) const
  {
    return MatrixView_type( val_ + block_i_[k], block_m_[k], block_n_[k]);
  }

  unsigned int block_m( unsigned int row ) const { return block_m_[row]; }
  unsigned int block_n( unsigned int col ) const { return block_n_[col]; }
  unsigned int block_i( unsigned int k   ) const { return block_i_[k];   }

  int m()                     const { return m_; }
  int n()                     const { return m_; }

protected:

  //Release all the memory and rest the pointers
  void deallocate()
  {
    m_ = 0;
    delete [] val_;     val_ = NULL;
    delete [] block_m_; block_m_ = NULL;
    delete [] block_n_; block_n_ = NULL;
    delete [] block_i_; block_i_ = NULL;
  }

protected:
  unsigned int m_;         //Number of rows in the diagonal matrix
  unsigned int *block_m_; //Block rows sizes
  unsigned int *block_n_; //Block column sizes
  unsigned int *block_i_; //Block starting index for a non-zero block
  TM *val_;               //The diagonal values
};

}
}

#endif //SPARSEMATRIX_DIAG_H
