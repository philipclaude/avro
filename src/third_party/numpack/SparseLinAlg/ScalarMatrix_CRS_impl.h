// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

// No include blocker here. This file should only be included once in a cpp file for explicit instantiation

#if !defined(SCALARMATRIX_CRS_INSTANTIATE) && !defined(SANS_HEADERCOMPILE_CHECK)
#error "This file should only be included in a cpp file for explicit instantiations"
#endif

#include <iomanip>
#include <fstream>

#include "ScalarMatrix_CRS.h"
#include "numpack/DenseLinAlg/StaticSize/MatrixS.h"
#include "numpack/DenseLinAlg/tools/index.h"
#include "numpack/DenseLinAlg/tools/MatrixSize.h"

#include "numpack/BlockLinAlg/MatrixBlock_2x2.h"
#include "numpack/BlockLinAlg/MatrixBlock_3x3.h"
#include "numpack/BlockLinAlg/MatrixBlock_4x4.h"



namespace numpack 
{

namespace SLA
{

//-----------------------------------------------------------------------------
inline int blockColumns( const Real& A ) { return 1; }

template<int M, class T>
inline int blockColumns( const DLA::MatrixS< M, M, T >& A ) { return M; }

template<class T>
inline int blockColumns( const DLA::MatrixDView< T >& A )
{
  int col = 0;
  for (int n = 0; n < A.n(); n++)
   col += blockColumns(A(0,n));

  return col;
}

inline int blockColumns( const SparseMatrix_CRS<Real>& A, const unsigned int n ) { return 1; }

inline int blockColumns( const SparseMatrix_CRS<DLA::MatrixD<Real> >& A, const unsigned int n ) { return A.block_n(n); }

template<int M>
inline int blockColumns( const SparseMatrix_CRS<DLA::MatrixD<DLA::MatrixS< M, M, Real > > >& A, const unsigned int n )
{
  return A.block_n(n)*M;
}



//-----------------------------------------------------------------------------
inline int scalarNonZeros( const Real& A ) { return 1; }

template<int M, int N, class T>
int scalarNonZeros( const DLA::MatrixS< M, N, T >& A ) { return M*N; }

template<class T>
int scalarNonZeros( const DLA::MatrixDView< T >& A )
{
  return scalarNonZeros(A(0,0))*A.size();
}

template<class TM>
int scalarNonZeros( const DLA::MatrixDView<SparseMatrix_CRS<TM> >& A )
{
  int nnz_ = 0;
  for ( int i = 0; i < A.m(); i++ )
    for ( int j = 0; j < A.n(); j++ )
      nnz_ += scalarNonZeros(A(i,j));

  return nnz_;
}

template<class TM>
int scalarNonZeros( const SparseMatrix_CRS<TM>& A )
{
  return A.getNumNonZero()*scalarNonZeros(A[0]);
}

template<class TM>
int scalarNonZeros( const SparseMatrix_CRS<DLA::MatrixD<TM> >& A )
{
  const int blocknnz = A.getNumNonZero();
  int nnz_ = 0;
  for ( int k = 0; k < blocknnz; k++ )
    nnz_ += scalarNonZeros(A[k]);

  return nnz_;
}

//-----------------------------------------------------------------------------
template<class INT, class T>
template<class MatrixQ>
void ScalarMatrix_CRS<INT, T>::getScalarRows(const DLA::MatrixDView<SparseMatrix_CRS<MatrixQ> >& A,
                                             DLA::VectorD< int >& row_m)
{
  SANS_ASSERT( A.m() == row_m.m() );

  const int M = DLA::MatrixSize<MatrixQ>::M;

  row_m = 0;

  for ( int i = 0; i < A.m(); i++ )
  {
    row_m[i] = std::max( row_m[i], A(i,0).m() )*M;

    m_ += row_m[i]; //add number of rows in scalar matrix
  }
}

//-----------------------------------------------------------------------------
template<class INT, class T>
template<class MatrixQ>
void ScalarMatrix_CRS<INT, T>::getScalarCols(const DLA::MatrixDView<SparseMatrix_CRS<MatrixQ> >& A,
                                             DLA::VectorD< int >& col_n)
{
  SANS_ASSERT( A.n() == col_n.m() );

  const int An = A.n();
  const int N = DLA::MatrixSize<MatrixQ>::N;

  col_n = 0;

  for ( int j = 0; j < An; j++ )
  {
    col_n[j] = std::max( col_n[j], A(0,j).n() )*N;

    n_ += col_n[j]; //add number of columns in scalar matrix
  }
}

//-----------------------------------------------------------------------------
// add to *this scalar matrix the row(s) di of the row(s) si of the row i of A
template<class INT, class T>
template<class MatrixQ>
void ScalarMatrix_CRS<INT, T>::addrow( const DLA::MatrixDView< SparseMatrix_CRS<MatrixQ> >& A,
                                       int i, int si, int di, int row, int col, int& nz)
{
  // i: row (of CRS sparse matrix) index in the dynamic dense matrix A
  // si: row (of dense static matrix blocks) index in a CRS sparse matrix block of A
  // di: row (of scalar values) index in a dense static matrix block of a CRS sparse matrix block of A
  // row/col: starting row/column (of scalar values) indices in *this scalar matrix to add to
  // nz: index of nonzero scalar entry  in *this scalar matrix

  const int N = DLA::MatrixSize<MatrixQ>::N;

  //Loop over columns in a dynamic matrix
  for ( int j = 0; j < A.n(); j++ )
  {
    //Update the column offset
    if ( j > 0 ) col += A(i,j-1).n()*N;

    addrow( A(i,j), si, di, row, col, nz); // NOTE: nz is passed in by reference and thus updated here
  }
}

//-----------------------------------------------------------------------------
// add to *this scalar matrix the row(s) di of the row(s) si of A
template<class INT, class T>
template<class MatrixQ>
void ScalarMatrix_CRS<INT, T>::addrow( const SparseMatrix_CRS<MatrixQ>& A,
                                       int si, int di, int row, int col, int& nz)
{
  // si: row (of dense static matrix blocks) index in the CRS sparse matrix A
  // di: row (of scalar values) index in a dense static matrix block of A
  // row/col: starting row/column (of scalar values) indices in *this scalar matrix to add to
  // nz: index of nonzero scalar entry  in *this scalar matrix

  //Do nothing if A is an empty matrix
  if ( A.getNumNonZero() == 0 ) return;

  const int N = DLA::MatrixSize<MatrixQ>::N;

  const int *col_ind = A.get_col_ind();
  const int *row_ptr = A.get_row_ptr();

  //Loop over the number of non-zero columns of the si row of the sparse block matrix A
  for ( int sk = row_ptr[si]; sk < row_ptr[si+1]; ++sk )
  {
    addrow(A[sk], di, row, col + col_ind[sk]*N, nz); // NOTE: nz is passed in by reference and thus updated here
  }
}

//-----------------------------------------------------------------------------
// add to *this scalar matrix the row di of A
template<class INT, class T>
template<class MatrixQ>
void ScalarMatrix_CRS<INT, T>::addrow(const MatrixQ& A,
                                      int di, int row, int col, int& nz)
{
  // di: row (of scalar values) index in a dense static matrix block of CRS sparse matrix A, i.e. DLA::MatrixS<M,N,Real>
  // row/col: starting row/column (of scalar values) indices in *this scalar matrix to add to
  // nz: index of nonzero scalar entry  in *this scalar matrix

  const int N = DLA::MatrixSize<MatrixQ>::N;

  //Loop over columns in the static matrix A
  for ( int dj = 0; dj < N; ++dj )
  {
    Ri_[nz] = col++; // Ri_[nz] is assigned the value of col, and then col is incremented
    Rx_[nz] = DLA::index(A,di,dj);
    nz++;
  }
  Rp_[row] += N; // add to column count of the current row of scalar matrix
}


//-----------------------------------------------------------------------------
template<class INT, class T>
template<class MatrixQ>
ScalarMatrix_CRS<INT, T>::ScalarMatrix_CRS( const SparseMatrix_CRS<MatrixQ>& A )
  : m_(0), n_(0), nnz_(0), Rp_(0), Ri_(0), Rx_(0)
{

  //Number of non-zeros in the scalar matrix
  nnz_ = scalarNonZeros(A);

  //Size of a matrix block MatrixQ (which can be a matrix or simply a scalar (i.e. 1x1 matrix))
  const int M = DLA::MatrixSize<MatrixQ>::M;
  const int N = DLA::MatrixSize<MatrixQ>::N;

  //Number of rows/columns in scalar matrix
  m_ = A.m()*M;
  n_ = A.n()*N;

  //Just in case the matrix is actually empty
  if ( nnz_ == 0 ) return;

  //Allocate temporary scalar matrix that will be transposed
  Rp_ = new INT[m_+1];
  Ri_ = new INT[nnz_];
  Rx_ = new double[nnz_];

  int nz = 0;
  int row = 0;
  int col = 0;
  Ri_[0] = 0;
  Rp_[0] = 0;

  //Loop over the rows of blocks
  for ( int si = 0; si < A.m(); si++ )
  {
    //Loop over the rows in a matrix block MatrixQ
    for ( int di = 0; di < M; di++ )
    {
      row++;
      Rp_[row] = Rp_[row-1];

      addrow( A, si, di, row, col, nz );
    }
  }
}


//-----------------------------------------------------------------------------
template<class INT, class T>
template<class MatrixQ>
ScalarMatrix_CRS<INT, T>::ScalarMatrix_CRS( const SparseMatrix_CRS< DLA::MatrixD<MatrixQ> >& A )
  : m_(0), n_(0), nnz_(0), Rp_(0), Ri_(0), Rx_(0)
{
  typedef typename SparseMatrix_CRS< DLA::MatrixD<MatrixQ> >::MatrixView_type MatrixView_type;

  const int *col_ind = A.get_col_ind();
  const int *row_ptr = A.get_row_ptr();

  //Yeah let's just make sure there is something in the matrix as (0,0) is used below
  SANS_ASSERT( A.getNumNonZero() > 0 );
  SANS_ASSERT( row_ptr[1] > 0 );

  //Number of non-zeros in the scalar matrix
  nnz_ = scalarNonZeros(A);

  //Compute the starting index for each column
  int *colBlock_i = new int[A.n()];
  colBlock_i[0] = 0;
  for ( int i = 1; i < A.n(); i++ )
    colBlock_i[i] = colBlock_i[i-1] + blockColumns(A,i-1);

  //Number of columns in scalar matrix
  n_ = colBlock_i[A.n()-1] + blockColumns(A,A.n()-1);

  //Size of a matrix block MatrixQ (which can be a matrix or simply a scalar (i.e. 1x1 matrix))
  const int M = DLA::MatrixSize<MatrixQ>::M;
  const int N = DLA::MatrixSize<MatrixQ>::N;

  //Compute the number of rows
  m_ = 0;
  for ( int i = 0; i < A.m(); i++ )
    m_ += A.block_m(i)*M;

  //Allocate temporary scalar matrix that will be transposed
  Rp_ = new INT[m_+1];
  Ri_ = new INT[nnz_];
  Rx_ = new double[nnz_];


  int nz = 0;
  int row = 0;
  int col = 0;
  Ri_[0] = 0;
  Rp_[0] = 0;

  //Loop over the rows of blocks
  for ( int bi = 0; bi < A.m(); bi++ )
  {
    //Loop over the rows in a dynamic block
    for ( int i = 0; i < A.block_m(bi); i++ )
    {
      //Loop over the rows in a static block
      for ( int si = 0; si < M; si++ )
      {
        row++;
        Rp_[row] = Rp_[row-1];

        //Loop over the number of non-zero columns of blocks
        for ( int bk = row_ptr[bi]; bk < row_ptr[bi+1]; bk++ )
        {
          MatrixView_type Block = A[bk];
          col = colBlock_i[col_ind[bk]];

          //Loop over columns in a dynamic block
          for ( int j = 0; j < Block.n(); j++ )
          {
            //Loop over columns in a static block
            for ( int sj = 0; sj < N; sj++ )
            {
              Ri_[nz] = col++;
              Rx_[nz] = DLA::index(Block(i,j),si,sj);
              nz++;
              Rp_[row]++;
            }
          }
        }
      }
    }
  }

  delete [] colBlock_i;
}


//-----------------------------------------------------------------------------
template<class INT, class T>
template<class MatrixQ>
ScalarMatrix_CRS<INT, T>::ScalarMatrix_CRS( const DLA::MatrixDView< SparseMatrix_CRS<MatrixQ> >& A)
  : m_(0), n_(0), nnz_(0), Rp_(0), Ri_(0), Rx_(0)
{
  //Yeah lets just make sure there is something in the matrix as (0,0) is used below
  SANS_ASSERT( A.m() > 0 );

  const int Am = A.m();
  const int An = A.n();

  //Number of non-zeros in the scalar matrix
  nnz_ = scalarNonZeros(A);

  DLA::VectorD< int > row_m(Am);
  DLA::VectorD< int > col_n(An);

  row_m = 0;
  col_n = 0;

  //Size of a matrix block MatrixQ (which can be a matrix or simply a scalar (i.e. 1x1 matrix))
  const int M = DLA::MatrixSize<MatrixQ>::M;
  const int N = DLA::MatrixSize<MatrixQ>::N;

  m_ = 0;
  for ( int i = 0; i < Am; i++ )
  {
    for ( int j = 0; j < An; j++ )
    {
      row_m[i] = std::max( row_m[i], A(i,j).m() );
      col_n[j] = std::max( col_n[j], A(i,j).n() );
    }

    //Number of rows in scalar matrix
    m_ += row_m[i]*M;
  }

  //Compute the number of columns
  n_ = 0;
  for ( int j = 0; j < An; j++ )
    n_ += col_n[j]*N;

  //Allocate temporary scalar matrix that will be transposed
  Rp_ = new INT[m_+1];
  Ri_ = new INT[nnz_];
  Rx_ = new double[nnz_];


  int nz = 0;
  int row = 0;
  int col = 0;
  Ri_[0] = 0;
  Rp_[0] = 0;

  //Loop over the rows in the dynamic outer matrix
  for ( int i = 0; i < Am; i++ )
  {
    //Loop over the rows of sparse matrices
    for ( int si = 0; si < row_m[i]; si++ )
    {
      //Loop over the rows of dense matrices
      for ( int di = 0; di < M; di++ )
      {
        row++;
        Rp_[row] = Rp_[row-1];

        addrow( A, i, si, di, row, col, nz );
      }
    }
  }

}


//-----------------------------------------------------------------------------
template<class INT, class T>
template<class MatrixQ>
ScalarMatrix_CRS<INT, T>::ScalarMatrix_CRS( const DLA::MatrixDView< SparseMatrix_CRS< DLA::MatrixD<MatrixQ> > >& A)
  : m_(0), n_(0), nnz_(0), Rp_(0), Ri_(0), Rx_(0)
{
  typedef typename SparseMatrix_CRS< DLA::MatrixD<MatrixQ> >::MatrixView_type MatrixView_type;

  //Yeah lets just make sure there is something in the matrix as (0,0) is used below
  SANS_ASSERT( A.m() > 0 );
  SANS_ASSERT( A.m() == A.n() );

  //Number of non-zeros in the scalar matrix
  nnz_ = scalarNonZeros(A);

  DLA::VectorD< int* > colBlock_i(A.n());
  DLA::VectorD< int* > block_m( A.m() );
  DLA::VectorD< int > row_m(A.m());
  DLA::VectorD< int > col_n(A.n());

  row_m = 0;
  col_n = 0;

  //Size of a matrix block MatrixQ (which can be a matrix or simply a scalar (i.e. 1x1 matrix))
  const int M = DLA::MatrixSize<MatrixQ>::M;
  const int N = DLA::MatrixSize<MatrixQ>::N;

  m_ = 0;
  for ( int i = 0; i < A.m(); i++ )
  {
    for ( int j = 0; j < A.n(); j++ )
    {
      row_m[i] = std::max( row_m[i], A(i,j).m() );
      col_n[j] = std::max( col_n[j], A(i,j).n() );
    }

    block_m[i] = new int[ row_m[i] ];
  }

  for ( int i = 0; i < A.m(); i++ )
  {
    for ( int k = 0; k < row_m[i]; k++ )
    {
      int bm = 0;
      for ( int j = 0; j < A.n(); j++ )
        if ( A(i,j).getNumNonZero() )
          bm = std::max( bm, A(i,j).block_m(k) );
      //Total number of rows in the scalar matrix
      m_ += bm*M;
      //Rows in each of the blocks
      block_m[i][k] = bm;
    }
  }



  //Compute the starting index for each column
  for ( int j = 0; j < A.n(); j++ )
  {
    colBlock_i[j] = new int[ col_n[j] ];
    colBlock_i[j][0] = 0;

    int k = j-1;
    if ( k >= 0)
    {
      colBlock_i[j][0] = colBlock_i[k][col_n[k]-1];
      int block_n = 0;
      for ( int i = 0; i < A.m(); i++ )
        if ( col_n[k]-1 < A(i,k).n() && A(i,k).getNumNonZero() > 0 )
          block_n = std::max( block_n, A(i,k).block_n(col_n[k]-1) );

      colBlock_i[j][0] += block_n*N;
    }

    for ( int n = 1; n < col_n[j]; n++ )
    {
      int block_n = 0;
      for ( int i = 0; i < A.m(); i++ )
        if ( n-1 < A(i,j).n() && A(i,j).getNumNonZero() > 0 )
          block_n = std::max( block_n, A(i,j).block_n(n-1) );

      colBlock_i[j][n] = colBlock_i[j][n-1] + block_n*N;

      //Compute the total number of columns excluding the last column of blocks
      n_ += block_n*N;
    }

    k = col_n[j]-1;
    int block_n = 0;
    for ( int i = 0; i < A.m(); i++ )
      if ( k < (int)A(i,j).n() && A(i,j).getNumNonZero() > 0 )
        block_n = std::max( block_n, A(i,j).block_n(k) );

    //Get the size of the last number of columns in the last block column
    n_ += block_n*N;
  }


  //Allocate temporary scalar matrix
  Rp_ = new INT[m_+1];
  Ri_ = new INT[nnz_];
  Rx_ = new double[nnz_];

  int nz = 0;
  int row = 0;
  int col = 0;
  Ri_[0] = 0;
  Rp_[0] = 0;

  //Loop over the rows in the dynamic outer matrix
  for ( int i = 0; i < A.m(); i++ )
  {
    //Loop over the rows of sparse matrices
    for ( int si = 0; si < row_m[i]; si++ )
    {
      //Loop over the rows of dense dynamic matrices
      for ( int di = 0; di < block_m[i][si]; di++ )
      {
        //Loop over the rows of dense static matrices
        for ( int bi = 0; bi < M; bi++ )
        {
          row++;
          Rp_[row] = Rp_[row-1];

          //Loop over columns in a dynamic
          for ( int j = 0; j < A.n(); j++ )
          {
            //Do nothing for an empty matrix
            if ( A(i,j).getNumNonZero() == 0 ) continue;

            const int *col_ind = A(i,j).get_col_ind();
            const int *row_ptr = A(i,j).get_row_ptr();

            //Loop over the number of non-zero columns of blocks
            for ( int sk = row_ptr[si]; sk < row_ptr[si+1]; sk++ )
            {
              col = colBlock_i[j][col_ind[sk]];
              MatrixView_type Block = A(i,j)[sk];

              //Loop over columns in a dynamic block
              for ( int dj = 0; dj < Block.n(); dj++ )
              {
                const DLA::MatrixS<M,N,Real>& SBlock = Block(di,dj);
                for ( int bj = 0; bj < N; bj++ )
                {
                  Ri_[nz] = col++;
                  Rx_[nz] = SBlock(bi,bj);
                  nz++;
                  Rp_[row]++;
                }
              }
            }
          }
        }
      }
    }
  }

  for ( int i = 0; i < A.m(); i++ )
    delete [] block_m[i];

  for ( int j = 0; j < A.n(); j++ )
    delete [] colBlock_i[j];
}


//-----------------------------------------------------------------------------
template<class INT, class T>
template<class MatrixQ00, class MatrixQ01,
         class MatrixQ10, class MatrixQ11>
ScalarMatrix_CRS<INT, T>::ScalarMatrix_CRS(
    const BLA::MatrixBlock_2x2< DLA::MatrixD<SparseMatrix_CRS<MatrixQ00> >,
                                DLA::MatrixD<SparseMatrix_CRS<MatrixQ01> >,
                                DLA::MatrixD<SparseMatrix_CRS<MatrixQ10> >,
                                DLA::MatrixD<SparseMatrix_CRS<MatrixQ11> > >& A )
  : m_(0), n_(0), nnz_(0), Rp_(0), Ri_(0), Rx_(0)
{
  //Number of non-zeros in the scalar matrix
  nnz_ = scalarNonZeros(A.m00) + scalarNonZeros(A.m01) +
         scalarNonZeros(A.m10) + scalarNonZeros(A.m11);

  //Just in case the matrix is actually empty
  if ( nnz_ == 0 ) return;

  DLA::VectorD< int > row_m0(A.m00.m());
  DLA::VectorD< int > row_m1(A.m10.m());

  DLA::VectorD< int > col_n0(A.m00.n());
  DLA::VectorD< int > col_n1(A.m01.n());

  //Number of rows/columns in scalar matrix
  getScalarRows(A.m00, row_m0);
  getScalarRows(A.m10, row_m1);

  getScalarCols(A.m00, col_n0);
  getScalarCols(A.m01, col_n1);

  int colsize_n0 = 0;
  for (int i = 0; i < col_n0.size(); ++i)
    colsize_n0 += col_n0[i];

  //Allocate temporary scalar matrix
  Rp_ = new INT[m_+1];
  Ri_ = new INT[nnz_];
  Rx_ = new double[nnz_];

  int nz = 0;
  int row = 0;
  int col = 0;
  Ri_[0] = 0;
  Rp_[0] = 0;

  SANS_ASSERT(A.m00.m() == A.m01.m());
  SANS_ASSERT(A.m10.m() == A.m11.m());

  SANS_ASSERT(A.m00.n() == A.m10.n());
  SANS_ASSERT(A.m01.n() == A.m11.n());

  // Deduce the sizes of MatrixQ types
  const int M0 = DLA::MatrixSize<MatrixQ00>::M;
  const int M1 = DLA::MatrixSize<MatrixQ10>::M;

  SANS_ASSERT( M0==DLA::MatrixSize<MatrixQ01>::M );

  SANS_ASSERT( M1==DLA::MatrixSize<MatrixQ11>::M );

  //Loop over the rows in the dynamic outer matrix
  for ( int i = 0; i < A.m00.m(); i++ )
  {
    //Loop over the rows of sparse matrices
    for ( int si = 0; si < A.m00(i,0).m(); si++ )
    {
      //Loop over the rows of dense matrices
      for ( int di = 0; di < M0; di++ )
      {
        row++;
        Rp_[row] = Rp_[row-1];

        col = 0;
        addrow( A.m00, i, si, di, row, col, nz );
        col = colsize_n0;
        addrow( A.m01, i, si, di, row, col, nz );
      }
    }
  }

  //Loop over the rows in the dynamic outer matrix
  for ( int i = 0; i < A.m11.m(); i++ )
  {
    //Loop over the rows of sparse matrices
    for ( int si = 0; si < A.m11(i,0).m(); si++ )
    {
      //Loop over the rows of dense matrices
      for ( int di = 0; di < M1; di++ )
      {
        row++;
        Rp_[row] = Rp_[row-1];

        col = 0;
        addrow( A.m10, i, si, di, row, col, nz );
        col = colsize_n0;
        addrow( A.m11, i, si, di, row, col, nz );
      }
    }
  }
}

//-----------------------------------------------------------------------------
template<class INT, class T>
template<class MatrixQ00, class MatrixQ01, class MatrixQ02,
         class MatrixQ10, class MatrixQ11, class MatrixQ12,
         class MatrixQ20, class MatrixQ21, class MatrixQ22>
ScalarMatrix_CRS<INT, T>::ScalarMatrix_CRS(
        const BLA::MatrixBlock_3x3< DLA::MatrixD<SparseMatrix_CRS<MatrixQ00> >,
                                    DLA::MatrixD<SparseMatrix_CRS<MatrixQ01> >,
                                    DLA::MatrixD<SparseMatrix_CRS<MatrixQ02> >,

                                    DLA::MatrixD<SparseMatrix_CRS<MatrixQ10> >,
                                    DLA::MatrixD<SparseMatrix_CRS<MatrixQ11> >,
                                    DLA::MatrixD<SparseMatrix_CRS<MatrixQ12> >,

                                    DLA::MatrixD<SparseMatrix_CRS<MatrixQ20> >,
                                    DLA::MatrixD<SparseMatrix_CRS<MatrixQ21> >,
                                    DLA::MatrixD<SparseMatrix_CRS<MatrixQ22> > >& A )
  : m_(0), n_(0), nnz_(0), Rp_(0), Ri_(0), Rx_(0)
{

  //Number of non-zeros in the scalar matrix
  nnz_ = scalarNonZeros(A.m00) + scalarNonZeros(A.m01) + scalarNonZeros(A.m02) +
         scalarNonZeros(A.m10) + scalarNonZeros(A.m11) + scalarNonZeros(A.m12) +
         scalarNonZeros(A.m20) + scalarNonZeros(A.m21) + scalarNonZeros(A.m22);

  //Just in case the matrix is actually empty
  if ( nnz_ == 0 ) return;

  DLA::VectorD< int > row_m0(A.m00.m());
  DLA::VectorD< int > row_m1(A.m10.m());
  DLA::VectorD< int > row_m2(A.m20.m());

  DLA::VectorD< int > col_n0(A.m00.n());
  DLA::VectorD< int > col_n1(A.m01.n());
  DLA::VectorD< int > col_n2(A.m02.n());

  //Number of rows/columns in scalar matrix
  getScalarRows(A.m00, row_m0);
  getScalarRows(A.m10, row_m1);
  getScalarRows(A.m20, row_m2);

  getScalarCols(A.m00, col_n0);
  getScalarCols(A.m01, col_n1);
  getScalarCols(A.m02, col_n2);

  int colsize_n0 = 0;
  for (int i = 0; i < col_n0.size(); i++)
    colsize_n0 += col_n0[i];

  int colsize_n1 = 0;
  for (int i = 0; i < col_n1.size(); i++)
    colsize_n1 += col_n1[i];

  //Allocate temporary scalar matrix that will be transposed
  Rp_ = new INT[m_+1];
  Ri_ = new INT[nnz_];
  Rx_ = new double[nnz_];

  int nz = 0;
  int row = 0;
  int col = 0;
  Ri_[0] = 0;
  Rp_[0] = 0;

  // Check that all the row sizes match up
  SANS_ASSERT(A.m00.m() == A.m01.m());
  SANS_ASSERT(A.m00.m() == A.m02.m());
  SANS_ASSERT(A.m10.m() == A.m11.m());
  SANS_ASSERT(A.m10.m() == A.m12.m());
  SANS_ASSERT(A.m20.m() == A.m21.m());
  SANS_ASSERT(A.m20.m() == A.m22.m());

  // Check that all column sizes match up
  SANS_ASSERT(A.m00.n() == A.m10.n());
  SANS_ASSERT(A.m00.n() == A.m20.n());
  SANS_ASSERT(A.m01.n() == A.m11.n());
  SANS_ASSERT(A.m01.n() == A.m21.n());
  SANS_ASSERT(A.m02.n() == A.m12.n());
  SANS_ASSERT(A.m02.n() == A.m22.n());

  // Deduce the sizes of MatrixQ types
  const int M0 = DLA::MatrixSize<MatrixQ00>::M;
  const int M1 = DLA::MatrixSize<MatrixQ10>::M;
  const int M2 = DLA::MatrixSize<MatrixQ20>::M;

  SANS_ASSERT( M0==DLA::MatrixSize<MatrixQ01>::M );
  SANS_ASSERT( M0==DLA::MatrixSize<MatrixQ02>::M );

  SANS_ASSERT( M1==DLA::MatrixSize<MatrixQ11>::M );
  SANS_ASSERT( M1==DLA::MatrixSize<MatrixQ12>::M );

  SANS_ASSERT( M2==DLA::MatrixSize<MatrixQ21>::M );
  SANS_ASSERT( M2==DLA::MatrixSize<MatrixQ22>::M );


  //Loop over the rows in the dynamic outer matrix in the first row of blocks
  for ( int i = 0; i < A.m00.m(); i++ )
  {
    //Loop over the rows of sparse matrices
    for ( int si = 0; si < A.m00(i,0).m(); si++ )
    {
      //Loop over the rows of dense matrices
      for ( int di = 0; di < M0; di++ )
      {
        row++;
        Rp_[row] = Rp_[row-1]; // move the row pointer the to the start of a new row

        col = 0;
        addrow( A.m00, i, si, di, row, col, nz ); // NOTE: nz is passed in by reference and thus updated here
        col += colsize_n0;
        addrow( A.m01, i, si, di, row, col, nz );
        col += colsize_n1;
        addrow( A.m02, i, si, di, row, col, nz );
      }
    }
  }

  //Loop over the rows in the dynamic outer matrix in the 2nd row of blocks
  for ( int i = 0; i < A.m10.m(); i++ )
  {
    //Loop over the rows of sparse matrices
    for ( int si = 0; si < A.m10(i,0).m(); si++ )
    {
      //Loop over the rows of dense matrices
      for ( int di = 0; di < M1; di++ )
      {
        row++;
        Rp_[row] = Rp_[row-1]; // move the row pointer the to the start of a new row

        col = 0;
        addrow( A.m10, i, si, di, row, col, nz ); // NOTE: nz is passed in by reference and thus updated here
        col += colsize_n0;
        addrow( A.m11, i, si, di, row, col, nz );
        col += colsize_n1;
        addrow( A.m12, i, si, di, row, col, nz );
      }
    }
  }

  //Loop over the rows in the dynamic outer matrix in the 3rd row of blocks
  for ( int i = 0; i < A.m20.m(); i++ )
  {
    //Loop over the rows of sparse matrices
    for ( int si = 0; si < A.m20(i,0).m(); si++ )
    {
      //Loop over the rows of dense matrices
      for ( int di = 0; di < M2; di++ )
      {
        row++;
        Rp_[row] = Rp_[row-1]; // move the row pointer the to the start of a new row

        col = 0;
        addrow( A.m20, i, si, di, row, col, nz ); // NOTE: nz is passed in by reference and thus updated here
        col += col_n0[0];
        addrow( A.m21, i, si, di, row, col, nz );
        col += col_n1[0];
        addrow( A.m22, i, si, di, row, col, nz );
      }
    }
  }
}

//-----------------------------------------------------------------------------//
// construct 4x4 block matrix from block components
template<class INT, class T>
template<class MatrixQ00, class MatrixQ01, class MatrixQ02, class MatrixQ03,
         class MatrixQ10, class MatrixQ11, class MatrixQ12, class MatrixQ13,
         class MatrixQ20, class MatrixQ21, class MatrixQ22, class MatrixQ23,
         class MatrixQ30, class MatrixQ31, class MatrixQ32, class MatrixQ33>
ScalarMatrix_CRS<INT, T>::ScalarMatrix_CRS(
    const BLA::MatrixBlock_4x4< DLA::MatrixD<SparseMatrix_CRS<MatrixQ00> >,
                                DLA::MatrixD<SparseMatrix_CRS<MatrixQ01> >,
                                DLA::MatrixD<SparseMatrix_CRS<MatrixQ02> >,
                                DLA::MatrixD<SparseMatrix_CRS<MatrixQ03> >,
                                //
                                DLA::MatrixD<SparseMatrix_CRS<MatrixQ10> >,
                                DLA::MatrixD<SparseMatrix_CRS<MatrixQ11> >,
                                DLA::MatrixD<SparseMatrix_CRS<MatrixQ12> >,
                                DLA::MatrixD<SparseMatrix_CRS<MatrixQ13> >,
                                //
                                DLA::MatrixD<SparseMatrix_CRS<MatrixQ20> >,
                                DLA::MatrixD<SparseMatrix_CRS<MatrixQ21> >,
                                DLA::MatrixD<SparseMatrix_CRS<MatrixQ22> >,
                                DLA::MatrixD<SparseMatrix_CRS<MatrixQ23> >,
                                //
                                DLA::MatrixD<SparseMatrix_CRS<MatrixQ30> >,
                                DLA::MatrixD<SparseMatrix_CRS<MatrixQ31> >,
                                DLA::MatrixD<SparseMatrix_CRS<MatrixQ32> >,
                                DLA::MatrixD<SparseMatrix_CRS<MatrixQ33> > >& A )
  : m_(0), n_(0), nnz_(0), Rp_(0), Ri_(0), Rx_(0)
{
  //Number of non-zeros in the scalar matrix
  nnz_ = scalarNonZeros(A.m00) + scalarNonZeros(A.m01) + scalarNonZeros(A.m02) + scalarNonZeros(A.m03)
       + scalarNonZeros(A.m10) + scalarNonZeros(A.m11) + scalarNonZeros(A.m12) + scalarNonZeros(A.m13)
       + scalarNonZeros(A.m20) + scalarNonZeros(A.m21) + scalarNonZeros(A.m22) + scalarNonZeros(A.m23)
       + scalarNonZeros(A.m30) + scalarNonZeros(A.m31) + scalarNonZeros(A.m32) + scalarNonZeros(A.m33);

  if ( nnz_ == 0 ) return; //Just in case the matrix is actually empty

  //rows/columns in scalar matrix
  // Check that all the row sizes match up
  SANS_ASSERT(A.m00.m() == A.m01.m()); SANS_ASSERT(A.m00.m() == A.m02.m()); SANS_ASSERT(A.m00.m() == A.m03.m());
  SANS_ASSERT(A.m10.m() == A.m11.m()); SANS_ASSERT(A.m10.m() == A.m12.m()); SANS_ASSERT(A.m10.m() == A.m13.m());
  SANS_ASSERT(A.m20.m() == A.m21.m()); SANS_ASSERT(A.m20.m() == A.m22.m()); SANS_ASSERT(A.m20.m() == A.m23.m());
  SANS_ASSERT(A.m30.m() == A.m31.m()); SANS_ASSERT(A.m30.m() == A.m32.m()); SANS_ASSERT(A.m30.m() == A.m33.m());
  // Check that all the column sizes match up
  SANS_ASSERT(A.m00.n() == A.m10.n()); SANS_ASSERT(A.m00.n() == A.m20.n()); SANS_ASSERT(A.m00.n() == A.m30.n());
  SANS_ASSERT(A.m01.n() == A.m11.n()); SANS_ASSERT(A.m01.n() == A.m21.n()); SANS_ASSERT(A.m01.n() == A.m31.n());
  SANS_ASSERT(A.m02.n() == A.m12.n()); SANS_ASSERT(A.m02.n() == A.m22.n()); SANS_ASSERT(A.m02.n() == A.m32.n());
  SANS_ASSERT(A.m03.n() == A.m13.n()); SANS_ASSERT(A.m03.n() == A.m23.n()); SANS_ASSERT(A.m03.n() == A.m33.n());

  DLA::VectorD<int> row_m0(A.m00.m()), row_m1(A.m10.m()), row_m2(A.m20.m()), row_m3(A.m30.m());
  DLA::VectorD<int> col_n0(A.m00.n()), col_n1(A.m01.n()), col_n2(A.m02.n()), col_n3(A.m03.n());

  // get row/col counts and update m_/n_ under the hood
  getScalarRows(A.m00, row_m0);
  getScalarRows(A.m10, row_m1);
  getScalarRows(A.m20, row_m2);
  getScalarRows(A.m30, row_m3);

  getScalarCols(A.m00, col_n0);
  getScalarCols(A.m01, col_n1);
  getScalarCols(A.m02, col_n2);
  getScalarCols(A.m03, col_n3);

  // count the scalar columns in MatrixD blocks
  int colsize_n0 = 0;
  for (int i = 0; i < col_n0.size(); ++i)
    colsize_n0 += col_n0[i];

  int colsize_n1 = 0;
  for (int i = 0; i < col_n1.size(); ++i)
    colsize_n1 += col_n1[i];

  int colsize_n2 = 0;
  for (int i = 0; i < col_n2.size(); ++i)
    colsize_n2 += col_n2[i];

  //Allocate memory for scalar matrix
  Rp_ = new INT[m_+1];
  Ri_ = new INT[nnz_];
  Rx_ = new double[nnz_];

  int nz = 0;
  int row = 0;
  int col = 0;
  Rp_[0] = 0;
  Ri_[0] = 0;

  // Deduce the sizes of MatrixQ types
  const int M0 = DLA::MatrixSize<MatrixQ00>::M;
  const int M1 = DLA::MatrixSize<MatrixQ10>::M;
  const int M2 = DLA::MatrixSize<MatrixQ20>::M;
  const int M3 = DLA::MatrixSize<MatrixQ30>::M;

  SANS_ASSERT( M0==DLA::MatrixSize<MatrixQ01>::M );
  SANS_ASSERT( M0==DLA::MatrixSize<MatrixQ02>::M );
  SANS_ASSERT( M0==DLA::MatrixSize<MatrixQ03>::M );

  SANS_ASSERT( M1==DLA::MatrixSize<MatrixQ11>::M );
  SANS_ASSERT( M1==DLA::MatrixSize<MatrixQ12>::M );
  SANS_ASSERT( M1==DLA::MatrixSize<MatrixQ13>::M );

  SANS_ASSERT( M2==DLA::MatrixSize<MatrixQ21>::M );
  SANS_ASSERT( M2==DLA::MatrixSize<MatrixQ22>::M );
  SANS_ASSERT( M2==DLA::MatrixSize<MatrixQ23>::M );

  SANS_ASSERT( M3==DLA::MatrixSize<MatrixQ31>::M );
  SANS_ASSERT( M3==DLA::MatrixSize<MatrixQ32>::M );
  SANS_ASSERT( M3==DLA::MatrixSize<MatrixQ33>::M );

  //Loop over the rows in the dynamic outer matrix for the first row of blocks
  for ( int i = 0; i < A.m00.m(); ++i )
  {
    //Loop over the rows of sparse matrices
    for ( int si = 0; si < A.m00(i,0).m(); ++si )
    {
      //Loop over the rows of static dense matrices
      for ( int di = 0; di < M0; ++di )
      {
        ++row;
        Rp_[row] = Rp_[row-1]; // move the row pointer the to the start of a new row

        col = 0;
        addrow( A.m00, i, si, di, row, col, nz ); // NOTE: nz is passed in by reference and thus updated here
        col += colsize_n0;
        addrow( A.m01, i, si, di, row, col, nz );
        col += colsize_n1;
        addrow( A.m02, i, si, di, row, col, nz );
        col += colsize_n2;
        addrow( A.m03, i, si, di, row, col, nz );
      }
    }
  }

  //Loop over the rows in the dynamic outer matrix for the second row of blocks
  for ( int i = 0; i < A.m10.m(); ++i )
  {
    //Loop over the rows of sparse matrices
    for ( int si = 0; si < A.m10(i,0).m(); ++si )
    {
      //Loop over the rows of static dense matrices
      for ( int di = 0; di < M1; ++di )
      {
        ++row;
        Rp_[row] = Rp_[row-1]; // move the row pointer the to the start of a new row

        col = 0;
        addrow( A.m10, i, si, di, row, col, nz ); // NOTE: nz is passed in by reference and thus updated here
        col += colsize_n0;
        addrow( A.m11, i, si, di, row, col, nz );
        col += colsize_n1;
        addrow( A.m12, i, si, di, row, col, nz );
        col += colsize_n2;
        addrow( A.m13, i, si, di, row, col, nz );
      }
    }
  }

  //Loop over the rows in the dynamic outer matrix for the second row of blocks
  for ( int i = 0; i < A.m20.m(); ++i )
  {
    //Loop over the rows of sparse matrices
    for ( int si = 0; si < A.m20(i,0).m(); ++si )
    {
      //Loop over the rows of static dense matrices
      for ( int di = 0; di < M2; ++di )
      {
        ++row;
        Rp_[row] = Rp_[row-1]; // move the row pointer the to the start of a new row

        col = 0;
        addrow( A.m20, i, si, di, row, col, nz ); // NOTE: nz is passed in by reference and thus updated here
        col += colsize_n0;
        addrow( A.m21, i, si, di, row, col, nz );
        col += colsize_n1;
        addrow( A.m22, i, si, di, row, col, nz );
        col += colsize_n2;
        addrow( A.m23, i, si, di, row, col, nz );
      }
    }
  }

  //Loop over the rows in the dynamic outer matrix for the second row of blocks
  for ( int i = 0; i < A.m30.m(); ++i )
  {
    //Loop over the rows of sparse matrices
    for ( int si = 0; si < A.m30(i,0).m(); ++si )
    {
      //Loop over the rows of static dense matrices
      for ( int di = 0; di < M3; ++di )
      {
        ++row;
        Rp_[row] = Rp_[row-1]; // move the row pointer the to the start of a new row

        col = 0;
        addrow( A.m30, i, si, di, row, col, nz ); // NOTE: nz is passed in by reference and thus updated here
        col += colsize_n0;
        addrow( A.m31, i, si, di, row, col, nz );
        col += colsize_n1;
        addrow( A.m32, i, si, di, row, col, nz );
        col += colsize_n2;
        addrow( A.m33, i, si, di, row, col, nz );
      }
    }
  }
}

//-----------------------------------------------------------------------------
template<class INT, class T>
void
ScalarMatrix_CRS<INT, T>::WriteMatrixMarketFile( const std::string& filename ) const
{
  std::ofstream file(filename);
  this->WriteMatrixMarketFile(file);
}

//-----------------------------------------------------------------------------
template<class INT, class T>
void
ScalarMatrix_CRS<INT, T>::WriteMatrixMarketFile( std::ostream& file ) const
{
  //Write the banner
  file << "%%MatrixMarket matrix coordinate real general" << std::endl;
  file << m_ << " " << n_ << " " << nnz_ << std::endl;

  //Write out the matrix data
  //Add one to the column index as the file format is 1-based
  for ( int row = 0; row < m_; row++ )
    for ( int k = Rp_[row]; k < Rp_[row+1]; k++ )
      file << std::scientific << std::setprecision( 16 ) << row+1 << " " << Ri_[k]+1 << " " << (std::abs(Rx_[k]) < 1e-12 ? 0 : Rx_[k]) << std::endl;
}


} // namespace SLA
} // namespace numpack 
