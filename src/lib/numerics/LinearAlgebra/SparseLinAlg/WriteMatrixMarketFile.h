// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef WRITEMATRIXMARKETFILE_H
#define WRITEMATRIXMARKETFILE_H

#include "SparseMatrix_CRS.h"
#include "SparseNonZeroPattern.h"
#include "ScalarMatrix_CRS.h"

#include <fstream>
#include <iomanip>

namespace SANS
{
namespace SLA
{

void WriteMatrixMarketFile( const DLA::MatrixDView<Real>& A, std::ostream& file );

template< int M, int N >
void WriteMatrixMarketFile( const DLA::MatrixS<M,N,Real>& A, std::ostream& file )
{
  //Write the banner
  file << "%%MatrixMarket matrix coordinate real general" << std::endl;
  file << M << " " << N << " " << M*N << std::endl;

  //Write out the matrix data
  //Add one to the row and column index as the file format is 1-based
  for ( int row = 0; row < M; row++ )
    for ( int col = 0; col < N; col++ )
      file << std::scientific << std::setprecision( 16 ) << row+1 << " " << col+1 << " " << A(row,col) << std::endl;
}

template< int M, int N >
void WriteMatrixMarketFile( const DLA::MatrixDView<DLA::MatrixS<M,N,Real>>& A, std::ostream& file )
{
  const int m = A.m()*M;
  const int n = A.n()*N;

  //Write the banner
  file << "%%MatrixMarket matrix coordinate real general" << std::endl;
  file << m << " " << n << " " << m*n << std::endl;

  //Write out the matrix data
  //Add one to the row and column index as the file format is 1-based
  for ( int i = 0; i < A.m(); i++ )
    for ( int m = 0; m < M; m++ )
      for ( int j = 0; j < A.n(); j++ )
        for ( int n = 0; n < N; n++ )
        {
          const int row = i*M+m + 1;
          const int col = j*N+n + 1;
          file << std::scientific << std::setprecision( 16 ) << row << " " << col << " " << A(i,j)(m,n) << std::endl;
        }
}


template< int M, int N >
void WriteMatrixMarketFile( const SparseMatrix_CRS<DLA::MatrixS<M,N,Real>>& A, std::ostream& file )
{
  const int m = A.m()*M;
  const int n = A.n()*N;

  //Write the banner
  file << "%%MatrixMarket matrix coordinate real general" << std::endl;
  file << m << " " << n << " " << m*n << std::endl;

  //Write out the matrix data
  //Add one to the row and column index as the file format is 1-based
  for ( int i = 0; i < A.m(); i++ )
    for ( int m = 0; m < M; m++ )
      for ( int j = 0; j < A.n(); j++ )
        for ( int n = 0; n < N; n++ )
        {
          const int row = i*M+m + 1;
          const int col = j*N+n + 1;
          file << std::scientific << std::setprecision( 16 ) << row << " " << col << " " << A(i,j)(m,n) << std::endl;
        }
}

template< class T >
void WriteMatrixMarketFile( const SparseMatrix_CRS<T>& A, std::ostream& file )
{
  ScalarMatrix_CRS<int> sA(A);
  sA.WriteMatrixMarketFile(file);
}

template< class T >
void WriteMatrixMarketFile( const SparseMatrix_CRS_Transpose<T>& A, std::ostream& file )
{
  SANS_DEVELOPER_EXCEPTION("Can't dump transpose yet");
//  ScalarMatrix_CRS<int> sA(A);
//  sA.WriteMatrixMarketFile(file);
}

template< class T >
void WriteMatrixMarketFile( const DLA::MatrixDView< SparseMatrix_CRS<T> >& A, std::ostream& file )
{
  ScalarMatrix_CRS<int> sA(A);
  sA.WriteMatrixMarketFile(file);
}

template< class T >
void WriteMatrixMarketFile( const SparseNonZeroPattern<T>& nz, std::ostream& file )
{
  SparseMatrix_CRS<T> A(nz);
  A = 1;
  ScalarMatrix_CRS<int> sA(A);
  sA.WriteMatrixMarketFile(file);
}

template< class T >
void WriteMatrixMarketFile( const DLA::DenseNonZeroPattern<T>& nz, std::ostream& file )
{
  SANS_DEVELOPER_EXCEPTION("Can't dump dense nz yet");
}

template< class T >
void WriteMatrixMarketFile( const SparseNonZeroPattern_Transpose<T>& nz, std::ostream& file )
{
  SANS_DEVELOPER_EXCEPTION("Can't dump transpose yet");
//  SparseMatrix_CRS<T> A(nz);
//  A = 1;
//  ScalarMatrix_CRS<int> sA(A);
//  sA.WriteMatrixMarketFile(file);
}

template< class T >
void WriteMatrixMarketFile( const DLA::MatrixDView< SparseNonZeroPattern<T> >& nz, std::ostream& file )
{
  DLA::MatrixD< SparseMatrix_CRS<T> > A(nz);
  A = 1;
  ScalarMatrix_CRS<int> sA(A);
  sA.WriteMatrixMarketFile(file);
}

// specialization for 2x2 block matrix
template< class T00, class T01, class T10, class T11 >
void WriteMatrixMarketFile(
  const BLA::MatrixBlock_2x2<DLA::MatrixD<SLA::SparseNonZeroPattern<T00> >, DLA::MatrixD<SLA::SparseNonZeroPattern<T01> >,
                             DLA::MatrixD<SLA::SparseNonZeroPattern<T10> >, DLA::MatrixD<SLA::SparseNonZeroPattern<T11> > >& nz,
  std::ostream& file )
{
  BLA::MatrixBlock_2x2<DLA::MatrixD<SLA::SparseMatrix_CRS<T00> >, DLA::MatrixD<SLA::SparseMatrix_CRS<T01> >,
                       DLA::MatrixD<SLA::SparseMatrix_CRS<T10> >, DLA::MatrixD<SLA::SparseMatrix_CRS<T11> > > A(nz);
  A = 1;
  ScalarMatrix_CRS<int> sA(A);
  sA.WriteMatrixMarketFile(file);
}

template< class T00, class T01, class T10, class T11 >
void WriteMatrixMarketFile(const BLA::MatrixBlock_2x2<DLA::MatrixD<SLA::SparseMatrix_CRS<T00> >, DLA::MatrixD<SLA::SparseMatrix_CRS<T01> >,
                                                      DLA::MatrixD<SLA::SparseMatrix_CRS<T10> >, DLA::MatrixD<SLA::SparseMatrix_CRS<T11> > >& A,
                           std::ostream& file )
{
  ScalarMatrix_CRS<int> sA(A);
  sA.WriteMatrixMarketFile(file);
}

// specialization for 4x4 block matrix
template< class T00, class T01, class T02, class T03,
          class T10, class T11, class T12, class T13,
          class T20, class T21, class T22, class T23,
          class T30, class T31, class T32, class T33 >
void WriteMatrixMarketFile(
  const BLA::MatrixBlock_4x4<DLA::MatrixD<SLA::SparseNonZeroPattern<T00> >, DLA::MatrixD<SLA::SparseNonZeroPattern<T01> >,
                             DLA::MatrixD<SLA::SparseNonZeroPattern<T02> >, DLA::MatrixD<SLA::SparseNonZeroPattern<T03> >,
                             DLA::MatrixD<SLA::SparseNonZeroPattern<T10> >, DLA::MatrixD<SLA::SparseNonZeroPattern<T11> >,
                             DLA::MatrixD<SLA::SparseNonZeroPattern<T12> >, DLA::MatrixD<SLA::SparseNonZeroPattern<T13> >,
                             DLA::MatrixD<SLA::SparseNonZeroPattern<T20> >, DLA::MatrixD<SLA::SparseNonZeroPattern<T21> >,
                             DLA::MatrixD<SLA::SparseNonZeroPattern<T22> >, DLA::MatrixD<SLA::SparseNonZeroPattern<T23> >,
                             DLA::MatrixD<SLA::SparseNonZeroPattern<T30> >, DLA::MatrixD<SLA::SparseNonZeroPattern<T31> >,
                             DLA::MatrixD<SLA::SparseNonZeroPattern<T32> >, DLA::MatrixD<SLA::SparseNonZeroPattern<T33> > >& nz,
  std::ostream& file )
{
  BLA::MatrixBlock_4x4< DLA::MatrixD<SLA::SparseNonZeroPattern<T00> >, DLA::MatrixD<SLA::SparseNonZeroPattern<T01> >,
                        DLA::MatrixD<SLA::SparseNonZeroPattern<T02> >, DLA::MatrixD<SLA::SparseNonZeroPattern<T03> >,
                        DLA::MatrixD<SLA::SparseNonZeroPattern<T10> >, DLA::MatrixD<SLA::SparseNonZeroPattern<T11> >,
                        DLA::MatrixD<SLA::SparseNonZeroPattern<T12> >, DLA::MatrixD<SLA::SparseNonZeroPattern<T13> >,
                        DLA::MatrixD<SLA::SparseNonZeroPattern<T20> >, DLA::MatrixD<SLA::SparseNonZeroPattern<T21> >,
                        DLA::MatrixD<SLA::SparseNonZeroPattern<T22> >, DLA::MatrixD<SLA::SparseNonZeroPattern<T23> >,
                        DLA::MatrixD<SLA::SparseNonZeroPattern<T30> >, DLA::MatrixD<SLA::SparseNonZeroPattern<T31> >,
                        DLA::MatrixD<SLA::SparseNonZeroPattern<T32> >, DLA::MatrixD<SLA::SparseNonZeroPattern<T33> > > A(nz);
  A = 1;
  ScalarMatrix_CRS<int> sA(A);
  sA.WriteMatrixMarketFile(file);
}

template< class T00, class T01, class T02, class T03,
          class T10, class T11, class T12, class T13,
          class T20, class T21, class T22, class T23,
          class T30, class T31, class T32, class T33 >
void WriteMatrixMarketFile(
    const BLA::MatrixBlock_4x4<DLA::MatrixD<SLA::SparseMatrix_CRS<T00> >, DLA::MatrixD<SLA::SparseMatrix_CRS<T01> >,
                               DLA::MatrixD<SLA::SparseMatrix_CRS<T02> >, DLA::MatrixD<SLA::SparseMatrix_CRS<T03> >,
                               DLA::MatrixD<SLA::SparseMatrix_CRS<T10> >, DLA::MatrixD<SLA::SparseMatrix_CRS<T11> >,
                               DLA::MatrixD<SLA::SparseMatrix_CRS<T12> >, DLA::MatrixD<SLA::SparseMatrix_CRS<T13> >,
                               DLA::MatrixD<SLA::SparseMatrix_CRS<T20> >, DLA::MatrixD<SLA::SparseMatrix_CRS<T21> >,
                               DLA::MatrixD<SLA::SparseMatrix_CRS<T22> >, DLA::MatrixD<SLA::SparseMatrix_CRS<T23> >,
                               DLA::MatrixD<SLA::SparseMatrix_CRS<T30> >, DLA::MatrixD<SLA::SparseMatrix_CRS<T31> >,
                               DLA::MatrixD<SLA::SparseMatrix_CRS<T32> >, DLA::MatrixD<SLA::SparseMatrix_CRS<T33> > >& A,
    std::ostream& file )
{
  ScalarMatrix_CRS<int> sA(A);
  sA.WriteMatrixMarketFile(file);
}

// write matrix given filename
template< class SparseMatrix >
void WriteMatrixMarketFile( const SparseMatrix& A, const std::string& filename )
{
  std::fstream file(filename, std::ios::out);
  WriteMatrixMarketFile(A, file);
}
}
}

#endif //WRITEMATRIXMARKETFILE_H
