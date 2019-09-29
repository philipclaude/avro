// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef SCALAR_MATRIX_CRS
#define SCALAR_MATRIX_CRS

#include "SparseMatrix_CRS.h"
#include "numpack/dense/dynamic/MatrixD.h"
#include "numpack/block/block_Type.h"

namespace numpack 
{
namespace SLA
{

//-----------------------------------------------------------------------------
// Scalar matrix data type
// - entries (i.e. smallest component) are scalar
// - structured in compressed row storage (CRS) format. Ref. http://netlib.org/linalg/html_templates/node91.html

template<class INT, class T = double> // class INT differentiates different integer types
class ScalarMatrix_CRS
{
public:
  // constructs a CRS scalar matrix from various custom matrix types
  template<class MatrixQ>
  ScalarMatrix_CRS( const SparseMatrix_CRS<MatrixQ>& A );
  template<class MatrixQ>
  ScalarMatrix_CRS( const SparseMatrix_CRS< DLA::MatrixD<MatrixQ> >& A );

  template<class MatrixQ>
  ScalarMatrix_CRS( const DLA::MatrixDView< SparseMatrix_CRS<MatrixQ> >& A );
  template<class MatrixQ>
  ScalarMatrix_CRS( const DLA::MatrixDView< SparseMatrix_CRS< DLA::MatrixD<MatrixQ> > >& A );

  template<class MatrixQ00, class MatrixQ01,
           class MatrixQ10, class MatrixQ11>
  ScalarMatrix_CRS(
    const BLA::MatrixBlock_2x2< DLA::MatrixD<SparseMatrix_CRS<MatrixQ00> >,
                                DLA::MatrixD<SparseMatrix_CRS<MatrixQ01> >,
                                DLA::MatrixD<SparseMatrix_CRS<MatrixQ10> >,
                                DLA::MatrixD<SparseMatrix_CRS<MatrixQ11> > >& A );

  template<class MatrixQ00, class MatrixQ01, class MatrixQ02,
           class MatrixQ10, class MatrixQ11, class MatrixQ12,
           class MatrixQ20, class MatrixQ21, class MatrixQ22>
  ScalarMatrix_CRS(
      const BLA::MatrixBlock_3x3< DLA::MatrixD<SparseMatrix_CRS<MatrixQ00> >,
                                  DLA::MatrixD<SparseMatrix_CRS<MatrixQ01> >,
                                  DLA::MatrixD<SparseMatrix_CRS<MatrixQ02> >,

                                  DLA::MatrixD<SparseMatrix_CRS<MatrixQ10> >,
                                  DLA::MatrixD<SparseMatrix_CRS<MatrixQ11> >,
                                  DLA::MatrixD<SparseMatrix_CRS<MatrixQ12> >,

                                  DLA::MatrixD<SparseMatrix_CRS<MatrixQ20> >,
                                  DLA::MatrixD<SparseMatrix_CRS<MatrixQ21> >,
                                  DLA::MatrixD<SparseMatrix_CRS<MatrixQ22> > >& A );

  template<class MatrixQ00, class MatrixQ01, class MatrixQ02, class MatrixQ03,
           class MatrixQ10, class MatrixQ11, class MatrixQ12, class MatrixQ13,
           class MatrixQ20, class MatrixQ21, class MatrixQ22, class MatrixQ23,
           class MatrixQ30, class MatrixQ31, class MatrixQ32, class MatrixQ33>
  ScalarMatrix_CRS(
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
                                  DLA::MatrixD<SparseMatrix_CRS<MatrixQ33> > >& A );


  ScalarMatrix_CRS(const ScalarMatrix_CRS&) = delete;
  ScalarMatrix_CRS& operator=(const ScalarMatrix_CRS&) = delete;

  ~ScalarMatrix_CRS()
  {
    delete [] Rp_;
    delete [] Ri_;
    delete [] Rx_;
  }

  void WriteMatrixMarketFile( const std::string& filename ) const;
  void WriteMatrixMarketFile( std::ostream& file ) const;

  int size() const { return m_*n_; }
  int m() const { return m_; }
  int n() const { return n_; }
  int nnz() const { return nnz_; }

  INT* Rp() const { return Rp_; }
  INT* Ri() const { return Ri_; }
  T* Rx() const { return Rx_; }

protected:
  template<class MatrixQ>
  void getScalarRows(const DLA::MatrixDView<SparseMatrix_CRS<MatrixQ> >& A, DLA::VectorD< int >& row_m);
  template<class MatrixQ>
  void getScalarCols(const DLA::MatrixDView<SparseMatrix_CRS<MatrixQ> >& A, DLA::VectorD< int >& col_n);

  template<class MatrixQ>
  void addrow( const DLA::MatrixDView<SparseMatrix_CRS<MatrixQ> >& A, int i, int si, int di, int row, int col, int& nz);
  template<class MatrixQ>
  void addrow( const SparseMatrix_CRS<MatrixQ>& A, int si, int di, int row, int col, int& nz);
  template<class MatrixQ>
  void addrow( const MatrixQ& A, int di, int row, int col, int& nz);

  INT m_, n_; //Matrix dimensions
  INT nnz_; //Total number of non-zero entries
  INT *Rp_; //Indices/locations (in non-zero entries) that start a row. size = m_+1 array
  INT *Ri_; //Column indexes of non-zero entries. size = nnz_ array
  T *Rx_; //Values non-zero entries. size = nnz_ array
};

} //namespace SLA


namespace DLA
{
//Constructs a dense matrix from a sparse matrix
template<class T>
MatrixD<T>::MatrixD( const SLA::ScalarMatrix_CRS<int, Real>& A )
  : MatrixDView<T>(new T[A.size()], A.m(), A.n())
{
  (*this) = 0;
  for ( int row = 0; row < m_; row++ )
    for ( int k = A.Rp()[row]; k < A.Rp()[row+1]; k++ )
      (*this)(row,  A.Ri()[k]) =  A.Rx()[k];
}
} // namespace DLA


} // namespace numpack 

#endif //SCALAR_MATRIX_CRS
