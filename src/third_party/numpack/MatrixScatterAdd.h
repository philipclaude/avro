// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXSCATTERADD_H
#define MATRIXSCATTERADD_H

#include <typeinfo>

#include "tools/SANSException.h"
#include "TransposeTraits.h"
#include "DenseLinAlg/DynamicSize/MatrixD_Type.h"
#include "SparseLinAlg/SparseLinAlg_Type.h"

namespace numpack 
{

namespace LA
{
enum DerivedMatrix
{
  eMatrixD,
  eMatrixDTranspose,
  eDenseNonZeroPattern,
  eSparseMatrix_CRS,
  eSparseMatrix_CRS_Transpose,
  eSparseNonZeroPattern,
  eSparseNonZeroPattern_Transpose
};
}

template<class T>
struct MatrixScatterAdd
{

  void scatterAdd( const DLA::MatrixDView< T >& M, const int Map[], const int nMap )
  {
    switch (eMatrix_)
    {
    case LA::eMatrixD:
      static_cast<DLA::MatrixD<T>*>(this)->scatterAdd(M, Map, nMap);
      break;
    case LA::eMatrixDTranspose:
      static_cast<DLA::MatrixDTranspose<typename TransposeTraits<T>::type>*>(this)->scatterAdd(M, Map, nMap);
      break;
    case LA::eDenseNonZeroPattern:
      static_cast<DLA::DenseNonZeroPattern<T>*>(this)->scatterAdd(M, Map, nMap);
      break;
    case LA::eSparseMatrix_CRS:
      static_cast<SLA::SparseMatrix_CRS<T>*>(this)->scatterAdd(M, Map, nMap);
      break;
    case LA::eSparseMatrix_CRS_Transpose:
      static_cast<SLA::SparseMatrix_CRS_Transpose<typename TransposeTraits<T>::type>*>(this)->scatterAdd(M, Map, nMap);
      break;
    case LA::eSparseNonZeroPattern:
      static_cast<SLA::SparseNonZeroPattern<T>*>(this)->scatterAdd(M, Map, nMap);
      break;
    case LA::eSparseNonZeroPattern_Transpose:
      static_cast<SLA::SparseNonZeroPattern_Transpose<typename TransposeTraits<T>::type>*>(this)->scatterAdd(M, Map, nMap);
      break;
    default:
      SANS_DEVELOPER_EXCEPTION("Unknown matrix type = %d", eMatrix_);
    }
  }

  void scatterAdd( const DLA::MatrixDView< T >& M, const int rowMap[], const int nrow, const int colMap[], const int ncol )
  {
    switch (eMatrix_)
    {
    case LA::eMatrixD:
      static_cast<DLA::MatrixD<T>*>(this)->scatterAdd(M, rowMap, nrow, colMap, ncol);
      break;
    case LA::eMatrixDTranspose:
      static_cast<DLA::MatrixDTranspose<typename TransposeTraits<T>::type>*>(this)->scatterAdd(M, rowMap, nrow, colMap, ncol);
      break;
    case LA::eDenseNonZeroPattern:
      static_cast<DLA::DenseNonZeroPattern<T>*>(this)->scatterAdd(M, rowMap, nrow, colMap, ncol);
      break;
    case LA::eSparseMatrix_CRS:
      static_cast<SLA::SparseMatrix_CRS<T>*>(this)->scatterAdd(M, rowMap, nrow, colMap, ncol);
      break;
    case LA::eSparseMatrix_CRS_Transpose:
      reinterpret_cast<SLA::SparseMatrix_CRS_Transpose<typename TransposeTraits<T>::type>*>(this)->scatterAdd(M, rowMap, nrow, colMap, ncol);
      break;
    case LA::eSparseNonZeroPattern:
      static_cast<SLA::SparseNonZeroPattern<T>*>(this)->scatterAdd(M, rowMap, nrow, colMap, ncol);
      break;
    case LA::eSparseNonZeroPattern_Transpose:
      static_cast<SLA::SparseNonZeroPattern_Transpose<typename TransposeTraits<T>::type>*>(this)->scatterAdd(M, rowMap, nrow, colMap, ncol);
      break;
    default:
      SANS_DEVELOPER_EXCEPTION("Unknown matrix type = %d", eMatrix_);
    }
  }

  // Overall matrix size
  int m() const { return m_; } // number of rows
  int n() const { return n_; } // number of columns

protected:

  explicit MatrixScatterAdd( const MatrixScatterAdd& M ) : m_(M.m_), n_(M.n_), eMatrix_(M.eMatrix_) {}
  MatrixScatterAdd& operator=( const MatrixScatterAdd& ) = delete;

  explicit MatrixScatterAdd( const LA::DerivedMatrix& eMatrix, const int m, const int n ) : m_(m), n_(n), eMatrix_(eMatrix) {}

  // private destructor as the base class in intended as an interface, and a matrix should never be deallocated using MatrixScatterAdd*
  ~MatrixScatterAdd() {}

  int m_; // number of rows
  int n_; // number of columns

private:
  // Type to indicate the derived type for this instance
  // This is private so a derived type can't accidentally modify it
  const LA::DerivedMatrix eMatrix_;
};


template<class T>
struct MatrixScatterAdd< DLA::MatrixD<T> >
{

  void scatterAdd( const DLA::MatrixDView< T >& M, const int row, const int col )
  {
    switch (eMatrix_)
    {
    case LA::eMatrixD:
      static_cast<DLA::MatrixD<DLA::MatrixD<T>>*>(this)->scatterAdd(M, row, col);
      break;
    case LA::eMatrixDTranspose:
      static_cast<DLA::MatrixDTranspose<typename TransposeTraits<DLA::MatrixD<T>>::type>*>(this)->scatterAdd(M, row, col);
      break;
    case LA::eDenseNonZeroPattern:
      static_cast<DLA::DenseNonZeroPattern<DLA::MatrixD<T>>*>(this)->scatterAdd(M, row, col);
      break;
    case LA::eSparseMatrix_CRS:
      static_cast<SLA::SparseMatrix_CRS<DLA::MatrixD<T>>*>(this)->scatterAdd(M, row, col);
      break;
    case LA::eSparseMatrix_CRS_Transpose:
      static_cast<SLA::SparseMatrix_CRS_Transpose<typename TransposeTraits<DLA::MatrixD<T>>::type>*>(this)->scatterAdd(M, row, col);
      break;
    case LA::eSparseNonZeroPattern:
      static_cast<SLA::SparseNonZeroPattern<DLA::MatrixD<T>>*>(this)->scatterAdd(M, row, col);
      break;
    case LA::eSparseNonZeroPattern_Transpose:
      static_cast<SLA::SparseNonZeroPattern_Transpose<typename TransposeTraits<DLA::MatrixD<T>>::type>*>(this)->scatterAdd(M, row, col);
      break;
    default:
      SANS_DEVELOPER_EXCEPTION("Unknown matrix type = %d", eMatrix_);
    }
  }

  // Overall matrix size
  int m() const { return m_; } // number of rows
  int n() const { return n_; } // number of columns

protected:

  explicit MatrixScatterAdd( const MatrixScatterAdd& M ) : m_(M.m_), n_(M.n_), eMatrix_(M.eMatrix_) {}
  MatrixScatterAdd& operator=( const MatrixScatterAdd& ) = delete;

  explicit MatrixScatterAdd( const LA::DerivedMatrix& eMatrix, const int m, const int n ) : m_(m), n_(n), eMatrix_(eMatrix) {}

  // private destructor as the base class in intended as an interface, and a matrix should never be deallocated using MatrixScatterAdd*
  ~MatrixScatterAdd() {}

  int m_; // number of rows
  int n_; // number of columns

private:
  // Type to indicate the derived type for this instance
  // This is private so a derived type can't accidentally modify it
  const LA::DerivedMatrix eMatrix_;
};

}

#endif //MATRIXSCATTERADD_H
