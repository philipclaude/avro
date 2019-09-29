// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef ALGEBRAICEQUATIONSET_TRAITS_H
#define ALGEBRAICEQUATIONSET_TRAITS_H

#include "AlgebraicEquationSetBase.h"

#include "numpack/sparse/SparseMatrix_CRS.h"
#include "numpack/sparse/SparseVector.h"
#include "numpack/sparse/SparseNonZeroPattern.h"

#include "numpack/dense/dynamic/MatrixD.h"
#include "numpack/dense/dynamic/VectorD.h"
#include "numpack/dense/dynamic/MatrixD_NonZeroPattern.h"

namespace numpack 
{

class AlgEqSetTraits_Sparse;
class AlgEqSetTraits_BlockSparse;
class AlgEqSetTraits_Dense;

template<class MatrixQ_, class ArrayQ_, class Type>
struct AlgebraicEquationSetTraits;

//---------------------------------------------------------------------------//
// Algebraic equation set traits for using Sparse matrix representation of the equation set

template<class MatrixQ_, class ArrayQ_>
struct AlgebraicEquationSetTraits<MatrixQ_, ArrayQ_, AlgEqSetTraits_Sparse>
{
  typedef ArrayQ_ ArrayQ;
  typedef MatrixQ_ MatrixQ;

  typedef SLA::SparseMatrix_CRS<MatrixQ>     MatrixClass;
  typedef SLA::SparseVector<ArrayQ>          VectorClass;
  typedef SLA::SparseNonZeroPattern<MatrixQ> NonZeroPatternClass;

  typedef DLA::MatrixD<MatrixClass> SystemMatrix;
  typedef DLA::VectorD<VectorClass> SystemVector;
  typedef DLA::MatrixD<NonZeroPatternClass> SystemNonZeroPattern;

  typedef DLA::MatrixDView<MatrixClass> SystemMatrixView;
  typedef DLA::VectorDView<VectorClass> SystemVectorView;
  typedef DLA::MatrixDView<NonZeroPatternClass> SystemNonZeroPatternView;

  typedef typename MatrixSizeType<SystemMatrix>::type MatrixSizeClass;
  typedef typename VectorSizeType<SystemVector>::type VectorSizeClass;

  template<class ArrayQT>
  using SystemVectorTemplate = DLA::VectorD< SLA::SparseVector<ArrayQT> >;

  typedef AlgebraicEquationSetBase<SystemMatrix> AlgebraicEquationSetBaseClass;
};

//---------------------------------------------------------------------------//
// Algebraic equation set traits for using block dynamic Sparse matrix representation of the equation set

template<class MatrixQ_, class ArrayQ_>
struct AlgebraicEquationSetTraits<MatrixQ_, ArrayQ_, AlgEqSetTraits_BlockSparse>
{
  typedef ArrayQ_ ArrayQ;
  typedef MatrixQ_ MatrixQ;

  typedef SLA::SparseMatrix_CRS< DLA::MatrixD<MatrixQ> >     MatrixClass;
  typedef SLA::SparseVector< DLA::VectorD<ArrayQ> >          VectorClass;
  typedef SLA::SparseNonZeroPattern< DLA::MatrixD<MatrixQ> > NonZeroPatternClass;

  typedef DLA::MatrixD<MatrixClass> SystemMatrix;
  typedef DLA::VectorD<VectorClass> SystemVector;
  typedef DLA::MatrixD<NonZeroPatternClass> SystemNonZeroPattern;

  typedef DLA::MatrixDView<MatrixClass> SystemMatrixView;
  typedef DLA::VectorDView<VectorClass> SystemVectorView;
  typedef DLA::MatrixDView<NonZeroPatternClass> SystemNonZeroPatternView;

  typedef typename MatrixSizeType<SystemMatrix>::type MatrixSizeClass;
  typedef typename VectorSizeType<SystemVector>::type VectorSizeClass;

  template<class ArrayQT>
  using SystemVectorTemplate = DLA::VectorD< SLA::SparseVector<ArrayQT> >;

  typedef AlgebraicEquationSetBase<SystemMatrix> AlgebraicEquationSetBaseClass;
};

//---------------------------------------------------------------------------//
// Algebraic equation set traits for using Dense matrix representation of the equation set

template<class MatrixQ_, class ArrayQ_>
struct AlgebraicEquationSetTraits<MatrixQ_, ArrayQ_, AlgEqSetTraits_Dense>
{
  typedef ArrayQ_ ArrayQ;
  typedef MatrixQ_ MatrixQ;

  typedef DLA::MatrixD<MatrixQ>             MatrixClass;
  typedef DLA::VectorD<ArrayQ>              VectorClass;
  typedef DLA::DenseNonZeroPattern<MatrixQ> NonZeroPatternClass;

  typedef DLA::MatrixD<MatrixClass> SystemMatrix;
  typedef DLA::VectorD<VectorClass> SystemVector;
  typedef DLA::MatrixD<NonZeroPatternClass> SystemNonZeroPattern;

  typedef DLA::MatrixDView<MatrixClass> SystemMatrixView;
  typedef DLA::VectorDView<VectorClass> SystemVectorView;
  typedef DLA::MatrixDView<NonZeroPatternClass> SystemNonZeroPatternView;

  typedef DLA::MatrixD< DLA::DenseMatrixSize > MatrixSizeClass;
  typedef DLA::VectorD< DLA::DenseVectorSize > VectorSizeClass;

  template<class ArrayQT>
  using SystemVectorTemplate = DLA::VectorD< DLA::VectorD<ArrayQT> >;

  typedef AlgebraicEquationSetBase<SystemMatrix> AlgebraicEquationSetBaseClass;
};

}

#endif //ALGEBRAICEQUATIONSET_TRAITS_H
