#ifndef ALGEBRAICEQUATIONSET_TRAITS_H
#define ALGEBRAICEQUATIONSET_TRAITS_H

#include <numpack/AlgebraicEquationSetBase.h>

#include "numpack/sparse/SparseMatrix_CRS.h"
#include "numpack/sparse/SparseVector.h"
#include "numpack/sparse/SparseNonZeroPattern.h"

#include "numpack/dense/dynamic/MatrixD.h"
#include "numpack/dense/dynamic/VectorD.h"
#include "numpack/dense/dynamic/MatrixD_NonZeroPattern.h"

namespace avro
{

class AlgEqSetTraits_Sparse;
class AlgEqSetTraits_BlockSparse;
class AlgEqSetTraits_Dense;

template<class MatrixQ_, class ArrayQ_, class Type>
struct EquationSetTraits;

//---------------------------------------------------------------------------//
// Algebraic equation set traits for using Sparse matrix representation of the equation set

template<class MatrixQ_, class ArrayQ_>
struct EquationSetTraits<MatrixQ_, ArrayQ_, AlgEqSetTraits_Sparse>
{
  typedef ArrayQ_ ArrayQ;
  typedef MatrixQ_ MatrixQ;

  typedef numpack::SLA::SparseMatrix_CRS<MatrixQ>     MatrixClass;
  typedef numpack::SLA::SparseVector<ArrayQ>          VectorClass;
  typedef numpack::SLA::SparseNonZeroPattern<MatrixQ> NonZeroPatternClass;

  typedef numpack::DLA::MatrixD<MatrixClass> SystemMatrix;
  typedef numpack::DLA::VectorD<VectorClass> SystemVector;
  typedef numpack::DLA::MatrixD<NonZeroPatternClass> SystemNonZeroPattern;

  typedef numpack::DLA::MatrixDView<MatrixClass> SystemMatrixView;
  typedef numpack::DLA::VectorDView<VectorClass> SystemVectorView;
  typedef numpack::DLA::MatrixDView<NonZeroPatternClass> SystemNonZeroPatternView;

  typedef typename numpack::MatrixSizeType<SystemMatrix>::type MatrixSizeClass;
  typedef typename numpack::VectorSizeType<SystemVector>::type VectorSizeClass;

  template<class ArrayQT>
  using SystemVectorTemplate = numpack::DLA::VectorD< numpack::SLA::SparseVector<ArrayQT> >;

  typedef numpack::AlgebraicEquationSetBase<SystemMatrix> AlgebraicEquationSetBaseClass;
};

//---------------------------------------------------------------------------//
// Algebraic equation set traits for using block dynamic Sparse matrix representation of the equation set

template<class MatrixQ_, class ArrayQ_>
struct EquationSetTraits<MatrixQ_, ArrayQ_, AlgEqSetTraits_BlockSparse>
{
  typedef ArrayQ_ ArrayQ;
  typedef MatrixQ_ MatrixQ;

  typedef numpack::SLA::SparseMatrix_CRS< numpack::DLA::MatrixD<MatrixQ> >     MatrixClass;
  typedef numpack::SLA::SparseVector< numpack::DLA::VectorD<ArrayQ> >          VectorClass;
  typedef numpack::SLA::SparseNonZeroPattern< numpack::DLA::MatrixD<MatrixQ> > NonZeroPatternClass;

  typedef numpack::DLA::MatrixD<MatrixClass> SystemMatrix;
  typedef numpack::DLA::VectorD<VectorClass> SystemVector;
  typedef numpack::DLA::MatrixD<NonZeroPatternClass> SystemNonZeroPattern;

  typedef numpack::DLA::MatrixDView<MatrixClass> SystemMatrixView;
  typedef numpack::DLA::VectorDView<VectorClass> SystemVectorView;
  typedef numpack::DLA::MatrixDView<NonZeroPatternClass> SystemNonZeroPatternView;

  typedef typename numpack::MatrixSizeType<SystemMatrix>::type MatrixSizeClass;
  typedef typename numpack::VectorSizeType<SystemVector>::type VectorSizeClass;

  template<class ArrayQT>
  using SystemVectorTemplate = numpack::DLA::VectorD< numpack::SLA::SparseVector<ArrayQT> >;

  typedef numpack::AlgebraicEquationSetBase<SystemMatrix> EquationSetBaseClass;
};

//---------------------------------------------------------------------------//
// Algebraic equation set traits for using Dense matrix representation of the equation set

template<class MatrixQ_, class ArrayQ_>
struct EquationSetTraits<MatrixQ_, ArrayQ_, AlgEqSetTraits_Dense>
{
  typedef ArrayQ_ ArrayQ;
  typedef MatrixQ_ MatrixQ;

  typedef numpack::DLA::MatrixD<MatrixQ>             MatrixClass;
  typedef numpack::DLA::VectorD<ArrayQ>              VectorClass;
  typedef numpack::DLA::DenseNonZeroPattern<MatrixQ> NonZeroPatternClass;

  typedef numpack::DLA::MatrixD<MatrixClass> SystemMatrix;
  typedef numpack::DLA::VectorD<VectorClass> SystemVector;
  typedef numpack::DLA::MatrixD<NonZeroPatternClass> SystemNonZeroPattern;

  typedef numpack::DLA::MatrixDView<MatrixClass> SystemMatrixView;
  typedef numpack::DLA::VectorDView<VectorClass> SystemVectorView;
  typedef numpack::DLA::MatrixDView<NonZeroPatternClass> SystemNonZeroPatternView;

  typedef numpack::DLA::MatrixD< numpack::DLA::DenseMatrixSize > MatrixSizeClass;
  typedef numpack::DLA::VectorD< numpack::DLA::DenseVectorSize > VectorSizeClass;

  template<class ArrayQT>
  using SystemVectorTemplate = numpack::DLA::VectorD< numpack::DLA::VectorD<ArrayQT> >;

  typedef numpack::AlgebraicEquationSetBase<SystemMatrix> EquationSetBaseClass;
};

} // avro

#endif //ALGEBRAICEQUATIONSET_TRAITS_H
