// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef LU_SGS_H
#define LU_SGS_H

#include "Python/PyDict.h"

#include "tools/noncopyable.h"

#include "numpack/AlgebraicEquationSetBase.h"
#include "numpack/sparse/LinearSolverBase.h"
#include "numpack/sparse/sparse_Inverse.h"
#include "numpack/sparse/SparseMatrix_Diag.h"

#include <vector>

namespace numpack 
{
namespace SLA
{


//=============================================================================
struct LU_SGSParam : noncopyable
{
  static void checkInputs(PyDict d);
  static LU_SGSParam params;
};


//=============================================================================
//
// Lower-Upper Symmetric Gauss-Seidel (LU-SGS)
//
// The matrix is decomposed into lower, L, diagonal, D, and upper, U, matrices as
//
//   A = L + D + U
//
// The LU-SGS preconditioner approximates A as
//
// M = (L + D)D^1(D + U) = L + D + U + L(D^-1)U
//
// Hence, M is a good approximation so long as L(D^-1)U is small.
//
// M^-1 is computed by solving the two system of equations
//
//   (L + D) x' = b
//   (I + D^-1 U) x = x'
//
template< class Matrix_type >
class LU_SGS : public LinearSolverBase< Matrix_type >
{
public:
  typedef LinearSolverBase< Matrix_type > Base_type;

  typedef typename Base_type::SparseVector_type SparseVector_type;
  typedef typename Base_type::SparseVectorView_type SparseVectorView_type;
  typedef typename NonZeroPatternType<Matrix_type>::type SystemNonZeroPattern;

  typedef typename Matrix_type::Ttype TM;
  typedef typename SparseVector_type::Ttype TV;

//-----------------------------------------------------------------------------
  LU_SGS( AlgebraicEquationSetBase<Matrix_type>& f, LinearSystemSolve solve = RegularSolve ) :
    Base_type(solve), f_(f)
  {
    init();
  }

  LU_SGS( const PyDict& d, AlgebraicEquationSetBase<Matrix_type>& f, LinearSystemSolve solve = RegularSolve ) :
    Base_type(solve), f_(f)
  {
    init();
  }

  virtual ~LU_SGS() { delete A_; }

//-----------------------------------------------------------------------------
  virtual void factorize() override;

//-----------------------------------------------------------------------------
  virtual LinearSolveStatus backsolve( const SparseVectorView_type& b, SparseVectorView_type& x ) const override;

protected:
  void init();

  using Base_type::A_;
  using Base_type::transpose_;
  AlgebraicEquationSetBase<Matrix_type>& f_;
  SparseMatrix_Diag<TM> Dinv_;
};



//=============================================================================
template< class Matrix_type >
class LU_SGS< DLA::MatrixD<Matrix_type> > : public LinearSolverBase< DLA::MatrixD<Matrix_type> >
{
public:
  typedef LinearSolverBase< DLA::MatrixD<Matrix_type> > Base_type;

  typedef typename Base_type::SparseVector_type SparseVector_type;
  typedef typename Base_type::SparseVectorView_type SparseVectorView_type;
  typedef typename NonZeroPatternType<DLA::MatrixD<Matrix_type>>::type SystemNonZeroPattern;

  typedef typename Matrix_type::Ttype TM;
  typedef typename SparseVector_type::node_type::Ttype TV;

//-----------------------------------------------------------------------------
  LU_SGS( AlgebraicEquationSetBase<DLA::MatrixD<Matrix_type>>& f, LinearSystemSolve solve = RegularSolve ) :
    Base_type(solve), f_(f)
  {
    init();
  }
  LU_SGS( const PyDict& d, AlgebraicEquationSetBase<DLA::MatrixD<Matrix_type>>& f, LinearSystemSolve solve = RegularSolve ) :
    Base_type(solve), f_(f)
  {
    init();
  }
  virtual ~LU_SGS() { delete A_; }

//-----------------------------------------------------------------------------
  virtual void factorize() override;

//-----------------------------------------------------------------------------
  virtual LinearSolveStatus backsolve( const SparseVectorView_type& b, SparseVectorView_type& x ) const override;

protected:
  void init();

  using Base_type::A_;
  using Base_type::transpose_;
  AlgebraicEquationSetBase<DLA::MatrixD<Matrix_type>>& f_;
  std::vector< SparseMatrix_Diag<TM> > Dinv_;
};


//typedef BLA::MatrixBlock_2x2<DLA::MatrixD<SLA::SparseMatrix_CRS<Real> >, DLA::MatrixD<SLA::SparseMatrix_CRS<Real> >,
//                             DLA::MatrixD<SLA::SparseMatrix_CRS<Real> >, DLA::MatrixD<SLA::SparseMatrix_CRS<Real> > > BlockMatrixRealReal;
//
//template< >
//class LU_SGS< BlockMatrixRealReal > : public LinearSolverBase< BlockMatrixRealReal >
//{
//public:
//  typedef LinearSolverBase< BlockMatrixRealReal > Base_type;
//
//  typedef typename Base_type::SparseVector_type SparseVector_type;
//  typedef typename Base_type::SparseVectorView_type SparseVectorView_type;
//
////-----------------------------------------------------------------------------
//  LU_SGS() {}
//  explicit LU_SGS( PyDict d )
//  {
//    SANS_DEVELOPER_EXCEPTION("LU_SGS<BlockMatrixRealReal> not implemented");
//  }
//  virtual ~LU_SGS() {}
//
////-----------------------------------------------------------------------------
//  virtual void factorize( const BlockMatrixRealReal& A ) override
//  {
//    SANS_DEVELOPER_EXCEPTION("LU_SGS<BlockMatrixRealReal> not implemented");
//  }
//
////-----------------------------------------------------------------------------
//  virtual LinearSolveStatus backsolve( const SparseVectorView_type& b, SparseVectorView_type& x ) const override
//  {
//    SANS_DEVELOPER_EXCEPTION("LU_SGS<BlockMatrixRealReal> not implemented");
//  }
//
//protected:
//};
//
//
//typedef BLA::MatrixBlock_2x2< DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<8,8,Real>> >,
//                              DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<8,8,Real>> >,
//                              DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<8,8,Real>> >,
//                              DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<8,8,Real>> > > BlockMatrix88;
//
//template< >
//class LU_SGS< BlockMatrix88 > : public LinearSolverBase< BlockMatrix88 >
//{
//public:
//  typedef LinearSolverBase< BlockMatrix88 > Base_type;
//
//  typedef typename Base_type::SparseVector_type SparseVector_type;
//  typedef typename Base_type::SparseVectorView_type SparseVectorView_type;
//
////-----------------------------------------------------------------------------
//  LU_SGS() {}
//  explicit LU_SGS( PyDict d )
//  {
//    SANS_DEVELOPER_EXCEPTION("LU_SGS<BlockMatrix88> not implemented");
//  }
//  virtual ~LU_SGS() {}
//
////-----------------------------------------------------------------------------
//  virtual void factorize( const BlockMatrix88& A ) override
//  {
//    SANS_DEVELOPER_EXCEPTION("LU_SGS<BlockMatrix88> not implemented");
//  }
//
////-----------------------------------------------------------------------------
//  virtual LinearSolveStatus backsolve( const SparseVectorView_type& b, SparseVectorView_type& x ) const override
//  {
//    SANS_DEVELOPER_EXCEPTION("LU_SGS<BlockMatrix88> not implemented");
//  }
//
//protected:
//};

template <class M00, class M01, class M10, class M11>
class LU_SGS< BLA::MatrixBlock_2x2<M00, M01, M10, M11> > :
  public LinearSolverBase< BLA::MatrixBlock_2x2<M00, M01, M10, M11> >
{
public:
  typedef BLA::MatrixBlock_2x2<M00, M01, M10, M11> Matrix_type;
  typedef LinearSolverBase< Matrix_type > Base_type;

  typedef typename Base_type::SparseVector_type SparseVector_type;
  typedef typename Base_type::SparseVectorView_type SparseVectorView_type;
  typedef typename NonZeroPatternType<Matrix_type>::type SystemNonZeroPattern;

//-----------------------------------------------------------------------------
  LU_SGS(AlgebraicEquationSetBase<Matrix_type>& f, LinearSystemSolve solve = RegularSolve) : Base_type(solve)
  {
    SANS_DEVELOPER_EXCEPTION("LU_SGS<BlockMatrix_2x2> not implemented");
  }
  explicit LU_SGS( const PyDict& d ) : Base_type(RegularSolve)
  {
    SANS_DEVELOPER_EXCEPTION("LU_SGS<BlockMatrix_2x2> not implemented");
  }
  virtual ~LU_SGS() {}

//-----------------------------------------------------------------------------
  virtual void factorize() override
  {
    SANS_DEVELOPER_EXCEPTION("LU_SGS<BlockMatrix_2x2> not implemented");
  }

//-----------------------------------------------------------------------------
  virtual LinearSolveStatus backsolve( const SparseVectorView_type& b, SparseVectorView_type& x ) const override
  {
    SANS_DEVELOPER_EXCEPTION("LU_SGS<BlockMatrix_2x2> not implemented");

    return LinearSolveStatus();
  }

protected:
};


template<class M00, class M01, class M02, class M03,
         class M10, class M11, class M12, class M13,
         class M20, class M21, class M22, class M23,
         class M30, class M31, class M32, class M33>
class LU_SGS< BLA::MatrixBlock_4x4<M00, M01, M02, M03,
                                   M10, M11, M12, M13,
                                   M20, M21, M22, M23,
                                   M30, M31, M32, M33> > :
  public LinearSolverBase< BLA::MatrixBlock_4x4<M00, M01, M02, M03,
                                                M10, M11, M12, M13,
                                                M20, M21, M22, M23,
                                                M30, M31, M32, M33> >
{
public:
  typedef BLA::MatrixBlock_4x4<M00, M01, M02, M03,
                               M10, M11, M12, M13,
                               M20, M21, M22, M23,
                               M30, M31, M32, M33> Matrix_type;
  typedef LinearSolverBase< Matrix_type > Base_type;

  typedef typename Base_type::SparseVector_type SparseVector_type;
  typedef typename Base_type::SparseVectorView_type SparseVectorView_type;
  typedef typename NonZeroPatternType<Matrix_type>::type SystemNonZeroPattern;

//-----------------------------------------------------------------------------
  LU_SGS(AlgebraicEquationSetBase<Matrix_type>& f, LinearSystemSolve solve = RegularSolve) : Base_type(solve)
  {
    SANS_DEVELOPER_EXCEPTION("LU_SGS<BlockMatrix_4x4> not implemented");
  }
  explicit LU_SGS( const PyDict& d ) : Base_type(RegularSolve)
  {
    SANS_DEVELOPER_EXCEPTION("LU_SGS<BlockMatrix_4x4> not implemented");
  }
  virtual ~LU_SGS() {}

//-----------------------------------------------------------------------------
  virtual void factorize() override
  {
    SANS_DEVELOPER_EXCEPTION("LU_SGS<BlockMatrix_4x4> not implemented");
  }

//-----------------------------------------------------------------------------
  virtual LinearSolveStatus backsolve( const SparseVectorView_type& b, SparseVectorView_type& x ) const override
  {
    SANS_DEVELOPER_EXCEPTION("LU_SGS<BlockMatrix_4x4> not implemented");

    return LinearSolveStatus();
  }

protected:
};

} //namespace SLA
} //namespace numpack 

#endif //LU_SGS_H
