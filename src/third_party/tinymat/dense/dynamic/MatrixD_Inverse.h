// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXD_INVERSE_H
#define MATRIXD_INVERSE_H

#include "tools/SANSnumerics.h"     // Real
#include "MatrixD.h"

//This is a general interfaces for adding matrix inverse operations. These classes are specialized for
//the specific type of inverse through the SolverType template argument

namespace tinymat 
{
namespace DLA
{

//This series of static functions generate the appropriate temporary variables if necessary
template< template<class, class> class SolverType, class T >
struct InverseSolver
{
  typedef typename SolverType<T, MatrixD<T> >::FactorType FactorType;

  // Inverse( A ) * (a + b)
  template< class RExpr, bool useRFR, class MatrixType >
  static void Solve( const FactorType& MatrixFac, const MatrixDType< RExpr, useRFR >& eR, const Real sgn, MatrixType& res )
  {
    res = eR;
    Solve( MatrixFac, sgn, res );
  }

  // Inverse( A ) * b
  template< class MatrixType >
  static void Solve( const FactorType& MatrixFac, const MatrixDView< T >& MatrixR, const Real sgn, MatrixType& res )
  {
    if ( &MatrixR(0,0) != &res(0,0)) res = MatrixR;
    Solve(MatrixFac, sgn, res);
  }

  // This calls that actual inverse algorithm
  template< class MatrixType >
  static void Solve( const FactorType& MatrixFac, const Real sgn, MatrixType& res )
  {
    SolverType<T, MatrixType>::Solve(MatrixFac, sgn, res);
  }
};

//Class for representing the inverse of a matrix
template< template<class, class> class SolverType, class Expr>
struct MatrixInverse : public MatrixDType< MatrixInverse<SolverType, Expr>, true >
{
  typedef typename Expr::node_type node_type;
  typedef node_type T;

  typedef typename SolverType<T, MatrixD<T> >::FactorType FactorType;

  const Expr& MatrixExpr; //The matrix expression that is solved

  MatrixInverse( const Expr& MatrixExpr ) : MatrixExpr(MatrixExpr) {}

  // A copy of MatrixExpr is created so that Inverse( A ) does not modify 'A'
  template<class Tres>
  inline void value(const Real sgn, MatrixDView<Tres>& res) const
  {
    FactorType AFac(MatrixExpr);
    res = Identity();
    InverseSolver<SolverType, T>::Solve(AFac, sgn, res);
  }

  // This allows for the expression (Ainv, Binv) = Inverse( C );
  template< class MatrixL >
  inline void value(const Real sgn, MatrixDTuple<MatrixL>& res) const
  {
    FactorType AFac(MatrixExpr);
    res = Identity();
    InverseSolver<SolverType, T>::Solve(AFac, sgn, res);
  }

  inline const MatrixInverse&
  operator+() const { return *this; }
  int m() const { return MatrixExpr.m(); }
  int n() const { return MatrixExpr.n(); }
  int size() const { return m()*n(); }
};

// Forward declare
template< template<class, class> class SolverType, class T, class RExpr>
struct MatrixBackSolve;

//Represents a factorized matrix
template< template<class, class> class SolverType,  class T>
struct MatrixFactorized
{
  typedef typename SolverType<T, MatrixD<T> >::FactorType FactorType;

  template<class Expr, bool useRF>
  MatrixFactorized( const MatrixDType< Expr, useRF >& Matrix)
    : MatrixFac(Matrix)
  {}

  //This allows for the expression x = factor.backsolve(a + b)
  template<class RExpr, bool useRFR>
  inline MatrixBackSolve<SolverType, T, RExpr >
  backsolve(const MatrixDType< RExpr, useRFR>& MatrixRExpr) const
  {
    return MatrixBackSolve<SolverType, T, RExpr >(MatrixFac, MatrixRExpr.cast());
  }

protected:
  FactorType MatrixFac;    //The factorized matrix
};

//Represents the solution of a matrix inverse multiplied with a matrix expression, i.e. Factor.backsolve(b);
template< template<class, class> class SolverType, class T, class RExpr>
struct MatrixBackSolve : public MatrixDType< MatrixBackSolve<SolverType, T, RExpr>, true >
{
  typedef T node_type;
  typedef typename SolverType<T, MatrixD<T> >::FactorType FactorType;

  friend struct MatrixFactorized<SolverType,T>;

protected:
  MatrixBackSolve( const FactorType& MatrixFac, const RExpr& MatrixRExpr )
    : MatrixFac(MatrixFac), MatrixRExpr(MatrixRExpr) {}

public:
  //This allows for the general expression x = s*Inverse( A + B )*(a + b)
  template<class Tres>
  inline void value(const Real sgn, MatrixDView<Tres>& res) const
  {
    InverseSolver<SolverType, T>::Solve(MatrixFac, MatrixRExpr, sgn, res);
  }

  //This allows for the general expression (x1, x2) = s*Inverse( A + B )*(a1 + a2, b1 + b2)
  template< class MatrixL >
  inline void value(const Real sgn, MatrixDTuple<MatrixL>& res) const
  {
    InverseSolver<SolverType, T>::Solve(MatrixFac, MatrixRExpr, sgn, res);
  }

  inline const MatrixBackSolve&
  operator+() const { return *this; }
  int m() const { return MatrixFac.m(); }
  int n() const { return MatrixRExpr.n(); }
  int size() const { return m()*n(); }

protected:
  const FactorType& MatrixFac; //The factorized matrix that is solved
  const RExpr& MatrixRExpr;    //The RHS matrix expression
};


//Represents the solution of a matrix inverse multiplied with a matrix expression, i.e. Solve(A,b);
template< template<class, class> class SolverType, class T, class RExpr>
struct MatrixSolve : public MatrixDType< MatrixSolve<SolverType, T, RExpr>, true >
{
  typedef T node_type;
  typedef typename SolverType<T, MatrixD<T> >::FactorType FactorType;

  template<class Expr, bool useRF>
  MatrixSolve( const MatrixDType< Expr, useRF >& Matrix, const RExpr& MatrixRExpr )
    : MatrixFac(Matrix), MatrixRExpr(MatrixRExpr) {}

public:
  //This allows for the general expression x = s*Inverse( A + B )*(a + b)
  template<class Tres>
  inline void value(const Real sgn, MatrixDView<Tres>& res) const
  {
    InverseSolver<SolverType, T>::Solve(MatrixFac, MatrixRExpr, sgn, res);
  }

  //This allows for the general expression (x1, x2) = s*Inverse( A + B )*(a1 + a2, b1 + b2)
  template< class MatrixL >
  inline void value(const Real sgn, MatrixDTuple<MatrixL>& res) const
  {
    InverseSolver<SolverType, T>::Solve(MatrixFac, MatrixRExpr, sgn, res);
  }

  inline const MatrixSolve&
  operator+() const { return *this; }
  int m() const { return MatrixFac.m(); }
  int n() const { return MatrixRExpr.n(); }
  int size() const { return m()*n(); }

protected:
  const FactorType MatrixFac; //The factorized matrix that is solved
  const RExpr& MatrixRExpr;   //The matrix expression multiplying the inverse operation
};


//=============================================================================
//Preclude multiplication between an inverse-matrix and a matrix, i.e. InverseLU(A)*b;
template< template<class, class> class SolverType, class InvExpr, class RExpr, bool useRFR >
void
operator*(const MatrixInverse< SolverType, InvExpr>& InvMatrix,
          const MatrixDType< RExpr, useRFR >& MatrixRExpr)
{
  //The syntax Inverse(A)*b is not permitted
}

//=============================================================================
//Preclude multiplication between a scalar, an inverse-matrix, and a matrix, i.e. 2*InverseLU(A)*b;
template< template<class, class> class SolverType, class InvExpr, bool useRFInv, class RExpr, bool useRFR >
void
operator*(const OpMulDScalar< MatrixInverse< SolverType, InvExpr>, useRFInv >& ScalMulInvMatrix,
          const MatrixDType< RExpr, useRFR >& MatrixRExpr)
{
  //The syntax Inverse(A)*b is not permitted
}

//=============================================================================
//Operator for multiplication between an inverse-matrix and a matrix, i.e. InverseLU(A)*b;
#if 0 //Deprecated
template< template<class, class> class SolverType, class InvExpr, class RExpr, bool useRFR >
inline MatrixInverseMul< SolverType, InvExpr, RExpr >
operator*(const MatrixInverse< SolverType, InvExpr>& InvMatrix,
          const MatrixDType< RExpr, useRFR >& MatrixRExpr)
{
  return MatrixInverseMul< SolverType, InvExpr, RExpr >(
         InvMatrix.MatrixExpr, MatrixRExpr.cast());
}
#endif

//=============================================================================
//Operator for multiplication between a scalar, an inverse-matrix, and a matrix, i.e. 2*InverseLU(A)*b;
#if 0 //Deprecated
template< template<class, class> class SolverType, class InvExpr, bool useRFInv, class RExpr, bool useRFR >
inline MatrixInverseMul< SolverType, InvExpr, RExpr >
operator*(const OpMulDScalar< MatrixInverse< SolverType, InvExpr>, useRFInv >& ScalMulInvMatrix,
          const MatrixDType< RExpr, useRFR >& MatrixRExpr)
{
  return MatrixInverseMul< SolverType, InvExpr, RExpr >(
         ScalMulInvMatrix.e.MatrixExpr, MatrixRExpr.cast(), ScalMulInvMatrix.s);
}
#endif

} //namespace DLA
} //namespace tinymat 



#endif //MATRIXD_INVERSE_H
