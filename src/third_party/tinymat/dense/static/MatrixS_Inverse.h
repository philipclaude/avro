// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXS_INVERSE_H
#define MATRIXS_INVERSE_H

#include "tools/SANSnumerics.h"     // Real
#include "MatrixS.h"
#include "MatrixS_Mul.h"

//This is a general interfaces for adding matrix inverse operations. These classes are specialized for
//the specific type of inverse through the SolverType template argument

namespace tinymat 
{
namespace DLA
{
namespace Fixed
{

  //This series of static functions generate the appropriate temporary variables if necessary
  template< template<int, int, class, class> class SolverType, int M, int N, class T >
  struct InverseSolver
  {
    // Inverse( A ) * (a + b)
    template< class RExpr, bool useRFR, bool FullR, class MatrixType >
    static void Solve( const typename SolverType<M, N, T, MatrixType>::FactorType& MatrixFac,
                       const MatrixSType< RExpr, useRFR, FullR >& eR, const Real sgn, MatrixType& res )
    {
      res = eR;
      Solve( MatrixFac, sgn, res );
    }

    // Inverse( A ) * b
    template< int MR, int NR, class MatrixType >
    static void Solve( const typename SolverType<M, N, T, MatrixType>::FactorType& MatrixFac,
                       const MatrixS< MR, NR, T >& MatrixR, const Real sgn, MatrixType& res )
    {
      res = MatrixR;
      Solve(MatrixFac, sgn, res);
    }

    // This calls that actual inverse algorithm
    template< class MatrixType >
    static void Solve( const typename SolverType<M, N, T, MatrixType>::FactorType& MatrixFac,
                       const Real sgn, MatrixType& res )
    {
      SolverType<M, N, T, MatrixType>::Solve(MatrixFac, sgn, res);
    }
  };

  //Class for representing the inverse of a matrix
  template< template<int, int, class, class> class SolverType, class Expr >
  struct MatrixInverse : public MatrixSType< MatrixInverse<SolverType, Expr>, true, true >
  {
    typedef typename Expr::Ttype Ttype;
    typedef Ttype T;
    static const int M = Expr::M;
    static const int N = Expr::N;
    typedef typename SolverType<M, N, T, MatrixS<M,N,T> >::FactorType FactorType;

    const Expr& MatrixExpr; //The matrix expression that is solved

    MatrixInverse( const Expr& MatrixExpr ) : MatrixExpr(MatrixExpr) {}

    // A copy of MatrixExpr is created so that Inverse( A ) does not modify 'A'
    template<class RExpr, bool useRFR, bool FullR>
    inline void value(const Real sgn, MatrixSType< RExpr, useRFR, FullR>& res) const
    {
      FactorType AFac(MatrixExpr);
      res.cast() = Identity();
      InverseSolver<SolverType, M, N, T>::Solve(AFac, sgn, res.cast());
    }
/*
    // This allows for the expression (Ainv, Binv) = Inverse( C );
    template< class MatrixL >
    inline void value(const Real sgn, DenseMatrixTuple<MatrixL>& res) const
    {
      DenseMatrix<T> A(MatrixExpr);
      res = Identity();
      InverseSolver<SolverType, T>::Solve(A, sgn, res);
    }
*/
    inline const MatrixInverse&
    operator+() const { return *this; }
  };


  // Forward declare
  template< template<int, int, class, class> class SolverType, int ML, int NL, class T, class RExpr>
  struct MatrixBackSolve;

  //Represents a factorized matrix
  template< template<int, int, class, class> class SolverType, int M, int N, class T>
  struct MatrixFactorized
  {
    typedef typename SolverType<M, N, T, MatrixS<M,N,T> >::FactorType FactorType;

    template<class Expr, bool useRF, bool Full>
    MatrixFactorized( const MatrixSType< Expr, useRF, Full >& Matrix)
      : MatrixFac(Matrix)
    {}

    // constructor for a no-op factorization
    explicit MatrixFactorized( const DLA::Identity& I)
      : MatrixFac(I)
    {}

    //This allows for the expression x = factor.backsolve(a + b)
    template<class RExpr, bool useRFR, bool FullR>
    inline MatrixBackSolve<SolverType, M, N, T, RExpr >
    backsolve(const MatrixSType< RExpr, useRFR, FullR>& MatrixRExpr) const
    {
      return MatrixBackSolve<SolverType, M, N, T, RExpr >(MatrixFac, MatrixRExpr.cast());
    }

  protected:
    FactorType MatrixFac;    //The factorized matrix
  };

  //Represents the solution of a matrix inverse multiplied with a matrix expression, i.e. Factor.backsolve(b);
  template< template<int, int, class, class> class SolverType, int ML, int NL, class T, class RExpr>
  struct MatrixBackSolve : public MatrixSType< MatrixBackSolve<SolverType, ML, NL, T, RExpr>, true, true >
  {
    typedef T Ttype;
    static const int M = NL; //This is particularly true for InverseQR
    static const int N = RExpr::N;
    typedef typename SolverType<ML, NL, T, MatrixS<M,N,T> >::FactorType FactorType;

    friend struct MatrixFactorized<SolverType, ML, NL, T>;

  protected:
    MatrixBackSolve( const FactorType& MatrixFac, const RExpr& MatrixRExpr )
      : MatrixFac(MatrixFac), MatrixRExpr(MatrixRExpr) {}

  public:
    //This allows for the general expression x = s*Inverse( A + B )*(a + b)
    template<class Tres>
    inline void value(const Real sgn, MatrixS<M, N, Tres>& res) const
    {
      InverseSolver<SolverType, ML, NL, T>::Solve(MatrixFac, MatrixRExpr, sgn, res);
    }
/*
    //This allows for the general expression (x1, x2) = s*Inverse( A + B )*(a1 + a2, b1 + b2)
    template< class MatrixL >
    inline void value(const Real sgn, MatrixDTuple<MatrixL>& res) const
    {
      InverseSolver<SolverType, T>::Solve(MatrixFac, MatrixRExpr, sgn, res);
    }
*/
    inline const MatrixBackSolve&
    operator+() const { return *this; }

  protected:
    const FactorType& MatrixFac; //The factorized matrix that is solved
    const RExpr& MatrixRExpr;    //The RHS matrix expression
  };


  //Represents the solution of a matrix inverse multiplied with a matrix expression, i.e. Solve(A,b);
  template< template<int,int,class, class> class SolverType, class InvExpr, class RExpr>
  struct MatrixSolve : public MatrixSType< MatrixSolve<SolverType, InvExpr, RExpr>, true, true >
  {
    typedef typename InvExpr::Ttype Ttype;
    typedef Ttype T;
    static const int M = InvExpr::N; //This is particularly true for InverseQR
    static const int N = RExpr::N;
    typedef typename SolverType<InvExpr::M, InvExpr::N, T, MatrixS<M,N,T> >::FactorType FactorType;

    template<class Expr, bool useRF, bool Full>
    MatrixSolve( const MatrixSType< Expr, useRF, Full >& Matrix, const RExpr& MatrixRExpr )
      : MatrixFac(Matrix), MatrixRExpr(MatrixRExpr) {}

    //This allows for the general expression x = s*Inverse( A + B )*(a + b)
    template<class Expr, bool useRF, bool Full>
    inline void value(const Real sgn, MatrixSType< Expr, useRF, Full>& res) const
    {
      InverseSolver<SolverType, InvExpr::M, InvExpr::N, T>::Solve(MatrixFac, MatrixRExpr, sgn, res.cast());
    }
/*
    //This allows for the general expression (x1, x2) = s*Inverse( A + B )*(a1 + a2, b1 + b2)
    template< class MatrixL >
    inline void value(const Real sgn, MatrixDTuple<MatrixL>& res) const
    {
      InverseSolver<SolverType, T>::Solve(MatrixFac, MatrixRExpr, sgn, res);
    }
*/
    inline const MatrixSolve&
    operator+() const { return *this; }

  protected:
    const FactorType MatrixFac; //The factorized matrix that is solved
    const RExpr& MatrixRExpr;   //The matrix expression multiplying the inverse operation
  };

#if 0
  //Represents the multiplication of a matrix inverse with a matrix expression, i.e. InverseLU(A)*b;
  template< template<int, class, class> class SolverType, class InvExpr, class RExpr>
  struct MatrixSolve : public MatrixSType< MatrixSolve<SolverType, InvExpr, RExpr>, true, true >
  {
    typedef typename InvExpr::Ttype Ttype;
    typedef Ttype T;
    static const int M = InvExpr::N; //This is particualry true for InverseQR
    static const int N = RExpr::N;

    const InvExpr& MatrixInvExpr; //The matrix expression that is solved
    const RExpr& MatrixRExpr;     //The matrix expression multiplying the inverse operation
    Real s;                       //A scalar quantity multiplying the inverse operation

    MatrixSolve( const InvExpr& MatrixInvExpr, const RExpr& MatrixRExpr, const Real s = 1 )
      : MatrixInvExpr(MatrixInvExpr), MatrixRExpr(MatrixRExpr), s(s) {}

    //This allows for the general expression x = s*Inverse( A + B )*(a + b)
    template<class Tres>
    inline void value(const Real sgn, MatrixS<M, N, Tres>& res) const
    {
      InverseSolver<SolverType, T>::Solve(MatrixInvExpr, MatrixRExpr, s*sgn, res);
    }
/*
    //This allows for the general expression (x1, x2) = s*Inverse( A + B )*(a1 + a2, b1 + b2)
    template< class MatrixL >
    inline void value(const Real sgn, DenseMatrixTuple<MatrixL>& res) const
    {
      InverseSolver<SolverType, T>::Solve(MatrixInvExpr, MatrixRExpr, s*sgn, res);
    }
*/
    inline const MatrixSolve&
    operator+() const { return *this; }
  };
#endif
} //namespace Fixed

//=============================================================================
//Preclude multiplication between an inverse-matrix and a matrix, i.e. InverseLU(A)*b;
template< template<int, int, class, class> class SolverType, class InvExpr, class RExpr, bool useRFR, bool FullR >
void
operator*(const Fixed::MatrixInverse< SolverType, InvExpr>& InvMatrix,
          const MatrixSType< RExpr, useRFR, FullR >& MatrixRExpr)
{
  //The syntax Inverse(A)*b is not permitted
}

//=============================================================================
//Preclude multiplication between a scalar, an inverse-matrix, and a matrix, i.e. 2*InverseLU(A)*b;
template< template<int, int, class, class> class SolverType, class InvExpr,
          class T, bool useRFInv, bool FullInv, class RExpr, bool useRFR, bool FullR >
void
operator*(const OpMulSScalar< Fixed::MatrixInverse< SolverType, InvExpr>, T, useRFInv, FullInv >& ScalMulInvMatrix,
          const MatrixSType< RExpr, useRFR, FullR >& MatrixRExpr)
{
  //The syntax Inverse(A)*b is not permitted
}


//=============================================================================
//Operator for multiplication between an inverse-matrix and a matrix, i.e. InverseLU(A)*b;
#if 0 // Deprecated
template< template<int, class, class> class SolverType, class InvExpr, class RExpr, bool useRFR, bool FullR >
inline Fixed::MatrixInverseMul< SolverType, InvExpr, RExpr >
operator*(const Fixed::MatrixInverse< SolverType, InvExpr>& InvMatrix,
          const MatrixSType< RExpr, useRFR, FullR >& MatrixRExpr)
{
  return Fixed::MatrixInverseMul< SolverType, InvExpr, RExpr >(
         InvMatrix.MatrixExpr, MatrixRExpr.cast());
}
#endif

//=============================================================================
//Operator for multiplication between a scalar, an inverse-matrix, and a matrix, i.e. 2*InverseLU(A)*b;
#if 0 // Deprecated
template< template<int, class, class> class SolverType, class InvExpr, class T, bool useRFInv, bool FullInv, class RExpr, bool useRFR, bool FullR >
inline Fixed::MatrixInverseMul< SolverType, InvExpr, RExpr >
operator*(const OpMulSScalar< Fixed::MatrixInverse< SolverType, InvExpr>, T, useRFInv, FullInv >& ScalMulInvMatrix,
          const MatrixSType< RExpr, useRFR, FullR >& MatrixRExpr)
{
  return Fixed::MatrixInverseMul< SolverType, InvExpr, RExpr >(
         ScalMulInvMatrix.e.MatrixExpr, MatrixRExpr.cast(), ScalMulInvMatrix.s);
}
#endif

} //namespace DLA
} //namespace tinymat 



#endif //MATRIXS_INVERSE_H
