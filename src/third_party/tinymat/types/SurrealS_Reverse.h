// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef SURREALS_REVERSE_H
#define SURREALS_REVERSE_H

//  overloaded derivative operator
//  ref: derivify.h (Google it)

#include <cmath>
#include <iostream>
#include <string>

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/type_traits/is_base_and_derived.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/not.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/range_c.hpp>


#include "tools/SANSnumerics.h"     // Real
#include "tools/SANSException.h"
#include "tools/AlignMem.h"
#include "tools/always_inline.h"

class SurrealSTypeBase {};

template< class Derived >
struct SurrealSType : SurrealSTypeBase
{
  //A convenient method for casting to the derived type
  ALWAYS_INLINE const Derived& cast() const { return static_cast<const Derived&>(*this); }

  //A simple way to call value without having to case first
  ALWAYS_INLINE Real value() const { return cast().value(); }
};

//Forward declarations
template<int N>
class SurrealS;

//#define SURREALS_LOOP_UNROLL

namespace SurrealSExpr
{
template<class L, class R > class OpMul;
class OpMul_impl;
}


//----------------------------------------------------------------------------//
// SurrealS:  value, N derivatives
//
// Operators with Lazy Expressions
//
// statically defined derivative array
//----------------------------------------------------------------------------//
template<int N_>
class SurrealS : public SurrealSType< SurrealS<N_> >
{
public:
  static const int N = N_;

  //The default constructor is intentionally empty here. This means Surreal is
  //not initialized when declared, which is consistent with regular numbers. This also
  //improves performance.
  ALWAYS_INLINE SurrealS() {}
  SurrealS( const SurrealS& z );
  SurrealS( const Real v0 );
  SurrealS( const Real v0, const Real d0[], int n );
  SurrealS( const Real v0, const Real& d0 );
  template<class Expr>
  ALWAYS_INLINE
  SurrealS( const SurrealSType<Expr>& r ) : v_(0) { operator=(r); }
  ALWAYS_INLINE ~SurrealS() {}

  ALWAYS_INLINE int size() const { return N; }

  // value accessor operators
  ALWAYS_INLINE       Real& value()       { return v_; }
  ALWAYS_INLINE const Real& value() const { return v_; }

  // derivative accessor operators
  ALWAYS_INLINE       Real& deriv( int i=0 )       { return d_[i]; }
  ALWAYS_INLINE const Real& deriv( int i=0 ) const { return d_[i]; }

  //Data and methods for reverse differentiation
  static const int nArgs = 1;
  void getPartials( const Real& bar, Real partials[] ) const
  {
    partials[0] = bar;
  }
  template<int arg>
  const Real& getArgDeriv(const int i) const
  {
    return d_[i];
  }

  // assignment
  SurrealS& operator=( const SurrealS& );
  SurrealS& operator=( const Real& );

  template<class Expr> SurrealS& operator= ( const SurrealSType<Expr>& );
  template<class Expr> SurrealS& operator+=( const SurrealSType<Expr>& );
  template<class Expr> SurrealS& operator-=( const SurrealSType<Expr>& );

  // unary operators; no side effects
  const SurrealS& operator+() const;

  // binary accumulation operators
  SurrealS& operator+=( const Real& );
  SurrealS& operator-=( const Real& );
  SurrealS& operator*=( const SurrealS& );
  SurrealS& operator*=( const Real& );
  SurrealS& operator/=( const SurrealS& );
  SurrealS& operator/=( const Real& );

#if 0
  // classification functions <cmath>
  friend bool isfinite( const SurrealS& );
  friend bool isinf( const SurrealS& );
  friend bool isnan( const SurrealS& );
#endif

  // input/output
  template<int M> friend std::istream& operator>>( std::istream&, SurrealS<M>& );

protected:
  ALIGN_MEM Real d_[N];   // derivative array
  ALIGN_MEM Real v_;      // value


  // Functor for boost::mpl::for_each to compute the accumulation
  template <typename Expr>
  struct AccumOp
  {
    static const int nArgs = Expr::nArgs;
    const Expr& e_;
    Real& sum;
    Real partials_[nArgs];
    int i;
    inline AccumOp(const Expr& e, Real& sum) : e_(e), sum(sum), i(0) { e.getPartials(1., partials_); }
    inline void getDeriv(int i) { this->i = i; sum = 0; }
    template <typename ArgT>
    inline void operator () (ArgT arg) const
    {
      const int Arg = ArgT::value;
      sum += partials_[Arg] * e_.template getArgDeriv<Arg>(i);
    }
  };

};

//Constructors

template<int N>
ALWAYS_INLINE
SurrealS<N>::SurrealS( const SurrealS& z )
{
  v_ = z.v_;
  for (int i = 0; i < N; i++)
    d_[i] = z.d_[i];
}
template<int N>
ALWAYS_INLINE
SurrealS<N>::SurrealS( const Real v0 )
{
  v_ = v0;
  for (int i = 0; i < N; i++)
    d_[i] = 0;
}
template<int N>
ALWAYS_INLINE
SurrealS<N>::SurrealS( const Real v0, const Real d0[], int n )
{
  SANS_ASSERT( n == N );
  v_ = v0;
  for (int i = 0; i < N; i++)
    d_[i] = d0[i];
}
template<int N>
ALWAYS_INLINE
SurrealS<N>::SurrealS( const Real v0, const Real& d0 )
{
  v_ = v0;
  for (int i = 0; i < N; i++)
    d_[i] = d0;
}


namespace SurrealSExpr
{

// Lazy expressions


// Addition and Subtraction

template<class ExprL, class ExprR>
class OpAdd : public SurrealSType< OpAdd<ExprL,ExprR> >
{
public:
  static const int N = ExprL::N;
  BOOST_MPL_ASSERT_RELATION( ExprL::N, ==, ExprR::N );

  ALWAYS_INLINE
  OpAdd(const ExprL& eL, const ExprR& eR) : eL(eL), eR(eR) {}

  ALWAYS_INLINE Real value() const { return eL.value() + eR.value(); }
  ALWAYS_INLINE Real deriv(const int& i) const { return eL.deriv(i) + eR.deriv(i); }

  static const int nArgsL = ExprL::nArgs;
  static const int nArgsR = ExprR::nArgs;
  static const int nArgs = nArgsL + nArgsR;

  void getPartials(const Real& bar, Real partials[]) const
  {
    eL.getPartials(bar, partials);
    eR.getPartials(bar, partials+nArgsL);
  }
  template<int arg>
  const Real& getArgDeriv(const int i) const
  {
    if ( arg < nArgsL )
      return eL.template getArgDeriv<arg>(i);
    else
      return eR.template getArgDeriv<arg-nArgsL>(i);
  }

  ALWAYS_INLINE const OpAdd&
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return N; }
private:
  const ExprL& eL;
  const ExprR& eR;
};

}

template<class L, class R>
ALWAYS_INLINE SurrealSExpr::OpAdd<L,R>
operator+(const SurrealSType<L>& Ll, const SurrealSType<R>& Rr)
{
  return SurrealSExpr::OpAdd<L,R>( Ll.cast(), Rr.cast() );
}

namespace SurrealSExpr
{

template<class ExprL, class ExprR>
class OpSub : public SurrealSType< OpSub<ExprL,ExprR> >
{
public:
  static const int N = ExprL::N;
  BOOST_MPL_ASSERT_RELATION( ExprL::N, ==, ExprR::N );

  ALWAYS_INLINE
  OpSub(const ExprL& eL, const ExprR& eR) : eL(eL), eR(eR) {}

  ALWAYS_INLINE Real value() const { return eL.value() - eR.value(); }
  ALWAYS_INLINE Real deriv(const int& i) const { return eL.deriv(i) - eR.deriv(i); }

  static const int nArgsL = ExprL::nArgs;
  static const int nArgsR = ExprR::nArgs;
  static const int nArgs = nArgsL + nArgsR;

  void getPartials(const Real& bar, Real partials[]) const
  {
    eL.getPartials(bar, partials);
    eR.getPartials(-bar, partials+nArgsL);
  }
  template<int arg>
  const Real& getArgDeriv(const int i) const
  {
    if ( arg < nArgsL )
      return eL.template getArgDeriv<arg>(i);
    else
      return eR.template getArgDeriv<arg-nArgsL>(i);
  }

  ALWAYS_INLINE const OpSub&
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return eL.size(); }
private:
  const ExprL& eL;
  const ExprR& eR;
};

}

template<class L, class R>
ALWAYS_INLINE SurrealSExpr::OpSub<L, R>
operator-(const SurrealSType<L>& Ll, const SurrealSType<R>& Rr)
{
  return SurrealSExpr::OpSub<L, R>( Ll.cast(), Rr.cast() );
}

//Addition and Subtraction with scalar quantities

namespace SurrealSExpr
{

template<class Expr>
class OpScalar : public SurrealSType< OpScalar<Expr> >
{
public:
  static const int N = Expr::N;

  ALWAYS_INLINE
  OpScalar(const Expr& e, const Real esgn, const Real s) : e(e), esgn(esgn), s(s) {}

  ALWAYS_INLINE Real value() const { return esgn*e.value() + s; }
  ALWAYS_INLINE Real deriv(const int& i) const { return esgn*e.deriv(i); }

  static const int nArgs = Expr::nArgs;

  void getPartials(const Real& bar, Real partials[]) const
  {
    e.getPartials(bar*esgn, partials);
  }
  template<int arg>
  const Real& getArgDeriv(const int i) const
  {
    return e.template getArgDeriv<arg>(i);
  }

  ALWAYS_INLINE const OpScalar&
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return e.size(); }
private:
  const Expr& e;
  const Real esgn;
  const Real s;
};

}

template<class Expr>
ALWAYS_INLINE SurrealSExpr::OpScalar<Expr>
operator+(const SurrealSType<Expr>& e, const Real& s)
{
  return SurrealSExpr::OpScalar<Expr>( e.cast(), 1, s );
}

template<class Expr>
ALWAYS_INLINE SurrealSExpr::OpScalar<Expr>
operator+(const Real& s, const SurrealSType<Expr>& e)
{
  return SurrealSExpr::OpScalar<Expr>( e.cast(), 1, s );
}

template<class Expr>
ALWAYS_INLINE SurrealSExpr::OpScalar<Expr>
operator-(const SurrealSType<Expr>& e, const Real& s)
{
  return SurrealSExpr::OpScalar<Expr>( e.cast(), 1, -s );
}

template<class Expr>
ALWAYS_INLINE SurrealSExpr::OpScalar<Expr>
operator-(const Real& s, const SurrealSType<Expr>& e)
{
  return SurrealSExpr::OpScalar<Expr>( e.cast(), -1, s );
}


//Multiplication with SurrealSs

namespace SurrealSExpr
{

template<class ExprL, class ExprR>
class OpMul : public SurrealSType< OpMul<ExprL, ExprR> >
{
public:
  static const int N = ExprL::N;
  BOOST_MPL_ASSERT_RELATION( ExprL::N, ==, ExprR::N );

  ALWAYS_INLINE
  OpMul(const ExprL& eL, const ExprR& eR) : eL(eL), eR(eR), eL_val(eL.value()), eR_val(eR.value()) {}

  ALWAYS_INLINE Real value() const { return eL_val*eR_val; }
  ALWAYS_INLINE Real deriv(const int& i) const { return eL_val*eR.deriv(i) + eL.deriv(i)*eR_val; }

  static const int nArgsL = ExprL::nArgs;
  static const int nArgsR = ExprR::nArgs;
  static const int nArgs = nArgsL + nArgsR;

  void getPartials(const Real& bar, Real partials[]) const
  {
    eL.getPartials(bar*eR_val, partials);
    eR.getPartials(bar*eL_val, partials+nArgsL);
  }
  template<int arg>
  const Real& getArgDeriv(const int i) const
  {
    if ( arg < nArgsL )
      return eL.template getArgDeriv<arg>(i);
    else
      return eR.template getArgDeriv<arg-nArgsL>(i);
  }

  ALWAYS_INLINE const OpMul&
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return eL.size(); }
private:
  const ExprL& eL;
  const ExprR& eR;
  const Real eL_val, eR_val;
};

template<class ExprL>
class OpMul<ExprL, Real> : public SurrealSType< OpMul<ExprL, Real> >
{
public:
  static const int N = ExprL::N;

  ALWAYS_INLINE
  OpMul(const ExprL& e, const Real s) : e(e), s(s) {}

  ALWAYS_INLINE Real value() const { return e.value()*s; }
  ALWAYS_INLINE Real deriv(const int& i) const { return e.deriv(i)*s; }

  static const int nArgs = ExprL::nArgs;

  void getPartials(const Real& bar, Real partials[]) const
  {
    e.getPartials(bar*s, partials);
  }
  template<int arg>
  const Real& getArgDeriv(const int i) const
  {
    return e.template getArgDeriv<arg>(i);
  }

  ALWAYS_INLINE const OpMul&
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return e.size(); }

  const ExprL& e;
  const Real s;
};
}

template<class T>
struct is_arithmetic_not_SurrealS : boost::mpl::and_< boost::is_arithmetic<T>,
                                                     boost::mpl::not_< boost::is_base_and_derived<SurrealSTypeBase, T> > >
{};


//=============================================================================
template<class ExprL, class ExprR>
ALWAYS_INLINE SurrealSExpr::OpMul<ExprL, ExprR>
operator*(const SurrealSType<ExprL>& z1, const SurrealSType<ExprR>& z2)
{
  return SurrealSExpr::OpMul<ExprL, ExprR>( z1.cast(), z2.cast() );
}

template<class Expr, typename T>
ALWAYS_INLINE typename boost::enable_if< is_arithmetic_not_SurrealS<T>,
                                         SurrealSExpr::OpMul<Expr, Real> >::type
operator*(const SurrealSType<Expr>& e, const T& s)
{
  return SurrealSExpr::OpMul<Expr, Real>( e.cast(), s );
}

template<class Expr, typename T>
ALWAYS_INLINE typename boost::enable_if< is_arithmetic_not_SurrealS<T>,
                                         SurrealSExpr::OpMul<Expr, Real> >::type
operator/(const SurrealSType<Expr>& e, const T& s)
{
  return SurrealSExpr::OpMul<Expr, Real>( e.cast(), Real(1)/s );
}

template<class Expr, typename T>
ALWAYS_INLINE typename boost::enable_if< is_arithmetic_not_SurrealS<T>,
                                         SurrealSExpr::OpMul<Expr, Real> >::type
operator*(const T& s, const SurrealSType<Expr>& e)
{
  return SurrealSExpr::OpMul<Expr, Real>( e.cast(), s );
}

//=============================================================================
// This is a special case when multiplies scalars are multiplying from two sides, i.e. B = 2*A*2;
// This reduces the complexity of the expression tree and hence reduces code bloat
template<class Expr, typename T>
ALWAYS_INLINE typename boost::enable_if< is_arithmetic_not_SurrealS<T>,
                                         SurrealSExpr::OpMul<Expr, Real> >::type
operator*(const SurrealSExpr::OpMul<Expr, Real>& MulScal, const T& s)
{
  return SurrealSExpr::OpMul<Expr, Real>( MulScal.e, MulScal.s*s );
}

template<class Expr, typename T>
ALWAYS_INLINE typename boost::enable_if< is_arithmetic_not_SurrealS<T>,
                                         SurrealSExpr::OpMul<Expr, Real> >::type
operator/(const SurrealSExpr::OpMul<Expr, Real>& MulScal, const T& s)
{
  return SurrealSExpr::OpMul<Expr, Real>( MulScal.e, MulScal.s/s );
}

template<class Expr, typename T>
ALWAYS_INLINE typename boost::enable_if< is_arithmetic_not_SurrealS<T>,
                                         SurrealSExpr::OpMul<Expr, Real> >::type
operator*(const T& s, const SurrealSExpr::OpMul<Expr, Real>& MulScal)
{
  return SurrealSExpr::OpMul<Expr, Real>( MulScal.e, MulScal.s*s );
}


//=============================================================================
//Division with SurrealSs

namespace SurrealSExpr
{

template<class ExprL, class ExprR>
class OpDiv : public SurrealSType< OpDiv<ExprL, ExprR> >
{
public:
  static const int N = ExprL::N;
  BOOST_MPL_ASSERT_RELATION( ExprL::N, ==, ExprR::N );

  ALWAYS_INLINE
  OpDiv(const ExprL& eL, const ExprR& eR) : eL(eL), eR(eR), eL_val(eL.value()), eR_val(eR.value())
                                          , vali(1/(eR_val*eR_val)) {}

  ALWAYS_INLINE Real value() const { return eL_val/eR_val; }
  ALWAYS_INLINE Real deriv(const int& i) const { return (eR_val*eL.deriv(i) - eR.deriv(i)*eL_val)*vali; }

  static const int nArgsL = ExprL::nArgs;
  static const int nArgsR = ExprR::nArgs;
  static const int nArgs = nArgsL + nArgsR;

  void getPartials(const Real& bar, Real partials[]) const
  {
    eL.getPartials( bar*eR_val*vali, partials);
    eR.getPartials(-bar*eL_val*vali, partials+nArgsL);
  }
  template<int arg>
  const Real& getArgDeriv(const int i) const
  {
    if ( arg < nArgsL )
      return eL.template getArgDeriv<arg>(i);
    else
      return eR.template getArgDeriv<arg-nArgsL>(i);
  }

  ALWAYS_INLINE const OpDiv&
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return eL.size(); }
private:
  const ExprL& eL;
  const ExprR& eR;
  const Real eL_val, eR_val, vali;
};

}

template< class ExprL, class ExprR >
ALWAYS_INLINE SurrealSExpr::OpDiv<ExprL, ExprR>
operator/(const SurrealSType<ExprL>& eL, const SurrealSType<ExprR>& eR)
{
  return SurrealSExpr::OpDiv<ExprL, ExprR>( eL.cast(), eR.cast() );
}


namespace SurrealSExpr
{

template<class Expr>
class OpDivScalarNumerator : public SurrealSType< OpDivScalarNumerator<Expr> >
{
public:
  static const int N = Expr::N;

  ALWAYS_INLINE
  OpDivScalarNumerator(const Expr& e, const Real& s) : e(e), s(s), e_val(e.value()), se_val2i(s/(e_val*e_val)) {}

  ALWAYS_INLINE Real value() const { return s/e_val; }
  ALWAYS_INLINE Real deriv(const int& i) const { return -se_val2i*e.deriv(i); }

  static const int nArgs = Expr::nArgs;

  void getPartials(const Real& bar, Real partials[]) const
  {
    e.getPartials( -bar*se_val2i, partials);
  }
  template<int arg>
  const Real& getArgDeriv(const int i) const
  {
    return e.template getArgDeriv<arg>(i);
  }

  ALWAYS_INLINE const OpDivScalarNumerator&
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return e.size(); }
private:
  const Expr& e;
  const Real s;
  const Real e_val, se_val2i;
};

}

template<class Expr, typename T>
ALWAYS_INLINE typename boost::enable_if< is_arithmetic_not_SurrealS<T>,
                                         SurrealSExpr::OpDivScalarNumerator<Expr> >::type
operator/(const T& s, const SurrealSType<Expr>& e)
{
  return SurrealSExpr::OpDivScalarNumerator<Expr>( e.cast(), s );
}


// assignment

template<int N>
ALWAYS_INLINE SurrealS<N>&
SurrealS<N>::operator=( const SurrealS& z )
{
  //Do nothing if assigning self to self
  if ( &z == this ) return *this;

  v_ = z.v_;
  for (int i = 0; i < N; i++)
    d_[i] = z.d_[i];
  return *this;
}

template<int N>
ALWAYS_INLINE SurrealS<N>&
SurrealS<N>::operator=( const Real& r )
{
  v_ = r;
  for (int i = 0; i < N; i++)
    d_[i] = 0;
  return *this;
}

template<int N>
template< class Expr >
ALWAYS_INLINE SurrealS<N>&
SurrealS<N>::operator=( const SurrealSType<Expr>& r )
{
  const Expr& Tree = r.cast();

  BOOST_MPL_ASSERT_RELATION( N, ==, Expr::N );

  Real sum;
  AccumOp<Expr> Accum(Tree, sum);

  for (int i = 0; i < N; ++i)
  {
    Accum.getDeriv(i);
    boost::mpl::for_each< boost::mpl::range_c< int, 0, Expr::nArgs > >(Accum);
    d_[i] = sum;
  }

  //Value must be set last as it might be used in the derivative calculation
  v_ = Tree.value();

  return *this;
}

template<int N>
template< class Expr >
ALWAYS_INLINE SurrealS<N>&
SurrealS<N>::operator+=( const SurrealSType<Expr>& r )
{
  const Expr& Tree = r.cast();

  BOOST_MPL_ASSERT_RELATION( N, ==, Expr::N );

  Real sum;
  AccumOp<Expr> Accum(Tree, sum);

  for (int i = 0; i < N; ++i)
  {
    Accum.getDeriv(i);
    boost::mpl::for_each< boost::mpl::range_c< int, 0, Expr::nArgs > >(Accum);
    d_[i] += sum;
  }

  //Value must be set last as it might be used in the derivative calculation
  v_ += Tree.value();

  return *this;
}

template<int N>
template< class Expr >
ALWAYS_INLINE SurrealS<N>&
SurrealS<N>::operator-=( const SurrealSType<Expr>& r )
{
  const Expr& Tree = r.cast();

  BOOST_MPL_ASSERT_RELATION( N, ==, Expr::N );

  Real sum;
  AccumOp<Expr> Accum(Tree, sum);

  for (int i = 0; i < N; ++i)
  {
    Accum.getDeriv(i);
    boost::mpl::for_each< boost::mpl::range_c< int, 0, Expr::nArgs > >(Accum);
    d_[i] -= sum;
  }

  //Value must be set last as it might be used in the derivative calculation
  v_ -= Tree.value();

  return *this;
}


// unary operators; no side effects

template<int N>
ALWAYS_INLINE const SurrealS<N>&
SurrealS<N>::operator+() const
{
  return *this;
}

template< class Expr >
ALWAYS_INLINE const SurrealSExpr::OpMul<Expr, Real>
operator-(SurrealSType<Expr> const& e)
{
  return SurrealSExpr::OpMul<Expr, Real>( e.cast(), -1 );
}

template< class Expr >
ALWAYS_INLINE const SurrealSExpr::OpMul<Expr, Real>
operator-(SurrealSExpr::OpMul<Expr, Real> const& Mul)
{
  return SurrealSExpr::OpMul<Expr, Real>( Mul.e, -1*Mul.s );
}

// binary accumulation operators


template<int N>
ALWAYS_INLINE SurrealS<N>&
SurrealS<N>::operator+=( const Real& r )
{
  v_ += r;
  return *this;
}

template<int N>
ALWAYS_INLINE SurrealS<N>&
SurrealS<N>::operator-=( const Real& r )
{
  v_ -= r;
  return *this;
}

template<int N>
ALWAYS_INLINE SurrealS<N>&
SurrealS<N>::operator*=( const SurrealS& z )
{
  for (int i = 0; i < N; i++)
    d_[i] = v_*z.d_[i] + d_[i]*z.v_;
  v_ *= z.v_;
  return *this;
}

template<int N>
ALWAYS_INLINE SurrealS<N>&
SurrealS<N>::operator*=( const Real& r )
{
  for (int i = 0; i < N; i++)
    d_[i] *= r;
  v_ *= r;
  return *this;
}

template<int N>
ALWAYS_INLINE SurrealS<N>&
SurrealS<N>::operator/=( const SurrealS& z)
{
  Real tmp = 1./(z.v_*z.v_);
  for (int i = 0; i < N; i++)
    d_[i] = (z.v_*d_[i] - v_*z.d_[i])*tmp;
  v_ /= z.v_;
  return *this;
}

template<int N>
ALWAYS_INLINE SurrealS<N>&
SurrealS<N>::operator/=( const Real& r )
{
  Real tmp = 1./r;
  for (int i = 0; i < N; i++)
    d_[i] *= tmp;
  v_ *= tmp;
  return *this;
}

// relational operators

template<class ExprL, class ExprR>
ALWAYS_INLINE bool
operator==( const SurrealSType<ExprL>& lhs, const SurrealSType<ExprR>& rhs )
{
  return lhs.value() == rhs.value();
}

template<class Expr>
ALWAYS_INLINE bool
operator==( const SurrealSType<Expr>& lhs, const Real& rhs )
{
  return lhs.value() == rhs;
}

template<class Expr>
ALWAYS_INLINE bool
operator==( const Real& lhs, const SurrealSType<Expr>& rhs )
{
  return lhs == rhs.value();
}

template<class ExprL, class ExprR>
ALWAYS_INLINE bool
operator!=( const SurrealSType<ExprL>& lhs, const SurrealSType<ExprR>& rhs )
{
  return lhs.value() != rhs.value();
}

template<class Expr>
ALWAYS_INLINE bool
operator!=( const SurrealSType<Expr>& lhs, const Real& rhs )
{
  return lhs.value() != rhs;
}

template<class Expr>
ALWAYS_INLINE bool
operator!=( const Real& lhs, const SurrealSType<Expr>& rhs )
{
  return lhs != rhs.value();
}

template<class ExprL, class ExprR>
ALWAYS_INLINE bool
operator>( const SurrealSType<ExprL>& lhs, const SurrealSType<ExprR>& rhs )
{
  return lhs.value() > rhs.value();
}

template<class Expr>
ALWAYS_INLINE bool
operator>( const SurrealSType<Expr>& lhs, const Real& rhs )
{
  return lhs.value() > rhs;
}

template<class Expr>
ALWAYS_INLINE bool
operator>( const Real& lhs, const SurrealSType<Expr>& rhs )
{
  return lhs > rhs.value();
}

template<class ExprL, class ExprR>
ALWAYS_INLINE bool
operator<( const SurrealSType<ExprL>& lhs, const SurrealSType<ExprR>& rhs )
{
  return lhs.value() < rhs.value();
}

template<class Expr>
ALWAYS_INLINE bool
operator<( const SurrealSType<Expr>& lhs, const Real& rhs )
{
  return lhs.value() < rhs;
}

template<class Expr>
ALWAYS_INLINE bool
operator<( const Real& lhs, const SurrealSType<Expr>& rhs )
{
  return lhs < rhs.value();
}

template<class ExprL, class ExprR>
ALWAYS_INLINE bool
operator>=( const SurrealSType<ExprL>& lhs, const SurrealSType<ExprR>& rhs )
{
  return lhs.value() >= rhs.value();
}

template<class Expr>
ALWAYS_INLINE bool
operator>=( const SurrealSType<Expr>& lhs, const Real& rhs )
{
  return lhs.value() >= rhs;
}

template<class Expr>
ALWAYS_INLINE bool
operator>=( const Real& lhs, const SurrealSType<Expr>& rhs )
{
  return lhs >= rhs.value();
}

template<class ExprL, class ExprR>
ALWAYS_INLINE bool
operator<=( const SurrealSType<ExprL>& lhs, const SurrealSType<ExprR>& rhs )
{
  return lhs.value() <= rhs.value();
}

template<class Expr>
ALWAYS_INLINE bool
operator<=( const SurrealSType<Expr>& lhs, const Real& rhs )
{
  return lhs.value() <= rhs;
}

template<class Expr>
ALWAYS_INLINE bool
operator<=( const Real& lhs, const SurrealSType<Expr>& rhs )
{
  return lhs <= rhs.value();
}


//Functions for SurrealSs
#define SURREALS_FUNC1( NAME, FUNC, DERIV ) \
namespace SurrealSExpr \
{  \
template<class Expr> \
class BOOST_PP_CAT(SurrealS_, NAME) : public SurrealSType< BOOST_PP_CAT(SurrealS_, NAME)<Expr> > \
{ /*This is for functions when the argument is an expression*/ \
public: \
  static const int N = Expr::N; \
  \
  ALWAYS_INLINE \
  BOOST_PP_CAT(SurrealS_, NAME)(const Expr& e) : e(e), z(e.value()), der(DERIV) {} \
  \
  ALWAYS_INLINE Real value() const { return FUNC; } \
  ALWAYS_INLINE Real deriv(const int& i) const { return der*e.deriv(i); } \
  \
  static const int nArgs = Expr::nArgs; \
  \
  void getPartials(const Real& bar, Real partials[]) const \
  { \
    e.getPartials( bar*der, partials); \
  } \
  template<int arg> \
  const Real& getArgDeriv(const int i) const \
  { \
    return e.template getArgDeriv<arg>(i); \
  } \
  \
  ALWAYS_INLINE const BOOST_PP_CAT(SurrealS_, NAME)& \
  operator+() const { return *this; } \
  ALWAYS_INLINE int size() const { return e.size(); } \
private: \
  const Expr& e; \
  const Real z, der; \
}; \
} \
\
template<class Expr> \
ALWAYS_INLINE SurrealSExpr::BOOST_PP_CAT(SurrealS_, NAME)<Expr> \
NAME(const SurrealSType<Expr>& z) { return SurrealSExpr::BOOST_PP_CAT(SurrealS_, NAME)<Expr>( z.cast() ); }


#define SURREALS_FUNC2( NAME, FUNC, DERIV ) \
namespace SurrealSExpr \
{  \
template<class ExprL, class ExprR> \
class BOOST_PP_CAT(SurrealS_, NAME) : public SurrealSType< BOOST_PP_CAT(SurrealS_, NAME)<ExprL, ExprR> > \
{ /*This is for functions when the argument is an expression*/ \
public: \
  static const int N = ExprL::N; \
  BOOST_MPL_ASSERT_RELATION(ExprL::N, ==, ExprR::N ); \
  \
  ALWAYS_INLINE \
  BOOST_PP_CAT(SurrealS_, NAME)(const ExprL& eL, const ExprR& eR) : eL(eL), eR(eR), z1(eL.value()), z2(eR.value()), \
                                                                 der(DERIV) {} \
  \
  ALWAYS_INLINE Real value() const { return FUNC; } \
  ALWAYS_INLINE Real deriv(const int& i) const { return der*(z2*eL.deriv(i) - z1*eR.deriv(i)); } \
  \
  static const int nArgsL = ExprL::nArgs; \
  static const int nArgsR = ExprR::nArgs; \
  static const int nArgs = nArgsL + nArgsR; \
  \
  void getPartials(const Real& bar, Real partials[]) const \
  { \
    eL.getPartials( bar*z2*der, partials); \
    eR.getPartials(-bar*z1*der, partials+nArgsL); \
  } \
  template<int arg> \
  const Real& getArgDeriv(const int i) const \
  { \
    if ( arg < nArgsL ) \
      return eL.template getArgDeriv<arg>(i); \
    else \
      return eR.template getArgDeriv<arg-nArgsL>(i); \
  } \
  \
  ALWAYS_INLINE const BOOST_PP_CAT(SurrealS_, NAME)& \
  operator+() const { return *this; } \
  ALWAYS_INLINE int size() const { return eL.size(); } \
private: \
  const ExprL& eL; \
  const ExprR& eR; \
  const Real z1, z2, der; \
}; \
  \
} \
\
template<class ExprL, class ExprR> \
ALWAYS_INLINE SurrealSExpr::BOOST_PP_CAT(SurrealS_, NAME)<ExprL, ExprR> \
NAME(const SurrealSType<ExprL>& z1, const SurrealSType<ExprR>& z2) \
{ return SurrealSExpr::BOOST_PP_CAT(SurrealS_, NAME)<ExprL, ExprR>( z1.cast(), z2.cast() ); }

// trig functions <cmath>

SURREALS_FUNC1( cos, std::cos(z), -std::sin(z) )
SURREALS_FUNC1( sin, std::sin(z),  std::cos(z) )
SURREALS_FUNC1( tan, std::tan(z),  Real(1)/(std::cos(z)*std::cos(z)) )
SURREALS_FUNC1( acos, std::acos(z), -Real(1)/std::sqrt(1 - z*z) )
SURREALS_FUNC1( asin, std::asin(z),  Real(1)/std::sqrt(1 - z*z) )
SURREALS_FUNC1( atan, std::atan(z),  Real(1)/(1 + z*z) )

SURREALS_FUNC2( atan2, std::atan2(z1, z2),  Real(1)/(z1*z1 + z2*z2) )

// hyperbolic functions <cmath>

SURREALS_FUNC1( cosh, std::cosh(z), std::sinh(z) )
SURREALS_FUNC1( sinh, std::sinh(z), std::cosh(z) )
SURREALS_FUNC1( tanh, std::tanh(z), Real(1)/(std::cosh(z)*std::cosh(z)) )

// exp and log functions <cmath>

SURREALS_FUNC1( exp, std::exp(z), std::exp(z) )
SURREALS_FUNC1( log, std::log(z), Real(1)/z )
SURREALS_FUNC1( log10, std::log10(z), Real(1)/(z*std::log(10.)) )

// power functions <cmath>

namespace SurrealSExpr
{

template<class ExprL, class ExprR>
class SurrealS_pow : public SurrealSType< SurrealS_pow<ExprL, ExprR> >
{ /*This is for functions when the argument is an expression*/
public:
  static const int N = ExprL::N;
  BOOST_MPL_ASSERT_RELATION(ExprL::N, ==, ExprR::N);

  ALWAYS_INLINE
  SurrealS_pow(const ExprL& eL, const ExprR& eR) : eL(eL), eR(eR), a(eL.value()), b(eR.value()),
                                                   powab(std::pow(a,b)),
                                                   tmp1( (a == 0) ? ((b == 1) ? 1 : 0) : b*std::pow(a, b - 1) ),
                                                   tmp2( (a == 0) ? 0 : powab*std::log(a) ) {}

  ALWAYS_INLINE Real value() const { return powab; }
  ALWAYS_INLINE Real deriv(const int& i) const { return tmp1*eL.deriv(i) + tmp2*eR.deriv(i); }

  static const int nArgsL = ExprL::nArgs;
  static const int nArgsR = ExprR::nArgs;
  static const int nArgs = nArgsL + nArgsR;

  void getPartials(const Real& bar, Real partials[]) const
  {
    eL.getPartials(bar*tmp1, partials);
    eR.getPartials(bar*tmp2, partials+nArgsL);
  }
  template<int arg>
  const Real& getArgDeriv(const int i) const
  {
    if ( arg < nArgsL )
      return eL.template getArgDeriv<arg>(i);
    else
      return eR.template getArgDeriv<arg-nArgsL>(i);
  }

  ALWAYS_INLINE const SurrealS_pow&
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return eL.size(); }
private:
  const ExprL& eL;
  const ExprR& eR;
  const Real a, b, powab, tmp1, tmp2;
};

template<class ExprL>
class SurrealS_pow<ExprL, Real> : public SurrealSType< SurrealS_pow<ExprL, Real> >
{ /*This is optimized when the argument is SurrealS and Real*/
public:
  static const int N = ExprL::N;

  ALWAYS_INLINE
  SurrealS_pow(const ExprL& eL, const Real& b) : eL(eL), a(eL.value()),
                                                 powab(std::pow(a,b)),
                                                 tmp1( (a == 0) ? ((b == 1) ? 1 : 0) : b*std::pow(a, b - 1) ) {}

  ALWAYS_INLINE Real value() const { return powab; }
  ALWAYS_INLINE Real deriv(const int& i) const { return tmp1*eL.deriv(i); }

  static const int nArgs = ExprL::nArgs;

  void getPartials(const Real& bar, Real partials[]) const
  {
    eL.getPartials(bar*tmp1, partials);
  }
  template<int arg>
  const Real& getArgDeriv(const int i) const
  {
    return eL.template getArgDeriv<arg>(i);
  }

  ALWAYS_INLINE const SurrealS_pow
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return eL.size(); }
private:
  const ExprL& eL;
  const Real a, powab, tmp1;
};


template<class ExprR>
class SurrealS_pow<Real, ExprR> : public SurrealSType< SurrealS_pow<Real, ExprR> >
{ /*This is optimized when the argument is a Real and SurrealS*/
public:
  static const int N = ExprR::N;

  ALWAYS_INLINE
  SurrealS_pow(const Real& a, const ExprR& eR) : eR(eR), b(eR.value()),
                                               powab(std::pow(a,b)),
                                               tmp2( (a == 0) ? 0 : powab*std::log(a) ) {}

  ALWAYS_INLINE Real value() const { return powab; }
  ALWAYS_INLINE Real deriv(const int& i) const { return tmp2*eR.deriv(i); }

  static const int nArgs = ExprR::nArgs;

  void getPartials(const Real& bar, Real partials[]) const
  {
    eR.getPartials(bar*tmp2, partials);
  }
  template<int arg>
  const Real& getArgDeriv(const int i) const
  {
    return eR.template getArgDeriv<arg>(i);
  }

  ALWAYS_INLINE const SurrealS_pow
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return eR.size(); }
private:
  const ExprR& eR;
  const Real b, powab, tmp2;

};

}

template<class ExprL, class ExprR>
ALWAYS_INLINE SurrealSExpr::SurrealS_pow<ExprL, ExprR>
pow(const SurrealSType<ExprL>& a, const SurrealSType<ExprR>& b)
{
  return SurrealSExpr::SurrealS_pow<ExprL, ExprR>( a.cast(), b.cast() );
}

template<class Expr, typename T>
ALWAYS_INLINE typename boost::enable_if< is_arithmetic_not_SurrealS<T>,
                                  SurrealSExpr::SurrealS_pow<Expr, Real> >::type
pow(const SurrealSType<Expr>& a, const T& b )
{
  return SurrealSExpr::SurrealS_pow<Expr, Real>( a.cast(), b );
}

template<class Expr, typename T>
ALWAYS_INLINE typename boost::enable_if< is_arithmetic_not_SurrealS<T>,
                                  SurrealSExpr::SurrealS_pow<Real, Expr> >::type
pow(const T& a, const SurrealSType<Expr>& b)
{
  return SurrealSExpr::SurrealS_pow<Real, Expr>( a, b.cast() );
}


namespace SurrealSExpr
{

template<class Expr>
class SurrealS_sqrt : public SurrealSType< SurrealS_sqrt<Expr> >
{ /*This is optimized when the argument is an Expression*/
public:
  static const int N = Expr::N;

  ALWAYS_INLINE
  SurrealS_sqrt(const Expr& e) : e(e), sqrtv( sqrt(e.value()) ), tmp( sqrtv == 0 ? 0 : 0.5/sqrtv ) {}

  ALWAYS_INLINE Real value() const { return sqrtv; }
  ALWAYS_INLINE Real deriv(const int& i) const { return tmp*e.deriv(i); }

  static const int nArgs = Expr::nArgs;

  void getPartials(const Real& bar, Real partials[]) const
  {
    e.getPartials(bar*tmp, partials);
  }
  template<int arg>
  const Real& getArgDeriv(const int i) const
  {
    return e.template getArgDeriv<arg>(i);
  }

  ALWAYS_INLINE const SurrealS_sqrt
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return e.size(); }
private:
  const Expr& e;
  const Real sqrtv, tmp;
};

}

template<class Expr>
ALWAYS_INLINE SurrealSExpr::SurrealS_sqrt<Expr>
sqrt(const SurrealSType<Expr>& z)
{
  return SurrealSExpr::SurrealS_sqrt<Expr>( z.cast() );
}


// rounding functions <cmath>

SURREALS_FUNC1( ceil, std::ceil(z), 0 )
SURREALS_FUNC1( floor, std::floor(z), 0 )

// misc functions <cmath>

template<int N>
ALWAYS_INLINE SurrealSExpr::OpMul<SurrealS<N>, Real>
abs( const SurrealS<N>& z )
{
  return (z.value() < 0) ?
         SurrealSExpr::OpMul<SurrealS<N>, Real>( z, -1 ) :
         SurrealSExpr::OpMul<SurrealS<N>, Real>( z,  1 );
}

template<int N>
ALWAYS_INLINE SurrealSExpr::OpMul<SurrealS<N>, Real>
fabs( const SurrealS<N>& z )
{
  return (z.value() < 0) ?
         SurrealSExpr::OpMul<SurrealS<N>, Real>( z, -1 ) :
         SurrealSExpr::OpMul<SurrealS<N>, Real>( z,  1 );
}


// I/O

template<int N>
std::istream&
operator>>( std::istream& is, SurrealS<N>& z )
{
  Real v = 0;
  Real d[10] = {0};
  char c = 0;
  int n = 0;

  is >> c;
  if (c == '(')
  {
    is >> v;

    is >> c;
    bool done = false;
    while (! done)
    {
      if (c != ')') is.clear(std::ios::badbit);
      if (c == ',')
      {
        is >> d[n]; n++;
      }
      else if (c == ')')
      {
        done = true;
      }
    }
  }
  else
  {
    is.putback(c);
    is >> v;
  }

  if (is) z = SurrealS<N>(v, d, n);
  return is;
}


template<class Expr>
std::ostream&
operator<<( std::ostream& os, const SurrealSType<Expr>& ztype )
{
  const Expr& z = ztype.cast();
  os << '(' << z.value() << ';';
  for (int i = 0; i < Expr::N - 1; i++)
    os << z.deriv(i) << ',';
  os << z.deriv(Expr::N - 1) << ')';
  return os;
}


//Created specialized version of fpt_abs to work with the boost unit testing framework

namespace boost
{

namespace test_tools
{

namespace tt_detail
{

template<typename FPT>
FPT
fpt_abs( FPT fpv );

template<typename FPT>
FPT
safe_fpt_division( FPT f1, FPT f2 );



template<int N>
inline Real
fpt_abs( const SurrealSExpr::OpMul<SurrealS<N>, Real>& fpv )
{
  Real val = fpv.value();
  return fpt_abs( val );
}

template<int N>
inline Real
fpt_abs( const SurrealSExpr::OpSub<SurrealS<N>, SurrealS<N> >& fpv )
{
  Real val = fpv.value();
  return fpt_abs( val );
}

// both f1 and f2 are unsigned here
template<int N>
inline Real
safe_fpt_division( const SurrealS<N>& f1, const SurrealS<N>& f2 )
{
  Real val1 = f1.value();
  Real val2 = f2.value();
  return safe_fpt_division( val1, val2 );
}

}
}
}


//Clean up macro definitions
#undef SURREALS_FUNC1
#undef SURREALS_FUNC2

#endif // SURREALS_REVERSE_H
