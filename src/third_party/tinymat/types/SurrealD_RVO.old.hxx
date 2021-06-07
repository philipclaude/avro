// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef SURREALD_RVO_H
#define SURREALD_RVO_H

//  overloaded derivative operator
//  ref: derivify.h (Google it)

#include <cmath>
#include <iostream>
#include <string>

#include "tools/SANSnumerics.h"     // Real
#include "tools/SANSException.h"
#include "tools/always_inline.h"

//----------------------------------------------------------------------------//
// SurrealD:  value, N derivatives
//
// Operators with Return Value Optimization (RVO)
//
// member functions:
//   .size          total number of derivatives
//   .value         value accessor
//   .deriv         nth derivative accessor (zero-based indexing)
//----------------------------------------------------------------------------//

class SurrealD
{
public:
  //The default constructor is intentionally not included here. This means Surral is
  //not initialized when declared, which is consistent with regular numbers. This also
  //improves performance.
  SurrealD( const SurrealD& z );
  explicit SurrealD( int n=0 );
  SurrealD( const Real& v0, const Real d0[], const int n );
  SurrealD( const Real& v0, const Real& d0, const int n );
  SurrealD( const int& v0, const int& d0, const int n );
  // cppcheck-suppress noExplicitConstructor
  SurrealD( const Real& v0 );
  ~SurrealD();

#if 0     // these are not needed and will be removed
  operator Real() const { return v_; }
  operator int() const { return int(v_); }
#endif

  int size() const { return N_; }

  // value accessor operators
  ALWAYS_INLINE       Real& value()       { return v_; }
  ALWAYS_INLINE const Real& value() const { return v_; }

  // derivative accessor operators
  ALWAYS_INLINE Real& deriv( int i=0 )       { SANS_ASSERT(N_ > 0); return d_[i]; }
  ALWAYS_INLINE Real  deriv( int i=0 ) const { return N_ > 0 ? d_[i] : 0; }

  // assignment
  SurrealD& operator=( const SurrealD& );
  SurrealD& operator=( const Real& );

  // unary operators; no side effects
  const SurrealD& operator+() const;
  const SurrealD  operator-() const;

  // binary accumulation operators
  SurrealD& operator+=( const SurrealD& );
  SurrealD& operator+=( const Real& );
  SurrealD& operator-=( const SurrealD& );
  SurrealD& operator-=( const Real& );
  SurrealD& operator*=( const SurrealD& );
  SurrealD& operator*=( const Real& );
  SurrealD& operator/=( const SurrealD& );
  SurrealD& operator/=( const Real& );

  // binary operators
  friend SurrealD operator+( const SurrealD&, const SurrealD& );
  friend SurrealD operator+( const SurrealD&, const Real& );
  friend SurrealD operator+( const Real&, const SurrealD& );
  friend SurrealD operator-( const SurrealD&, const SurrealD& );
  friend SurrealD operator-( const SurrealD&, const Real& );
  friend SurrealD operator-( const Real&, const SurrealD& );
  friend SurrealD operator*( const SurrealD&, const SurrealD& );
  friend SurrealD operator*( const SurrealD&, const Real& );
  friend SurrealD operator*( const Real&, const SurrealD& );
  friend SurrealD operator/( const SurrealD&, const SurrealD& );
  friend SurrealD operator/( const SurrealD&, const Real& );
  friend SurrealD operator/( const Real&, const SurrealD& );

  // relational operators
  friend bool operator==( const SurrealD&, const SurrealD& );
  friend bool operator==( const SurrealD&, const Real& );
  friend bool operator==( const Real&, const SurrealD& );
  friend bool operator!=( const SurrealD&, const SurrealD& );
  friend bool operator!=( const SurrealD&, const Real& );
  friend bool operator!=( const Real&, const SurrealD& );
  friend bool operator>( const SurrealD&, const SurrealD& );
  friend bool operator>( const SurrealD&, const Real& );
  friend bool operator>( const Real&, const SurrealD& );
  friend bool operator<( const SurrealD&, const SurrealD& );
  friend bool operator<( const SurrealD&, const Real& );
  friend bool operator<( const Real&, const SurrealD& );
  friend bool operator>=( const SurrealD&, const SurrealD& );
  friend bool operator>=( const SurrealD&, const Real& );
  friend bool operator>=( const Real&, const SurrealD& );
  friend bool operator<=( const SurrealD&, const SurrealD& );
  friend bool operator<=( const SurrealD&, const Real& );
  friend bool operator<=( const Real&, const SurrealD& );

  // trig functions <cmath>
  friend SurrealD cos( const SurrealD& );
  friend SurrealD sin( const SurrealD& );
  friend SurrealD tan( const SurrealD& );
  friend SurrealD acos( const SurrealD& );
  friend SurrealD asin( const SurrealD& );
  friend SurrealD atan( const SurrealD& );
  friend SurrealD atan2( const SurrealD&, const SurrealD& );

  // hyperbolic functions <cmath>
  friend SurrealD cosh( const SurrealD& );
  friend SurrealD sinh( const SurrealD& );
  friend SurrealD tanh( const SurrealD& );

  // exp and log functions <cmath>
  friend SurrealD exp( const SurrealD& );
  friend SurrealD log( const SurrealD& );
  friend SurrealD log10( const SurrealD& );

  // power functions <cmath>
  friend SurrealD pow( const SurrealD&, const SurrealD& );
  friend SurrealD pow( const SurrealD&, const Real& );
  friend SurrealD pow( const SurrealD&, const int& );
  friend SurrealD pow( const Real&, const SurrealD& );
  friend SurrealD pow( const int&, const SurrealD& );
  friend SurrealD sqrt( const SurrealD& );

  // rounding functions <cmath>
  friend SurrealD ceil( const SurrealD& );
  friend SurrealD floor( const SurrealD& );

  // misc functions <cmath>
  friend SurrealD abs( const SurrealD& );
  friend SurrealD fabs( const SurrealD& );

  // classification functions <cmath>
  friend bool isfinite( const SurrealD& );
  friend bool isinf( const SurrealD& );
  friend bool isnan( const SurrealD& );

  // input/output
  friend std::istream& operator>>( std::istream&, SurrealD& );
  friend std::ostream& operator<<( std::ostream&, const SurrealD& );

  void dump( int indentSize=0, std::ostream& out = std::cout ) const;

private:
  Real v_;          // value
  Real *d_;         // derivative array
  int N_;           // size of derivative array
};


// constructors

ALWAYS_INLINE
SurrealD::SurrealD( const SurrealD& z ) : v_(z.v_), d_(NULL), N_(z.N_)
{
  if (N_ > 0)
  {
    d_ = new Real[N_];
    for (int i = 0; i < N_; i++)
      d_[i] = z.d_[i];
  }
}

ALWAYS_INLINE
SurrealD::SurrealD( int n ) : v_(n), d_(NULL), N_(0) {}

ALWAYS_INLINE
SurrealD::SurrealD( const Real& v0, const Real d0[], const int n ) : v_(v0), d_(NULL), N_(n)
{
  SANS_ASSERT( N_ > 0 );

  d_ = new Real[N_];
  for (int i = 0; i < N_; i++)
    d_[i] = d0[i];
}

ALWAYS_INLINE
SurrealD::SurrealD( const Real& v0, const Real& d0, const int n ) : v_(v0), d_(NULL), N_(n)
{
  SANS_ASSERT( N_ >= 0 );

  if (N_ > 0)
  {
    d_ = new Real[N_];
    for (int i = 0; i < N_; i++)
      d_[i] = d0;
  }
}

ALWAYS_INLINE
SurrealD::SurrealD( const int& v0, const int& d0, const int n ) : v_(v0), d_(NULL), N_(n)
{
  SANS_ASSERT( N_ >= 0 );

  if (N_ > 0)
  {
    d_ = new Real[N_];
    for (int i = 0; i < N_; i++)
      d_[i] = d0;
  }
}

ALWAYS_INLINE
SurrealD::SurrealD( const Real& v0 ) : v_(v0), d_(NULL), N_(0) {}

ALWAYS_INLINE
SurrealD::~SurrealD()
{
  delete [] d_;
}


// assignment

ALWAYS_INLINE  SurrealD&
SurrealD::operator=( const SurrealD& z )
{
  //Do nothing if assigning self to self
  if ( &z == this ) return *this;

  if ((N_ == 0) && (z.N_ > 0))
  {
    N_ = z.N_;
    d_ = new Real[N_];
  }
  else if ( z.N_ == 0 )
  {
    //z has no derivatives, so it's a scalar
    (*this) = z.v_;
    return *this;
  }

  v_ = z.v_;
  for (int i = 0; i < N_; i++)
    d_[i] = z.d_[i];
  return *this;
}

ALWAYS_INLINE  SurrealD&
SurrealD::operator=( const Real& r )
{
  v_ = r;
  delete [] d_;
  d_ = NULL;
  N_ = 0;
//  if (N_ > 0)
//  {
//    for (int i = 0; i < N_; i++)
//      d_[i] = 0;
//  }
  return *this;
}


// unary operators; no side effects

ALWAYS_INLINE  const SurrealD&
SurrealD::operator+() const
{
  return *this;
}

ALWAYS_INLINE  const SurrealD
SurrealD::operator-() const
{
  if (N_ == 0)
  {
    return SurrealD( -v_ );
  }
  else
  {
    Real d0[N_];
    for (int i = 0; i < N_; i++)
      d0[i] = -d_[i];
    return SurrealD(-v_, d0, N_);
  }
}


// binary accumulation operators

ALWAYS_INLINE  SurrealD&
SurrealD::operator+=( const SurrealD& z )
{
  if ((N_ == 0) && (z.N_ > 0))
  {
    N_ = z.N_;
    d_ = new Real[N_]();    // NOTE: d_ is value-initialized here
  }
  else if ( z.N_ == 0 )
  {
    //z has no derivatives, so it is a scalar
    (*this) += z.v_;
    return *this;
  }

  v_ += z.v_;
  for (int i = 0; i < N_; i++)
    d_[i] += z.d_[i];

  return *this;
}

ALWAYS_INLINE  SurrealD&
SurrealD::operator+=( const Real& r )
{
  v_ += r;
  return *this;
}

ALWAYS_INLINE  SurrealD&
SurrealD::operator-=( const SurrealD& z )
{
  if ((N_ == 0) && (z.N_ > 0))
  {
    N_ = z.N_;
    d_ = new Real[N_]();    // NOTE: d_ is value-initialized here
  }
  else if ( z.N_ == 0 )
  {
    //z has no derivatives, so it is a scalar
    (*this) -= z.v_;
    return *this;
  }

  v_ -= z.v_;
  for (int i = 0; i < N_; i++)
    d_[i] -= z.d_[i];

  return *this;
}

ALWAYS_INLINE  SurrealD&
SurrealD::operator-=( const Real& r )
{
  v_ -= r;
  return *this;
}

ALWAYS_INLINE  SurrealD&
SurrealD::operator*=( const SurrealD& z )
{
  if ((N_ == 0) && (z.N_ > 0))
  {
    N_ = z.N_;
    d_ = new Real[N_]();    // NOTE: d_ is value-initialized here

  }
  else if ( z.N_ == 0 )
  {
    //z has no derivatives, so it is a scalar
    (*this) *= z.v_;
    return *this;
  }

  for (int i = 0; i < N_; i++)
    d_[i] = v_*z.d_[i] + d_[i]*z.v_;
  v_ *= z.v_;

  return *this;
}

ALWAYS_INLINE  SurrealD&
SurrealD::operator*=( const Real& r )
{
  for (int i = 0; i < N_; i++)
    d_[i] *= r;
  v_ *= r;
  return *this;
}

ALWAYS_INLINE  SurrealD&
SurrealD::operator/=( const SurrealD& z)
{
  if ((N_ == 0) && (z.N_ > 0))
  {
    N_ = z.N_;
    d_ = new Real[N_]();    // NOTE: d_ is value-initialized here
  }
  else if ( z.N_ == 0 )
  {
    //z has no derivatives, so it is a scalar
    (*this) /= z.v_;
    return *this;
  }


  Real tmp = 1./(z.v_*z.v_);
  for (int i = 0; i < N_; i++)
    d_[i] = (z.v_*d_[i] - v_*z.d_[i])*tmp;
  v_ /= z.v_;

  return *this;
}

ALWAYS_INLINE  SurrealD&
SurrealD::operator/=( const Real& r )
{
  Real tmp = 1./r;
  for (int i = 0; i < N_; i++)
    d_[i] *= tmp;
  v_ *= tmp;
  return *this;
}


// binary operators

ALWAYS_INLINE  SurrealD
operator+( const SurrealD& a, const SurrealD& b )
{
  if (a.N_ == 0 && b.N_ == 0)
    return SurrealD(a.v_ + b.v_);
  else if (a.N_ > 0 && b.N_ == 0)
    return SurrealD(a.v_ + b.v_, a.d_, a.N_);
  else if (a.N_ == 0 && b.N_ > 0)
    return SurrealD(a.v_ + b.v_, b.d_, b.N_);

  SANS_ASSERT( a.N_ == b.N_ );

  const int N = b.N_;
  Real d0[N];
  for (int i = 0; i < N; i++)
    d0[i] = a.d_[i] + b.d_[i];
  return SurrealD(a.v_ + b.v_, d0, N);
}

ALWAYS_INLINE  SurrealD
operator+( const SurrealD& a, const Real& b )
{
  if (a.N_ == 0 )
    return SurrealD(a.v_ + b);

  return SurrealD(a.v_ + b, a.d_, a.N_);
}

ALWAYS_INLINE  SurrealD
operator+( const Real& a, const SurrealD& b )
{
  if (b.N_ == 0 )
    return SurrealD(a + b.v_);

  return SurrealD(a + b.v_, b.d_, b.N_);
}

ALWAYS_INLINE  SurrealD
operator-( const SurrealD& a, const SurrealD& b )
{
  if (a.N_ == 0 && b.N_ == 0)
    return SurrealD(a.v_ - b.v_);
  else if (a.N_ > 0 && b.N_ == 0)
    return SurrealD(a.v_ - b.v_, a.d_, a.N_);
  else if (a.N_ == 0 && b.N_ > 0)
  {
    const int N = b.N_;
    Real d0[N];
    for (int i = 0; i < N; i++)
      d0[i] = -b.d_[i];
    return SurrealD(a.v_ - b.v_, d0, N);
  }

  SANS_ASSERT( a.N_ == b.N_ );

  const int N = b.N_;
  Real d0[N];
  for (int i = 0; i < N; i++)
    d0[i] = a.d_[i] - b.d_[i];
  return SurrealD(a.v_ - b.v_, d0, N);
}

ALWAYS_INLINE  SurrealD
operator-( const SurrealD& a, const Real& b )
{
  if (a.N_ == 0)
    return SurrealD(a.v_ - b);

  return SurrealD(a.v_ - b, a.d_, a.N_);
}

ALWAYS_INLINE  SurrealD
operator-( const Real& a, const SurrealD& b )
{
  if (b.N_ == 0)
    return SurrealD(a - b.v_);

  const int N = b.N_;
  Real d0[N];
  for (int i = 0; i < N; i++)
    d0[i] = -b.d_[i];
  return SurrealD(a - b.v_, d0, N);
}

ALWAYS_INLINE  SurrealD
operator*( const SurrealD& a, const SurrealD& b )
{
  if (a.N_ == 0 && b.N_ == 0)
    return SurrealD(a.v_ * b.v_);
  else if (a.N_ > 0 && b.N_ == 0)
  {
    const int N = a.N_;
    Real d0[N];
    for (int i = 0; i < N; i++)
      d0[i] = a.d_[i]*b.v_;
    return SurrealD(a.v_ * b.v_, d0, N);
  }
  else if (a.N_ == 0 && b.N_ > 0)
  {
    const int N = b.N_;
    Real d0[N];
    for (int i = 0; i < N; i++)
      d0[i] = a.v_*b.d_[i];
    return SurrealD(a.v_ * b.v_, d0, N);
  }

  const int N = b.N_;
  Real d0[N];
  for (int i = 0; i < N; i++)
    d0[i] = a.d_[i]*b.v_ + a.v_*b.d_[i];
  return SurrealD(a.v_ * b.v_, d0, N);
}

ALWAYS_INLINE  SurrealD
operator*( const SurrealD& a, const Real& b )
{
  if (a.N_ == 0)
    return SurrealD(a.v_ * b);

  const int N = a.N_;
  Real d0[N];
  for (int i = 0; i < N; i++)
    d0[i] = a.d_[i]*b;
  return SurrealD(a.v_ * b, d0, N);
}

ALWAYS_INLINE  SurrealD
operator*( const Real& a, const SurrealD& b )
{
  if (b.N_ == 0)
    return SurrealD(a * b.v_);

  const int N = b.N_;
  Real d0[N];
  for (int i = 0; i < N; i++)
    d0[i] = a*b.d_[i];
  return SurrealD(a * b.v_, d0, N);
}

ALWAYS_INLINE  SurrealD
operator/( const SurrealD& a, const SurrealD& b )
{
  Real tmp = 1./(b.v_*b.v_);

  if (a.N_ == 0 && b.N_ == 0)
    return SurrealD(a.v_ / b.v_);
  else if (a.N_ > 0 && b.N_ == 0)
  {
    const int N = a.N_;
    Real d0[N];
    for (int i = 0; i < N; i++)
      d0[i] = (b.v_*a.d_[i])*tmp;
    return SurrealD(a.v_ / b.v_, d0, N);
  }
  else if (a.N_ == 0 && b.N_ > 0)
  {
    const int N = b.N_;
    Real d0[N];
    for (int i = 0; i < N; i++)
      d0[i] = (-a.v_*b.d_[i])*tmp;
    return SurrealD(a.v_ / b.v_, d0, N);
  }

  SANS_ASSERT( a.N_ == b.N_ );

  const int N = b.N_;
  Real d0[N];
  for (int i = 0; i < N; i++)
    d0[i] = (b.v_*a.d_[i] - a.v_*b.d_[i])*tmp;
  return SurrealD(a.v_ / b.v_, d0, N);
}

ALWAYS_INLINE  SurrealD operator/( const SurrealD& a, const Real& b )
{
  if (a.N_ == 0)
    return SurrealD(a.v_ / b);

  Real tmp = 1./(b);
  const int N = a.N_;
  Real d0[N];
  for (int i = 0; i < N; i++)
    d0[i] = a.d_[i]*tmp;
  return SurrealD(a.v_ / b, d0, N);
}

ALWAYS_INLINE  SurrealD
operator/( const Real& a, const SurrealD& b )
{
  if (b.N_ == 0)
    return SurrealD(a / b.v_);

  const int N = b.N_;
  Real tmpv = a/(b.v_);
  Real tmpd = -1./(b.v_);
  Real d0[N];
  for (int i = 0; i < N; i++)
    d0[i] = tmpv*tmpd*b.d_[i];
  return SurrealD(tmpv, d0, N);
}


// relational operators

ALWAYS_INLINE  bool
operator==( const SurrealD& lhs, const SurrealD& rhs )
{
  return lhs.v_ == rhs.v_;
}

ALWAYS_INLINE  bool
operator==( const SurrealD& lhs, const Real& rhs )
{
  return lhs.v_ == rhs;
}

ALWAYS_INLINE  bool
operator==( const Real& lhs, const SurrealD& rhs )
{
  return lhs == rhs.v_;
}

ALWAYS_INLINE  bool
operator!=( const SurrealD& lhs, const SurrealD& rhs )
{
  return lhs.v_ != rhs.v_;
}

ALWAYS_INLINE  bool
operator!=( const SurrealD& lhs, const Real& rhs )
{
  return lhs.v_ != rhs;
}

ALWAYS_INLINE  bool
operator!=( const Real& lhs, const SurrealD& rhs )
{
  return lhs != rhs.v_;
}

ALWAYS_INLINE  bool
operator>( const SurrealD& lhs, const SurrealD& rhs )
{
  return lhs.v_ > rhs.v_;
}

ALWAYS_INLINE  bool
operator>( const SurrealD& lhs, const Real& rhs )
{
  return lhs.v_ > rhs;
}

ALWAYS_INLINE  bool
operator>( const Real& lhs, const SurrealD& rhs )
{
  return lhs > rhs.v_;
}

ALWAYS_INLINE  bool
operator<( const SurrealD& lhs, const SurrealD& rhs )
{
  return lhs.v_ < rhs.v_;
}

ALWAYS_INLINE  bool
operator<( const SurrealD& lhs, const Real& rhs )
{
  return lhs.v_ < rhs;
}

ALWAYS_INLINE  bool
operator<( const Real& lhs, const SurrealD& rhs )
{
  return lhs < rhs.v_;
}

ALWAYS_INLINE  bool
operator>=( const SurrealD& lhs, const SurrealD& rhs )
{
  return lhs.v_ >= rhs.v_;
}

ALWAYS_INLINE  bool
operator>=( const SurrealD& lhs, const Real& rhs )
{
  return lhs.v_ >= rhs;
}

ALWAYS_INLINE  bool
operator>=( const Real& lhs, const SurrealD& rhs )
{
  return lhs >= rhs.v_;
}

ALWAYS_INLINE  bool
operator<=( const SurrealD& lhs, const SurrealD& rhs )
{
  return lhs.v_ <= rhs.v_;
}

ALWAYS_INLINE  bool
operator<=( const SurrealD& lhs, const Real& rhs )
{
  return lhs.v_ <= rhs;
}

ALWAYS_INLINE  bool
operator<=( const Real& lhs, const SurrealD& rhs )
{
  return lhs <= rhs.v_;
}


// trig functions <cmath>

ALWAYS_INLINE  SurrealD
cos( const SurrealD& z )
{
  Real tmp = -std::sin(z.v_);
  Real d0[z.N_];
  for (int i = 0; i < z.N_; i++)
    d0[i] = tmp*z.d_[i];
  return SurrealD(std::cos(z.v_), d0, z.N_);
}

ALWAYS_INLINE  SurrealD
sin( const SurrealD& z )
{
  Real tmp = std::cos(z.v_);
  Real d0[z.N_];
  for (int i = 0; i < z.N_; i++)
    d0[i] = tmp*z.d_[i];
  return SurrealD(std::sin(z.v_), d0, z.N_);
}

ALWAYS_INLINE  SurrealD
tan( const SurrealD& z )
{
  Real cosv=std::cos(z.v_);
  Real tmp = 1./(cosv*cosv);
  Real d0[z.N_];
  for (int i = 0; i < z.N_; i++)
    d0[i] = tmp*z.d_[i];
  return SurrealD(std::tan(z.v_), d0, z.N_);
}

ALWAYS_INLINE  SurrealD
acos( const SurrealD& z )
{
  // derivative infinite at z.v_ = +/- 1.0
//  if ((z.v_ == -1) || (z.v_ == 1))
//  {
//    throw MathError( "acos arg = +/- 1" );
//  }

  Real tmp = -1./std::sqrt(1 - z.v_*z.v_);
  Real d0[z.N_];
  for (int i = 0; i < z.N_; i++)
    d0[i] = tmp*z.d_[i];
  return SurrealD(std::acos(z.v_), d0, z.N_);
}

ALWAYS_INLINE  SurrealD
asin( const SurrealD& z )
{
  // derivative infinite at z.v_ = +/- 1.0
//  if ((z.v_ == -1) || (z.v_ == 1))
//  {
//    throw MathError( "asin arg = +/- 1" );
//
//  }

  Real tmp = 1./std::sqrt(1 - z.v_*z.v_);
  Real d0[z.N_];
  for (int i = 0; i < z.N_; i++)
    d0[i] = tmp*z.d_[i];
  return SurrealD(std::asin(z.v_), d0, z.N_);
}

ALWAYS_INLINE  SurrealD
atan( const SurrealD& z )
{
  Real tmp = 1./(1 + z.v_*z.v_);
  Real d0[z.N_];
  for (int i = 0; i < z.N_; i++)
    d0[i] = tmp*z.d_[i];
  return SurrealD(std::atan(z.v_), d0, z.N_);
}

ALWAYS_INLINE  SurrealD
atan2( const SurrealD& z1, const SurrealD& z2)
{
  SANS_ASSERT(z1.N_ == z2.N_);
  Real tmp = 1./(z1.v_*z1.v_ + z2.v_*z2.v_);
  Real d0[z1.N_];
  for (int i = 0; i < z1.N_; i++)
    d0[i] = tmp*(z2.v_*z1.d_[i] - z1.v_*z2.d_[i]);
  return SurrealD(std::atan2(z1.v_,z2.v_), d0, z1.N_);
}


// hyperbolic functions <cmath>

ALWAYS_INLINE  SurrealD
cosh( const SurrealD& z )
{
  Real tmp = std::sinh(z.v_);
  Real d0[z.N_];
  for (int i = 0; i < z.N_; i++)
    d0[i] = tmp*z.d_[i];
  return SurrealD(std::cosh(z.v_), d0, z.N_);
}

ALWAYS_INLINE  SurrealD
sinh( const SurrealD& z )
{
  Real tmp = std::cosh(z.v_);
  Real d0[z.N_];
  for (int i = 0; i < z.N_; i++)
    d0[i] = tmp*z.d_[i];
  return SurrealD(std::sinh(z.v_), d0, z.N_);
}

ALWAYS_INLINE  SurrealD
tanh( const SurrealD& z )
{
  Real coshv=std::cosh(z.v_);
  Real tmp = 1./(coshv*coshv);
  Real d0[z.N_];
  for (int i = 0; i < z.N_; i++)
    d0[i] = tmp*z.d_[i];
  return SurrealD(std::tanh(z.v_), d0, z.N_);
}


// exp and log functions <cmath>

ALWAYS_INLINE  SurrealD
exp( const SurrealD& z )
{
  Real expv=std::exp(z.v_);
  Real d0[z.N_];
  for (int i = 0; i < z.N_; i++)
    d0[i] = expv*z.d_[i];
  return SurrealD(expv, d0, z.N_);
}

ALWAYS_INLINE  SurrealD
log( const SurrealD& z )
{
  Real tmp = 1./z.v_;
  Real d0[z.N_];
  for (int i = 0; i < z.N_; i++)
    d0[i] = tmp*z.d_[i];
  return SurrealD(std::log(z.v_), d0, z.N_);
}

ALWAYS_INLINE  SurrealD
log10( const SurrealD& z )
{
  Real tmp = 1./(z.v_*std::log(10.));
  Real d0[z.N_];
  for (int i = 0; i < z.N_; i++)
    d0[i] = tmp*z.d_[i];
  return SurrealD(std::log10(z.v_), d0, z.N_);
}


// power functions <cmath>

ALWAYS_INLINE  SurrealD
pow( const SurrealD& a, const SurrealD& b)
{
  SANS_ASSERT(a.N_ == b.N_);
  // many sticky points were derivative is undefined or infinite
  // badness if 0 <= b < 1 and a == 0
  Real powab=std::pow(a.v_,b.v_);
  Real tmp1 = b.v_*std::pow(a.v_, b.v_ - 1);
  Real tmp2 = powab*std::log(a.v_);
  Real d0[a.N_];
  for (int i = 0; i < a.N_; i++)
    d0[i] = tmp1*a.d_[i] + tmp2*b.d_[i];
  return SurrealD(powab, d0, a.N_);
}

ALWAYS_INLINE  SurrealD
pow( const SurrealD& a, const Real& b)
{
  Real tmp = b*std::pow(a.v_, b-1);
  Real d0[a.N_];
  for (int i = 0; i < a.N_; i++)
    d0[i] = tmp*a.d_[i];
  return SurrealD(std::pow(a.v_,b), d0, a.N_);
}

ALWAYS_INLINE  SurrealD
pow( const SurrealD& a, const int& b)
{
  Real tmp = static_cast<Real>(b)*std::pow(a.v_, static_cast<Real>(b - 1));
  Real d0[a.N_];
  for (int i = 0; i < a.N_; i++)
    d0[i] = tmp*a.d_[i];
  return SurrealD(std::pow(a.v_, static_cast<Real>(b)), d0, a.N_);
}

ALWAYS_INLINE  SurrealD
pow( const Real& a, const SurrealD& b)
{
  Real powab=std::pow(a, b.v_);
  Real tmp = powab*std::log(a);
  Real d0[b.N_];
  for (int i = 0; i < b.N_; i++)
    d0[i] = tmp*b.d_[i];
  return SurrealD(powab, d0, b.N_);
}

ALWAYS_INLINE  SurrealD
pow( const int& a, const SurrealD& b)
{
  Real powab=std::pow(static_cast<Real>(a), b.v_);
  Real tmp = powab*std::log(static_cast<Real>(a));
  Real d0[b.N_];
  for (int i = 0; i < b.N_; i++)
    d0[i] = tmp*b.d_[i];
  return SurrealD(powab, d0, b.N_);
}

ALWAYS_INLINE  SurrealD
sqrt( const SurrealD& z )
{
  Real sqrtv=std::sqrt(z.v_);
  if (sqrtv == 0)
  {
    return SurrealD(0., 0., z.N_);
  }
  else
  {
    Real tmp = 0.5/sqrtv;
    Real d0[z.N_];
    for (int i = 0; i < z.N_; i++)
      d0[i] = tmp*z.d_[i];
    return SurrealD(sqrtv, d0, z.N_);
  }
}


// rounding functions <cmath>

ALWAYS_INLINE  SurrealD
ceil( const SurrealD& z )
{
  Real d0[z.N_];
  for (int i = 0; i < z.N_; i++)
    d0[i] = 0;
  return SurrealD(std::ceil(z.v_), d0, z.N_);
}

ALWAYS_INLINE  SurrealD
floor( const SurrealD& z )
{
  Real d0[z.N_];
  for (int i = 0; i < z.N_; i++)
    d0[i] = 0;
  return SurrealD(std::floor(z.v_), d0, z.N_);
}


// misc functions <cmath>

ALWAYS_INLINE  SurrealD
abs( const SurrealD& z )
{
  return (z.v_ < 0) ? -z : z;
}

ALWAYS_INLINE  SurrealD
fabs( const SurrealD& z )
{
  return (z.v_ < 0) ? -z : z;
}


// output format: (v;d0,d1,d2,...,dN)

inline std::ostream&
operator<<( std::ostream& os, const SurrealD& z )
{
  os << '(' << z.v_;
  if (z.N_ > 0)
  {
    os << ';';
    for (int i = 0; i < z.N_ - 1; i++)
      os << z.d_[i] << ',';
    os << z.d_[z.N_ - 1];
  }
  os << ')';
  return os;
}

#endif // SURREALD_RVO_H
