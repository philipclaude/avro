// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef SURREALD_TRAD_H
#define SURREALD_TRAD_H

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
// Operators with traditional implementation
//
//----------------------------------------------------------------------------//

class SurrealD
{
public:
  //The default constructor is intentionally not included here. This means Surreal is
  //not initialized when declared, which is consistent with regular numbers. This also
  //improves performance.
  SurrealD( const SurrealD& z );
  explicit SurrealD( int n=0 );
  SurrealD( const Real& v0, const Real d0[], const int n );
  SurrealD( const Real& v0, const Real& d0, const int n );
  SurrealD( const int& v0, const int& d0, const int n );
  // cppcheck-suppress noExplicitConstructor
  SurrealD( const Real& v0 );
  SurrealD( const Real& v0, const int n );
  ~SurrealD();

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
  friend SurrealD log1p( const SurrealD& );

  // power functions <cmath>
  friend SurrealD pow( const SurrealD&, const SurrealD& );
  friend SurrealD pow( const SurrealD&, const Real& );
  friend SurrealD pow( const Real&, const SurrealD& );
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
  unsigned int N_;  // size of derivative array
};


// constructors

ALWAYS_INLINE
SurrealD::SurrealD( const SurrealD& z ) : v_(z.v_), d_(NULL), N_(z.N_)
{
  if (N_ > 0)
  {
    d_ = new Real[N_];
    for (unsigned int i = 0; i < N_; i++)
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
  for (unsigned int i = 0; i < N_; i++)
    d_[i] = d0[i];
}

ALWAYS_INLINE
SurrealD::SurrealD( const Real& v0, const Real& d0, const int n ) : v_(v0), d_(NULL), N_(n)
{
  if (N_ > 0)
  {
    d_ = new Real[N_];
    for (unsigned int i = 0; i < N_; i++)
      d_[i] = d0;
  }
}

ALWAYS_INLINE
SurrealD::SurrealD( const int& v0, const int& d0, const int n ) : v_(v0), d_(NULL), N_(n)
{
  if (N_ > 0)
  {
    d_ = new Real[N_];
    for (unsigned int i = 0; i < N_; i++)
      d_[i] = d0;
  }
}

ALWAYS_INLINE
SurrealD::SurrealD( const Real& v0 ) : v_(v0), d_(NULL), N_(0) {}

ALWAYS_INLINE
SurrealD::SurrealD( const Real& v0, const int n ) : v_(v0), d_(NULL), N_(n)
{
  d_ = new Real[N_];
}

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
  for (unsigned int i = 0; i < N_; i++)
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
    SurrealD c(-v_, N_);
    for (unsigned int i = 0; i < N_; i++)
      c.d_[i] = -d_[i];
    return c;
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
  for (unsigned int i = 0; i < N_; i++)
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
  for (unsigned int i = 0; i < N_; i++)
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

  for (unsigned int i = 0; i < N_; i++)
    d_[i] = v_*z.d_[i] + d_[i]*z.v_;
  v_ *= z.v_;

  return *this;
}

ALWAYS_INLINE  SurrealD&
SurrealD::operator*=( const Real& r )
{
  for (unsigned int i = 0; i < N_; i++)
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
  for (unsigned int i = 0; i < N_; i++)
    d_[i] = (z.v_*d_[i] - v_*z.d_[i])*tmp;
  v_ /= z.v_;

  return *this;
}

ALWAYS_INLINE  SurrealD&
SurrealD::operator/=( const Real& r )
{
  Real tmp = 1./r;
  for (unsigned int i = 0; i < N_; i++)
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

  const unsigned int N = b.N_;
  SurrealD c(a.v_ + b.v_, N);
  for (unsigned int i = 0; i < N; i++)
    c.d_[i] = a.d_[i] + b.d_[i];
  return c;
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
    const unsigned int N = b.N_;
    SurrealD c(a.v_ - b.v_, N);
    for (unsigned int i = 0; i < N; i++)
      c.d_[i] = -b.d_[i];
    return c;
  }

  SANS_ASSERT( a.N_ == b.N_ );

  const unsigned int N = b.N_;
  SurrealD c(a.v_ - b.v_, N);
  for (unsigned int i = 0; i < N; i++)
    c.d_[i] = a.d_[i] - b.d_[i];
  return c;
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

  const unsigned int N = b.N_;
  SurrealD c(a - b.v_, N);
  for (unsigned int i = 0; i < N; i++)
    c.d_[i] = -b.d_[i];
  return c;
}

ALWAYS_INLINE  SurrealD
operator*( const SurrealD& a, const SurrealD& b )
{
  if (a.N_ == 0 && b.N_ == 0)
    return SurrealD(a.v_ * b.v_);
  else if (a.N_ > 0 && b.N_ == 0)
  {
    const unsigned int N = a.N_;
    SurrealD c(a.v_ * b.v_, N);
    for (unsigned int i = 0; i < N; i++)
      c.d_[i] = a.d_[i]*b.v_;
    return c;
  }
  else if (a.N_ == 0 && b.N_ > 0)
  {
    const unsigned int N = b.N_;
    SurrealD c(a.v_ * b.v_, N);
    for (unsigned int i = 0; i < N; i++)
      c.d_[i] = a.v_*b.d_[i];
    return c;
  }

  const unsigned int N = b.N_;
  SurrealD c(a.v_ * b.v_, N);
  for (unsigned int i = 0; i < N; i++)
    c.d_[i] = a.d_[i]*b.v_ + a.v_*b.d_[i];
  return c;
}

ALWAYS_INLINE  SurrealD
operator*( const SurrealD& a, const Real& b )
{
  if (a.N_ == 0)
    return SurrealD(a.v_ * b);

  const unsigned int N = a.N_;
  SurrealD c(a.v_ * b, N);
  for (unsigned int i = 0; i < N; i++)
    c.d_[i] = a.d_[i]*b;
  return c;
}

ALWAYS_INLINE  SurrealD
operator*( const Real& a, const SurrealD& b )
{
  if (b.N_ == 0)
    return SurrealD(a * b.v_);

  const unsigned int N = b.N_;
  SurrealD c(a * b.v_, N);
  for (unsigned int i = 0; i < N; i++)
    c.d_[i] = a*b.d_[i];
  return c;
}

ALWAYS_INLINE  SurrealD
operator/( const SurrealD& a, const SurrealD& b )
{
  Real tmp = 1./(b.v_*b.v_);

  if (a.N_ == 0 && b.N_ == 0)
    return SurrealD(a.v_ / b.v_);
  else if (a.N_ > 0 && b.N_ == 0)
  {
    const unsigned int N = a.N_;
    SurrealD c(a.v_ / b.v_, N);
    for (unsigned int i = 0; i < N; i++)
      c.d_[i] = (b.v_*a.d_[i])*tmp;
    return c;
  }
  else if (a.N_ == 0 && b.N_ > 0)
  {
    const unsigned int N = b.N_;
    SurrealD c(a.v_ / b.v_, N);
    for (unsigned int i = 0; i < N; i++)
      c.d_[i] = (-a.v_*b.d_[i])*tmp;
    return c;
  }

  SANS_ASSERT( a.N_ == b.N_ );

  const unsigned int N = b.N_;
  SurrealD c(a.v_ / b.v_, N);
  for (unsigned int i = 0; i < N; i++)
    c.d_[i] = (b.v_*a.d_[i] - a.v_*b.d_[i])*tmp;
  return c;
}

ALWAYS_INLINE  SurrealD operator/( const SurrealD& a, const Real& b )
{
  if (a.N_ == 0)
    return SurrealD(a.v_ / b);

  Real tmp = 1./(b);
  const unsigned int N = a.N_;
  SurrealD c(a.v_ / b, N);
  for (unsigned int i = 0; i < N; i++)
    c.d_[i] = a.d_[i]*tmp;
  return c;
}

ALWAYS_INLINE  SurrealD
operator/( const Real& a, const SurrealD& b )
{
  if (b.N_ == 0)
    return SurrealD(a / b.v_);

  const unsigned int N = b.N_;
  Real tmpv = a/(b.v_);
  Real tmpd = -1./(b.v_);
  SurrealD c(tmpv, N);
  for (unsigned int i = 0; i < N; i++)
    c.d_[i] = tmpv*tmpd*b.d_[i];
  return c;
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


//Macros for functions

#define SURREALD_FUNC1( NAME, FUNC, DERIV ) \
ALWAYS_INLINE  SurrealD \
NAME( const SurrealD& z ) \
{ \
  if ( z.N_ == 0 ) \
    return SurrealD(FUNC); \
  \
  const unsigned int N = z.N_; \
  Real tmp = DERIV; \
  SurrealD c(FUNC, N); \
  for (unsigned int i = 0; i < N; i++) \
    c.d_[i] = tmp*z.d_[i]; \
  return c; \
}

#define SURREALD_FUNC2( NAME, FUNC, DERIV ) \
ALWAYS_INLINE  SurrealD \
NAME( const SurrealD& z1, const SurrealD& z2) \
{ \
  if ( z1.N_ == 0 && z2.N_ == 0 ) \
    return SurrealD(FUNC); \
  else if ( z1.N_ > 0 && z2.N_ == 0 ) \
  { \
    const unsigned int N = z1.N_; \
    Real tmp = DERIV; \
    SurrealD c(FUNC, N); \
    for (unsigned int i = 0; i < N; i++) \
      c.d_[i] = tmp*(z2.v_*z1.d_[i]); \
    return c; \
  } \
  else if ( z1.N_ == 0 && z2.N_ > 0 ) \
  { \
    const unsigned int N = z2.N_; \
    Real tmp = DERIV; \
    SurrealD c(FUNC, N); \
    for (unsigned int i = 0; i < N; i++) \
      c.d_[i] = tmp*(-z1.v_*z2.d_[i]); \
    return c; \
  } \
  \
  SANS_ASSERT( z1.N_ == z2.N_ ); \
  \
  const unsigned int N = z1.N_; \
  Real tmp = DERIV; \
  SurrealD c(FUNC, N); \
  for (unsigned int i = 0; i < N; i++) \
    c.d_[i] = tmp*(z2.v_*z1.d_[i] - z1.v_*z2.d_[i]); \
  return c; \
}

// trig functions <cmath>

SURREALD_FUNC1( cos, cos(z.v_), -sin(z.v_) )
SURREALD_FUNC1( sin, sin(z.v_),  cos(z.v_) )
SURREALD_FUNC1( tan, tan(z.v_),  Real(1)/(cos(z.v_)*cos(z.v_)) )
SURREALD_FUNC1( acos, acos(z.v_), -Real(1)/sqrt(1 - z.v_*z.v_) )
SURREALD_FUNC1( asin, asin(z.v_),  Real(1)/sqrt(1 - z.v_*z.v_) )
SURREALD_FUNC1( atan, atan(z.v_),  Real(1)/(1 + z.v_*z.v_) )

SURREALD_FUNC2( atan2, atan2(z1.v_, z2.v_),  Real(1)/(z1.v_*z1.v_ + z2.v_*z2.v_) )

// hyperbolic functions <cmath>

SURREALD_FUNC1( cosh, cosh(z.v_), sinh(z.v_) )
SURREALD_FUNC1( sinh, sinh(z.v_), cosh(z.v_) )
SURREALD_FUNC1( tanh, tanh(z.v_), Real(1)/(cosh(z.v_)*cosh(z.v_)) )

// exp and log functions <cmath>

SURREALD_FUNC1( exp, exp(z.v_), exp(z.v_) )
SURREALD_FUNC1( log, log(z.v_), Real(1)/z.v_ )
SURREALD_FUNC1( log10, log10(z.v_), Real(1)/(z.v_*log(10.)) )
SURREALD_FUNC1( log1p, log1p(z.v_), Real(1)/( 1 + z.v_ ) )


// power functions <cmath>

ALWAYS_INLINE  SurrealD
pow( const SurrealD& a, const SurrealD& b)
{
  Real powab=pow(a.v_,b.v_);

  if (a.N_ == 0 && b.N_ == 0)
    return SurrealD(powab);
  else if (a.N_ > 0 && b.N_ == 0)
  {
    Real tmp1 = (a.v_ == 0) ? ((b.v_ == 1) ? 1 : 0) : b.v_*pow(a.v_, b.v_ - 1);

    const unsigned int N = a.N_;
    SurrealD c(powab, N);
    for (unsigned int i = 0; i < N; i++)
      c.d_[i] = tmp1*a.d_[i];
    return c;
  }
  else if (a.N_ == 0 && b.N_ > 0)
  {
    Real tmp2 = (a.v_ == 0) ? 0 : powab*log(a.v_);

    const unsigned int N = b.N_;
    SurrealD c(powab, N);
    for (unsigned int i = 0; i < N; i++)
      c.d_[i] = tmp2*b.d_[i];
    return c;
  }

  SANS_ASSERT(a.N_ == b.N_);
  // many sticky points were derivative is undefined or infinite
  // badness if 0 <= b < 1 and a == 0
  Real tmp1 = (a.v_ == 0) ? ((b.v_ == 1) ? 1 : 0) : b.v_*pow(a.v_, b.v_ - 1);
  Real tmp2 = (a.v_ == 0) ? 0 : powab*log(a.v_);

  SurrealD c(powab, a.N_);
  for (unsigned int i = 0; i < a.N_; i++)
    c.d_[i] = tmp1*a.d_[i] + tmp2*b.d_[i];
  return c;
}

ALWAYS_INLINE  SurrealD
pow( const SurrealD& a, const Real& b)
{
  Real powab=pow(a.v_,b);

  if (a.N_ == 0)
    return SurrealD(powab);

  Real tmp = (a.v_ == 0) ? ((b == 1) ? 1 : 0) : b*pow(a.v_, b - 1);
  SurrealD c(powab, a.N_);
  for (unsigned int i = 0; i < a.N_; i++)
    c.d_[i] = tmp*a.d_[i];
  return c;
}

ALWAYS_INLINE  SurrealD
pow( const Real& a, const SurrealD& b)
{
  Real powab=pow(a, b.v_);

  if (b.N_ == 0)
    return SurrealD(powab);

  Real tmp = (a == 0) ? 0 : powab*log(a);
  SurrealD c(powab, b.N_);
  for (unsigned int i = 0; i < b.N_; i++)
    c.d_[i] = tmp*b.d_[i];
  return c;
}

ALWAYS_INLINE  SurrealD
sqrt( const SurrealD& z )
{
  Real sqrtv=sqrt(z.v_);
  if (sqrtv == 0)
  {
    return SurrealD(0., 0., z.N_);
  }
  else
  {
    Real tmp = 0.5/sqrtv;
    SurrealD c(sqrtv, z.N_);
    for (unsigned int i = 0; i < z.N_; i++)
      c.d_[i] = tmp*z.d_[i];
    return c;
  }
}


// rounding functions <cmath>

ALWAYS_INLINE  SurrealD
ceil( const SurrealD& z )
{
  if (z.N_ == 0)
    return SurrealD(ceil(z.v_));

  return SurrealD(ceil(z.v_), 0., z.N_);
}

ALWAYS_INLINE  SurrealD
floor( const SurrealD& z )
{
  if (z.N_ == 0)
    return SurrealD(floor(z.v_));

  return SurrealD(floor(z.v_), 0., z.N_);
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

//Clean up macro definitions
#undef SURREALD_FUNC1
#undef SURREALD_FUNC2


// output format: (v;d0,d1,d2,...,dN)

inline std::ostream&
operator<<( std::ostream& os, const SurrealD& z )
{
  os << '(' << z.v_;
  if (z.N_ > 0)
  {
    os << ';';
    for (unsigned int i = 0; i < z.N_ - 1; i++)
      os << z.d_[i] << ',';
    os << z.d_[z.N_ - 1];
  }
  os << ')';
  return os;
}

#endif // SURREALD_TRAD_H
