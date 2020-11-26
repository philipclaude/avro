// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef SURREALS_RVO_H
#define SURREALS_RVO_H

//  overloaded derivative operator
//  ref: derivify.h (Google it)

#include <cmath>
#include <iostream>
#include <string>

#include "tools/SANSnumerics.h"
#include "tools/SANSException.h"
#include "tools/AlignMem.h"
#include "tools/always_inline.h"

//----------------------------------------------------------------------------//
// SurrealS<N>:  value, N derivatives
//
// Operators with Return Value Optimization (RVO)
//
// statically defined derivative array
//----------------------------------------------------------------------------//

template<int N>
class SurrealS
{
public:
  //The default constructor is intentionally empty here. This means Surreal is
  //not initialized when declared, which is consistent with regular numbers. This also
  //improves performance.
  SurrealS() {}
  SurrealS( const SurrealS& z );
  SurrealS( const Real& v0, const Real d0[], int n );
  SurrealS( const Real& v0, const Real& d0 );
  SurrealS( const int& v0, const int& d0 );
  explicit SurrealS( const Real& v0 );

  ~SurrealS() {}

#if 0
  operator Real() const { return v_; }
  operator int() const { return int(v_); }
#endif

  int size() const { return N; }

  // value accessor operators
  Real& value()       { return v_; }
  Real  value() const { return v_; }

  // derivative accessor operators
  Real& deriv( int i=0 )       { return d_[i]; }
  Real  deriv( int i=0 ) const { return d_[i]; }

  // assignment
  SurrealS& operator=( const SurrealS& );
  SurrealS& operator=( const Real& );

  // unary operators; no side effects
  const SurrealS& operator+() const;
  const SurrealS  operator-() const;

  // binary accumulation operators
  SurrealS& operator+=( const SurrealS& );
  SurrealS& operator+=( const Real& );
  SurrealS& operator-=( const SurrealS& );
  SurrealS& operator-=( const Real& );
  SurrealS& operator*=( const SurrealS& );
  SurrealS& operator*=( const Real& );
  SurrealS& operator/=( const SurrealS& );
  SurrealS& operator/=( const Real& );

  // binary operators
  template<int M> friend SurrealS<M> operator+( const SurrealS<M>&, const SurrealS<M>& );
  template<int M> friend SurrealS<M> operator+( const SurrealS<M>&, const Real& );
  template<int M> friend SurrealS<M> operator+( const Real&, const SurrealS<M>& );
  template<int M> friend SurrealS<M> operator-( const SurrealS<M>&, const SurrealS<M>& );
  template<int M> friend SurrealS<M> operator-( const SurrealS<M>&, const Real& );
  template<int M> friend SurrealS<M> operator-( const Real&, const SurrealS<M>& );
  template<int M> friend SurrealS<M> operator*( const SurrealS<M>&, const SurrealS<M>& );
  template<int M> friend SurrealS<M> operator*( const SurrealS<M>&, const Real& );
  template<int M> friend SurrealS<M> operator*( const Real&, const SurrealS<M>& );
  template<int M> friend SurrealS<M> operator/( const SurrealS<M>&, const SurrealS<M>& );
  template<int M> friend SurrealS<M> operator/( const SurrealS<M>&, const Real& );
  template<int M> friend SurrealS<M> operator/( const Real&, const SurrealS<M>& );

  // relational operators
  template<int M> friend bool operator==( const SurrealS<M>&, const SurrealS<M>& );
  template<int M> friend bool operator==( const SurrealS<M>&, const Real& );
  template<int M> friend bool operator==( const Real&, const SurrealS<M>& );
  template<int M> friend bool operator!=( const SurrealS<M>&, const SurrealS<M>& );
  template<int M> friend bool operator!=( const SurrealS<M>&, const Real& );
  template<int M> friend bool operator!=( const Real&, const SurrealS<M>& );
  template<int M> friend bool operator>( const SurrealS<M>&, const SurrealS<M>& );
  template<int M> friend bool operator>( const SurrealS<M>&, const Real& );
  template<int M> friend bool operator>( const Real&, const SurrealS<M>& );
  template<int M> friend bool operator<( const SurrealS<M>&, const SurrealS<M>& );
  template<int M> friend bool operator<( const SurrealS<M>&, const Real& );
  template<int M> friend bool operator<( const Real&, const SurrealS<M>& );
  template<int M> friend bool operator>=( const SurrealS<M>&, const SurrealS<M>& );
  template<int M> friend bool operator>=( const SurrealS<M>&, const Real& );
  template<int M> friend bool operator>=( const Real&, const SurrealS<M>& );
  template<int M> friend bool operator<=( const SurrealS<M>&, const SurrealS<M>& );
  template<int M> friend bool operator<=( const SurrealS<M>&, const Real& );
  template<int M> friend bool operator<=( const Real&, const SurrealS<M>& );

  // trig functions <cmath>
  template<int M> friend SurrealS<M> cos( const SurrealS<M>& );
  template<int M> friend SurrealS<M> sin( const SurrealS<M>& );
  template<int M> friend SurrealS<M> tan( const SurrealS<M>& );
  template<int M> friend SurrealS<M> acos( const SurrealS<M>& );
  template<int M> friend SurrealS<M> asin( const SurrealS<M>& );
  template<int M> friend SurrealS<M> atan( const SurrealS<M>& );
  template<int M> friend SurrealS<M> atan2( const SurrealS<M>&, const SurrealS<M>& );

  // hyperbolic functions <cmath>
  template<int M> friend SurrealS<M> cosh( const SurrealS<M>& );
  template<int M> friend SurrealS<M> sinh( const SurrealS<M>& );
  template<int M> friend SurrealS<M> tanh( const SurrealS<M>& );

  // exp and log functions <cmath>
  template<int M> friend SurrealS<M> exp( const SurrealS<M>& );
  template<int M> friend SurrealS<M> log( const SurrealS<M>& );
  template<int M> friend SurrealS<M> log10( const SurrealS<M>& );

  // power functions <cmath>
  template<int M> friend SurrealS<M> pow( const SurrealS<M>&, const SurrealS<M>& );
  template<int M> friend SurrealS<M> pow( const SurrealS<M>&, const Real& );
  template<int M> friend SurrealS<M> pow( const SurrealS<M>&, const int& );
  template<int M> friend SurrealS<M> pow( const Real&, const SurrealS<M>& );
  template<int M> friend SurrealS<M> pow( const int&, const SurrealS<M>& );
  template<int M> friend SurrealS<M> sqrt( const SurrealS<M>& );

  // rounding functions <cmath>
  template<int M> friend SurrealS<M> ceil( const SurrealS<M>& );
  template<int M> friend SurrealS<M> floor( const SurrealS<M>& );

  // misc functions <cmath>
  template<int M> friend SurrealS<M> abs( const SurrealS<M>& );
  template<int M> friend SurrealS<M> fabs( const SurrealS<M>& );

  // classification functions <cmath>
  template<int M> friend bool isfinite( const SurrealS<M>& );
  template<int M> friend bool isinf( const SurrealS<M>& );
  template<int M> friend bool isnan( const SurrealS<M>& );

  // input/output
  template<int M> friend std::istream& operator>>( std::istream&, SurrealS<M>& );
  template<int M> friend std::ostream& operator<<( std::ostream&, const SurrealS<M>& );

  void dump( int indentSize=0 ) const;

private:
  ALIGN_MEM Real v_;      // value
  ALIGN_MEM Real d_[N];   // derivative array
};




// constructors

template<int N>
ALWAYS_INLINE
SurrealS<N>::SurrealS( const SurrealS& z ) : v_(z.v_)
{
  for (int i = 0; i < N; i++)
    d_[i] = z.d_[i];
}

template<int N>
ALWAYS_INLINE
SurrealS<N>::SurrealS( const Real& v0, const Real d0[], int n ) : v_(v0)
{
  SANS_ASSERT( n == N );
  for (int i = 0; i < N; i++)
    d_[i] = d0[i];
}

template<int N>
ALWAYS_INLINE
SurrealS<N>::SurrealS( const Real& v0, const Real& d0 ) : v_(v0)
{
  for (int i = 0; i < N; i++)
    d_[i] = d0;
}

template<int N>
ALWAYS_INLINE
SurrealS<N>::SurrealS( const int& v0, const int& d0 ) : v_(v0)
{
  for (int i = 0; i < N; i++)
    d_[i] = d0;
}

template<int N>
ALWAYS_INLINE
SurrealS<N>::SurrealS( const Real& v0 ) : v_(v0)
{
  for (int i = 0; i < N; i++)
    d_[i] = 0;
}


// assignment

template<int N>
ALWAYS_INLINE  SurrealS<N>&
SurrealS<N>::operator=( const SurrealS<N>& z )
{
  //Do nothing if assigning self to self
  if ( &z == this ) return *this;

  v_ = z.v_;
  for (int i = 0; i < N; i++)
    d_[i] = z.d_[i];
  return *this;
}

template<int N>
ALWAYS_INLINE  SurrealS<N>&
SurrealS<N>::operator=( const Real& r )
{
  v_ = r;
  for (int i = 0; i < N; i++)
    d_[i] = 0;
  return *this;
}


// unary operators; no side effects

template<int N>
ALWAYS_INLINE  const SurrealS<N>&
SurrealS<N>::operator+() const
{
  return *this;
}

template<int N>
ALWAYS_INLINE  const SurrealS<N>
SurrealS<N>::operator-() const
{
  Real d0[N];
  for (int i = 0; i < N; i++)
    d0[i] = -d_[i];
  return SurrealS<N>(-v_, d0, N);
}


// binary accumulation operators

template<int N>
ALWAYS_INLINE  SurrealS<N>&
SurrealS<N>::operator+=( const SurrealS<N>& z )
{
  v_ += z.v_;
  for (int i = 0; i < N; i++)
    d_[i] += z.d_[i];
  return *this;
}

template<int N>
ALWAYS_INLINE  SurrealS<N>&
SurrealS<N>::operator+=( const Real& r )
{
  v_ += r;
  return *this;
}

template<int N>
ALWAYS_INLINE  SurrealS<N>&
SurrealS<N>::operator-=( const SurrealS<N>& z )
{
  v_ -= z.v_;
  for (int i = 0; i < N; i++)
    d_[i] -= z.d_[i];
  return *this;
}

template<int N>
ALWAYS_INLINE  SurrealS<N>&
SurrealS<N>::operator-=( const Real& r )
{
  v_ -= r;
  return *this;
}


template<int N>
ALWAYS_INLINE  SurrealS<N>&
SurrealS<N>::operator*=( const SurrealS<N>& z )
{
  for (int i = 0; i < N; i++)
    d_[i] = v_*z.d_[i] + d_[i]*z.v_;
  v_ *= z.v_;
  return *this;
}

template<int N>
ALWAYS_INLINE  SurrealS<N>&
SurrealS<N>::operator*=( const Real& r )
{
  for (int i = 0; i < N; i++)
    d_[i] *= r;
  v_ *= r;
  return *this;
}

template<int N>
ALWAYS_INLINE  SurrealS<N>&
SurrealS<N>::operator/=( const SurrealS<N>& z)
{
  Real tmp = 1./(z.v_*z.v_);
  for (int i = 0; i < N; i++)
    d_[i] = (z.v_*d_[i] - v_*z.d_[i])*tmp;
  v_ /= z.v_;
  return *this;
}

template<int N>
ALWAYS_INLINE  SurrealS<N>&
SurrealS<N>::operator/=( const Real& r )
{
  Real tmp = 1./r;
  for (int i = 0; i < N; i++)
    d_[i] *= tmp;
  v_ *= tmp;
  return *this;
}


// debug dump of private data
template<int N>
void
SurrealS<N>::dump( int indentSize ) const
{
  std::string indent(indentSize, ' ');
  std::cout << indent << "SurrealS<" << N << ">: v_ = " << v_;
  std::cout << "  d_[" << N << "] = (";
  for (int n = 0; n < N-1; n++)
    std::cout << d_[n] << ",";
  std::cout << d_[N-1] << ")" << std::endl;
}


// binary operators

template<int N>
ALWAYS_INLINE  SurrealS<N>
operator+( const SurrealS<N>& a, const SurrealS<N>& b )
{
  Real d0[N];
  for (int i = 0; i < N; i++)
    d0[i] = a.d_[i] + b.d_[i];
  return SurrealS<N>(a.v_ + b.v_, d0, N);
}

template<int N>
ALWAYS_INLINE  SurrealS<N>
operator+( const SurrealS<N>& a, const Real& b )
{
  Real d0[N];
  for (int i = 0; i < N; i++)
    d0[i] = a.d_[i];
  return SurrealS<N>(a.v_ + b, d0, N);

}

template<int N>
ALWAYS_INLINE  SurrealS<N>
operator+( const Real& a, const SurrealS<N>& b )
{
  Real d0[N];
  for (int i = 0; i < N; i++)
    d0[i] = b.d_[i];
  return SurrealS<N>(a + b.v_, d0, N);
}

template<int N>
ALWAYS_INLINE  SurrealS<N>
operator-( const SurrealS<N>& a, const SurrealS<N>& b )
{
  Real d0[N];
  for (int i = 0; i < N; i++)
    d0[i] = a.d_[i] - b.d_[i];
  return SurrealS<N>(a.v_ - b.v_, d0, N);
}

template<int N>
ALWAYS_INLINE  SurrealS<N>
operator-( const SurrealS<N>& a, const Real& b )
{
  Real d0[N];
  for (int i = 0; i < N; i++)
    d0[i] = a.d_[i];
  return SurrealS<N>(a.v_ - b, d0, N);
}

template<int N>
ALWAYS_INLINE  SurrealS<N>
operator-( const Real& a, const SurrealS<N>& b )
{
  Real d0[N];
  for (int i = 0; i < N; i++)
    d0[i] = -b.d_[i];
  return SurrealS<N>(a - b.v_, d0, N);
}

template<int N>
ALWAYS_INLINE  SurrealS<N>
operator*( const SurrealS<N>& a, const SurrealS<N>& b )
{
  Real d0[N];
  for (int i = 0; i < N; i++)
    d0[i] = a.d_[i]*b.v_ + a.v_*b.d_[i];
  return SurrealS<N>(a.v_ * b.v_, d0, N);
}

template<int N>
ALWAYS_INLINE  SurrealS<N>
operator*( const SurrealS<N>& a, const Real& b )
{
  Real d0[N];
  for (int i = 0; i < N; i++)
    d0[i] = a.d_[i]*b;
  return SurrealS<N>(a.v_ * b, d0, N);
}

template<int N>
ALWAYS_INLINE  SurrealS<N>
operator*( const Real& a, const SurrealS<N>& b )
{
  Real d0[N];
  for (int i = 0; i < N; i++)
    d0[i] = a*b.d_[i];
  return SurrealS<N>(a * b.v_, d0, N);
}

template<int N>
ALWAYS_INLINE  SurrealS<N>
operator/( const SurrealS<N>& a, const SurrealS<N>& b )
{
  Real tmp = 1./(b.v_*b.v_);
  Real d0[N];
  for (int i = 0; i < N; i++)
    d0[i] = (b.v_*a.d_[i] - a.v_*b.d_[i])*tmp;
  return SurrealS<N>(a.v_ / b.v_, d0, N);
}

template<int N>
ALWAYS_INLINE  SurrealS<N> operator/( const SurrealS<N>& a, const Real& b )
{
  Real tmp = 1./(b);
  Real d0[N];
  for (int i = 0; i < N; i++)
    d0[i] = a.d_[i]*tmp;
  return SurrealS<N>(a.v_ / b, d0, N);
}

template<int N>
ALWAYS_INLINE  SurrealS<N>
operator/( const Real& a, const SurrealS<N>& b )
{
  Real tmpv = a/(b.v_);
  Real tmpd = -1./(b.v_);
  Real d0[N];
  for (int i = 0; i < N; i++)
    d0[i] = tmpv*tmpd*b.d_[i];
  return SurrealS<N>(tmpv, d0, N);
}


// relational operators

template<int N>
ALWAYS_INLINE  bool
operator==( const SurrealS<N>& lhs, const SurrealS<N>& rhs )
{
  return lhs.v_ == rhs.v_;
}

template<int N>
ALWAYS_INLINE  bool
operator==( const SurrealS<N>& lhs, const Real& rhs )
{
  return lhs.v_ == rhs;
}

template<int N>
ALWAYS_INLINE  bool
operator==( const Real& lhs, const SurrealS<N>& rhs )
{
  return lhs == rhs.v_;
}

template<int N>
ALWAYS_INLINE  bool
operator!=( const SurrealS<N>& lhs, const SurrealS<N>& rhs )
{
  return lhs.v_ != rhs.v_;
}

template<int N>
ALWAYS_INLINE  bool
operator!=( const SurrealS<N>& lhs, const Real& rhs )
{
  return lhs.v_ != rhs;
}

template<int N>
ALWAYS_INLINE  bool
operator!=( const Real& lhs, const SurrealS<N>& rhs )
{
  return lhs != rhs.v_;
}

template<int N>
ALWAYS_INLINE  bool
operator>( const SurrealS<N>& lhs, const SurrealS<N>& rhs )
{
  return lhs.v_ > rhs.v_;
}

template<int N>
ALWAYS_INLINE  bool
operator>( const SurrealS<N>& lhs, const Real& rhs )
{
  return lhs.v_ > rhs;
}

template<int N>
ALWAYS_INLINE  bool
operator>( const Real& lhs, const SurrealS<N>& rhs )
{
  return lhs > rhs.v_;
}

template<int N>
ALWAYS_INLINE  bool
operator<( const SurrealS<N>& lhs, const SurrealS<N>& rhs )
{
  return lhs.v_ < rhs.v_;
}

template<int N>
ALWAYS_INLINE  bool
operator<( const SurrealS<N>& lhs, const Real& rhs )
{
  return lhs.v_ < rhs;
}

template<int N>
ALWAYS_INLINE  bool
operator<( const Real& lhs, const SurrealS<N>& rhs )
{
  return lhs < rhs.v_;
}

template<int N>
ALWAYS_INLINE  bool
operator>=( const SurrealS<N>& lhs, const SurrealS<N>& rhs )
{
  return lhs.v_ >= rhs.v_;
}

template<int N>
ALWAYS_INLINE  bool
operator>=( const SurrealS<N>& lhs, const Real& rhs )
{
  return lhs.v_ >= rhs;
}

template<int N>
ALWAYS_INLINE  bool
operator>=( const Real& lhs, const SurrealS<N>& rhs )
{
  return lhs >= rhs.v_;
}

template<int N>
ALWAYS_INLINE  bool
operator<=( const SurrealS<N>& lhs, const SurrealS<N>& rhs )
{
  return lhs.v_ <= rhs.v_;
}

template<int N>
ALWAYS_INLINE  bool
operator<=( const SurrealS<N>& lhs, const Real& rhs )
{
  return lhs.v_ <= rhs;
}

template<int N>
ALWAYS_INLINE  bool
operator<=( const Real& lhs, const SurrealS<N>& rhs )
{
  return lhs <= rhs.v_;
}


// trig functions <cmath>

template<int N>
ALWAYS_INLINE  SurrealS<N>
cos( const SurrealS<N>& z )
{
  Real tmp = -std::sin(z.v_);
  Real d0[N];
  for (int i = 0; i < N; i++)
    d0[i] = tmp*z.d_[i];
  return SurrealS<N>(std::cos(z.v_), d0, N);
}

template<int N>
ALWAYS_INLINE  SurrealS<N>
sin( const SurrealS<N>& z )
{
  Real tmp = std::cos(z.v_);
  Real d0[N];
  for (int i = 0; i < N; i++)
    d0[i] = tmp*z.d_[i];
  return SurrealS<N>(std::sin(z.v_), d0, N);
}

template<int N>
ALWAYS_INLINE  SurrealS<N>
tan( const SurrealS<N>& z )
{
  Real cosv=std::cos(z.v_);
  Real tmp = 1./(cosv*cosv);
  Real d0[N];
  for (int i = 0; i < N; i++)
    d0[i] = tmp*z.d_[i];
  return SurrealS<N>(std::tan(z.v_), d0, N);
}

template<int N>
ALWAYS_INLINE  SurrealS<N>
acos( const SurrealS<N>& z )
{
  // derivative infinite at z.v_ = +/- 1.0
//  if ((z.v_ == -1) || (z.v_ == 1))
//  {
//    throw MathError( "acos arg = +/- 1" );
//  }

  Real tmp = -1./std::sqrt(1 - z.v_*z.v_);
  Real d0[N];
  for (int i = 0; i < N; i++)
    d0[i] = tmp*z.d_[i];
  return SurrealS<N>(std::acos(z.v_), d0, N);
}

template<int N>
ALWAYS_INLINE  SurrealS<N>
asin( const SurrealS<N>& z )
{
  // derivative infinite at z.v_ = +/- 1.0
//  if ((z.v_ == -1) || (z.v_ == 1))
//  {
//    throw MathError( "asin arg = +/- 1" );
//
//  }

  Real tmp = 1./std::sqrt(1 - z.v_*z.v_);
  Real d0[N];
  for (int i = 0; i < N; i++)
    d0[i] = tmp*z.d_[i];
  return SurrealS<N>(std::asin(z.v_), d0, N);
}

template<int N>
ALWAYS_INLINE  SurrealS<N>
atan( const SurrealS<N>& z )
{
  Real tmp = 1./(1 + z.v_*z.v_);
  Real d0[N];
  for (int i = 0; i < N; i++)
    d0[i] = tmp*z.d_[i];
  return SurrealS<N>(std::atan(z.v_), d0, N);
}

template<int N>
ALWAYS_INLINE  SurrealS<N>
atan2( const SurrealS<N>& z1, const SurrealS<N>& z2)
{
  Real tmp = 1./(z1.v_*z1.v_ + z2.v_*z2.v_);
  Real d0[N];
  for (int i = 0; i < N; i++)
    d0[i] = tmp*(z2.v_*z1.d_[i] - z1.v_*z2.d_[i]);
  return SurrealS<N>(std::atan2(z1.v_,z2.v_), d0, N);
}


// hyperbolic functions <cmath>

template<int N>
ALWAYS_INLINE  SurrealS<N>
cosh( const SurrealS<N>& z )
{
  Real tmp = std::sinh(z.v_);
  Real d0[N];
  for (int i = 0; i < N; i++)
    d0[i] = tmp*z.d_[i];
  return SurrealS<N>(std::cosh(z.v_), d0, N);
}

template<int N>
ALWAYS_INLINE  SurrealS<N>
sinh( const SurrealS<N>& z )
{
  Real tmp = std::cosh(z.v_);
  Real d0[N];
  for (int i = 0; i < N; i++)
    d0[i] = tmp*z.d_[i];
  return SurrealS<N>(std::sinh(z.v_), d0, N);
}

template<int N>
ALWAYS_INLINE  SurrealS<N>
tanh( const SurrealS<N>& z )
{
  Real coshv=std::cosh(z.v_);
  Real tmp = 1./(coshv*coshv);
  Real d0[N];
  for (int i = 0; i < N; i++)
    d0[i] = tmp*z.d_[i];
  return SurrealS<N>(std::tanh(z.v_), d0, N);
}


// exp and log functions <cmath>

template<int N>
ALWAYS_INLINE  SurrealS<N>
exp( const SurrealS<N>& z )
{
  Real expv=std::exp(z.v_);
  Real d0[N];
  for (int i = 0; i < N; i++)
    d0[i] = expv*z.d_[i];
  return SurrealS<N>(expv, d0, N);
}

template<int N>
ALWAYS_INLINE  SurrealS<N>
log( const SurrealS<N>& z )
{
  Real tmp = 1./z.v_;
  Real d0[N];
  for (int i = 0; i < N; i++)
    d0[i] = tmp*z.d_[i];
  return SurrealS<N>(std::log(z.v_), d0, N);
}

template<int N>
ALWAYS_INLINE  SurrealS<N>
log10( const SurrealS<N>& z )
{
  Real tmp = 1./(z.v_*std::log(10.));
  Real d0[N];
  for (int i = 0; i < N; i++)
    d0[i] = tmp*z.d_[i];
  return SurrealS<N>(std::log10(z.v_), d0, N);
}


// power functions <cmath>

template<int N>
ALWAYS_INLINE  SurrealS<N>
pow( const SurrealS<N>& a, const SurrealS<N>& b)
{
  // many sticky points were derivative is undefined or infinite
  // badness if 0 <= b < 1 and a == 0
  Real powab=std::pow(a.v_,b.v_);
  Real tmp1 = b.v_*std::pow(a.v_, b.v_ - 1);
  Real tmp2 = powab*std::log(a.v_);
  Real d0[N];
  for (int i = 0; i < N; i++)
    d0[i] = tmp1*a.d_[i] + tmp2*b.d_[i];
  return SurrealS<N>(powab, d0, N);
}

template<int N>
ALWAYS_INLINE  SurrealS<N>
pow( const SurrealS<N>& a, const Real& b)
{
  Real tmp = b*std::pow(a.v_, b-1);
  Real d0[N];
  for (int i = 0; i < N; i++)
    d0[i] = tmp*a.d_[i];
  return SurrealS<N>(std::pow(a.v_,b), d0, N);
}

template<int N>
ALWAYS_INLINE  SurrealS<N>
pow( const SurrealS<N>& a, const int& b)
{
  Real tmp = static_cast<Real>(b)*std::pow(a.v_, static_cast<Real>(b - 1));
  Real d0[N];
  for (int i = 0; i < N; i++)
    d0[i] = tmp*a.d_[i];
  return SurrealS<N>(std::pow(a.v_, static_cast<Real>(b)), d0, N);
}

template<int N>
ALWAYS_INLINE  SurrealS<N>
pow( const Real& a, const SurrealS<N>& b)
{
  Real powab=std::pow(a, b.v_);
  Real tmp = powab*std::log(a);
  Real d0[N];
  for (int i = 0; i < N; i++)
    d0[i] = tmp*b.d_[i];
  return SurrealS<N>(powab, d0, N);
}

template<int N>
ALWAYS_INLINE  SurrealS<N>
pow( const int& a, const SurrealS<N>& b)
{
  Real powab=std::pow(static_cast<Real>(a), b.v_);
  Real tmp = powab*std::log(static_cast<Real>(a));
  Real d0[N];
  for (int i = 0; i < N; i++)
    d0[i] = tmp*b.d_[i];
  return SurrealS<N>(powab, d0, N);
}

template<int N>
ALWAYS_INLINE  SurrealS<N>
sqrt( const SurrealS<N>& z )
{
  Real sqrtv=std::sqrt(z.v_);
  if (sqrtv == 0)
  {
    return SurrealS<N>(0, 0);
  }
  else
  {
    Real tmp = 0.5/sqrtv;
    Real d0[N];
    for (int i = 0; i < N; i++)
      d0[i] = tmp*z.d_[i];
    return SurrealS<N>(sqrtv, d0, N);
  }
}


// rounding functions <cmath>

template<int N>
ALWAYS_INLINE  SurrealS<N>
ceil( const SurrealS<N>& z )
{
  Real d0[N];
  for (int i = 0; i < N; i++)
    d0[i] = 0;
  return SurrealS<N>(std::ceil(z.v_), d0, N);
}

template<int N>
ALWAYS_INLINE  SurrealS<N>
floor( const SurrealS<N>& z )
{
  Real d0[N];
  for (int i = 0; i < N; i++)
    d0[i] = 0;
  return SurrealS<N>(std::floor(z.v_), d0, N);
}


// misc functions <cmath>

template<int N>
ALWAYS_INLINE  SurrealS<N>
abs( const SurrealS<N>& z )
{
  return (z.v_ < 0) ? -z : z;
}

template<int N>
ALWAYS_INLINE  SurrealS<N>
fabs( const SurrealS<N>& z )
{
  return (z.v_ < 0) ? -z : z;
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

template<int N>
std::ostream&
operator<<( std::ostream& os, const SurrealS<N>& z )
{
  os << '(' << z.value() << ';';
  for (int i = 0; i < N - 1; i++)
    os << z.deriv(i) << ',';
  os << z.deriv(N - 1) << ')';
  return os;
}


#endif // SURREALS_RVO_H
