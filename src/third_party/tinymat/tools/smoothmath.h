// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef SMOOTHMATH_H
#define SMOOTHMATH_H

#include <type_traits>
#include <cmath> //log, exp, sin, tan

#include "SANSnumerics.h"

#include <boost/type_traits/is_arithmetic.hpp>

// Smooth abs that goes through 0 at x = 0
// The curve is shifted away from x otherwise, i.e. smoothabs0(x) < x
// C1 continuous
template<typename T>
inline typename std::enable_if< boost::is_arithmetic<T>::value, T>::type
smoothabs0(const T x, const Real eps)
{
  Real signx = x > 0 ? 1 : -1; // Sign function

  return x*x / (signx*x + eps);
}

// Smooth abs that goes through eps/2 at x = 0
// For abs(x) > eps the function is equal to x, i.e. smoothabsP(x) == x
template<typename T>
inline typename std::enable_if< boost::is_arithmetic<T>::value, T>::type
smoothabsP(const T x, const Real eps)
{
  Real signx = x > 0 ? 1 : -1; // Sign function

  if (signx*x < eps)
    return 0.5*(eps + x*x/eps);

  return signx*x;
}


//-------------------------------------------------------------------------------------------------------

// Larger values of alpha approach a sharp max
template<typename T>
inline typename std::enable_if< boost::is_arithmetic<T>::value, T>::type
smoothmax(const T x, const T y, const Real alpha)
{
  // The if statements below is a more precise version of
  // T m = max(x,y);
  // return m + log( exp(alpha*(x-m)) + exp(alpha*(y-m)) )/alpha;

  if ( x > y )
    return x + log1p( exp(alpha*(y-x)) )/alpha;
  else
    return y + log1p( exp(alpha*(x-y)) )/alpha;
}

// Larger values of alpha approach a sharp min
template<typename T>
inline typename std::enable_if< boost::is_arithmetic<T>::value, T>::type
smoothmin(const T x, const T y, const Real alpha)
{
  // The if statements below is a more precise version of
  // T m = min(x,y);
  // return m - log( exp(-alpha*(x-m)) + exp(-alpha*(y-m)) )/alpha;

  if ( x < y )
    return x - log1p( exp(-alpha*(y-x)) )/alpha;
  else
    return y - log1p( exp(-alpha*(x-y)) )/alpha;
}



// Smooth max with C2 continuity that matches at x = +- eps/2
template<typename T>
inline typename std::enable_if< boost::is_arithmetic<T>::value, T>::type
smoothmaxC2(const T x, const T y, const Real eps = 1e-2 )
{
  const T m = x > y ? y : x;
  const T M = x > y ? x : y;
  const T d = M-m; // max(x,y) = y + max(0,x-y);

  // from solving linear system with 6 boundary conditions
  return ( d > eps / 2 ) ? M : m + (pow(eps+2*d,3)*(3*eps - 2*d))/(32*pow(eps,3));
}

// Smooth max with C2 continuity that matches at x = +- eps/2
template<typename T>
inline typename std::enable_if< boost::is_arithmetic<T>::value, T>::type
smoothminC2(const T x, const T y, const Real eps = 1e-2 )
{
  // min(a,b) = -max(-a,-b)
  return - smoothmaxC2(-x,-y,eps);
}

// smoooth max with C1 continuity that matches at x = +- eps/2
template<typename T>
inline typename std::enable_if< boost::is_arithmetic<T>::value, T>::type
smoothmaxC1(const T x, const T y, const Real eps = 1e-2 )
{
  const T m = x > y ? y : x;
  const T M = x > y ? x : y;
  const T d = M-m; // max(x,y) = y + max(0,x-y);

  // from solving linear system with 4 boundary conditions
  return ( d > eps / 2 ) ? M : m + pow(eps + 2*d,2)/(8*eps);
}

// smoooth max with C1 continuity that matches at x = +- eps/2
template<typename T>
inline typename std::enable_if< boost::is_arithmetic<T>::value, T>::type
smoothminC1(const T x, const T y, const Real eps = 1e-2 )
{
  // min(a,b) = -max(-a,-b)
  return - smoothmaxC1(-x,-y,eps);
}

//-------------------------------------------------------------------------------------------------------

// C1 continuous smooth ramp function that goes through 0 at x = 0
// If x < 0, the function returns 0, else it returns an approximation to x
template<typename T>
inline typename std::enable_if< boost::is_arithmetic<T>::value, T>::type
smoothRamp0(const T x, const Real eps)
{
  if (x >= 0)
    return x*x / (x + eps);
  else
    return 0;
}

//-------------------------------------------------------------------------------------------------------

//Smooth activation functions - a function that smoothly ramps from 0 to 1 for inputs from 0 to 1
// y = 0    for x < 0
// y = f(x) for 0 < x < 1
// y = 1    for x > 1

//C1 continuous
template<typename T>
inline typename std::enable_if< boost::is_arithmetic<T>::value, T>::type
smoothActivation_sine(const T x)
{
  if (x < 0.0)
    return 0.0;
  else if (x > 1.0)
    return 1.0;
  else
    return 0.5*(1.0 + sin(PI*(x-0.5)));
}

//C1 continuous
template<typename T>
inline typename std::enable_if< boost::is_arithmetic<T>::value, T>::type
smoothActivation_cubic(const T x)
{
  if (x < 0.0)
    return 0.0;
  else if (x > 1.0)
    return 1.0;
  else
    return 3.0*(x*x) - 2.0*(x*x*x);
}

//C-infinity continuous
template<typename T>
inline typename std::enable_if< boost::is_arithmetic<T>::value, T>::type
smoothActivation_tanh(const T x, const Real k = 4.0)
{
  return 0.5*(1.0 + tanh(k*(x - 0.5)));
}

// Special smooth activation function for shock capturing
template<typename T>
inline typename std::enable_if< boost::is_arithmetic<T>::value, T>::type
smoothActivation_sine(const T x, const Real xL, const Real xH)
{
  if (x < xL)
    return 0.0;
  else if (x > xH)
    return xH;
  else
    return 0.5*xH*(sin(PI*((x - xL)/(xH - xL) - 0.5)) + 1.0);
}

// Another special smooth activation function for shock capturing
template<typename T>
inline typename std::enable_if< boost::is_arithmetic<T>::value, T>::type
smoothActivation_exp(const T x, const Real alpha, const Real eps = 0)
{
  // Compute the the offset 'b' so smoothActivation_exp(0) == eps, if eps != 0
  Real b = 0;
  if (eps != 0) b = -std::log(std::expm1(alpha*eps))/alpha;

  // Avoid overflow/underflow in exp function
  if (-alpha*(x - b) > 100)
    return 0;
  else if (-alpha*(1 - x) > 100)
    return 1;
  else
    return x - b + ( log1p( exp(-alpha*(x - b)) ) - log1p( exp(-alpha*(1 - (x - b))) ) )/alpha;
}


#endif //SMOOTHMATH_H
