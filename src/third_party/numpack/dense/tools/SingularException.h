// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef DENSEMATRIX_SINGULAR_EXCEPTION_H
#define DENSEMATRIX_SINGULAR_EXCEPTION_H

#include "tools/SANSException.h"
#include "tools/SANSnumerics.h" // Real

#include <string>

template<int N, class T>
class SurrealS;

namespace numpack 
{

namespace DLA
{
//An exception class representing the attempt to invert a singular matrix
struct SingularMatrixException : public BackTraceException
{
  SingularMatrixException(const Real denom)
  {
    errMessage(denom);
  }

  template<int N>
  SingularMatrixException(const SurrealS<N,Real>& denom)
  {
    errMessage(denom.value());
  }

  void add(const std::string& msg)
  {
    errString += msg;
  }

  virtual ~SingularMatrixException() throw() {}

protected:
  void errMessage(const Real denom);
};

// check if denom is zero...
#define SANS_ASSERT_NONSINGULAR( denom ) \
  if ( unlikely( denom == 0 ) ) \
  {}
    //BOOST_THROW_EXCEPTION( SingularMatrixException( denom ) )

//if ( unlikely( (denom) < std::numeric_limits<Real>::min()*100 && (denom) > -std::numeric_limits<Real>::min()*100 ) )

} //namespace DLA
} //namespace numpack 


#endif //DENSEMATRIX_SINGULAR_EXCEPTION_H
