// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef IDENTITY_H
#define IDENTITY_H

#include "Python/PyDict.h" //This must be first

#include "tools/noncopyable.h"

#include "numpack/AlgebraicEquationSetBase.h"
#include "numpack/sparse/LinearSolverBase.h"
#include "numpack/sparse/sparse_Inverse.h"

namespace numpack 
{
namespace SLA
{


//=============================================================================
struct IdentityParam : noncopyable
{
  static void checkInputs(PyDict d);
  static IdentityParam params;
};


//=============================================================================
template< class Matrix_type >
class Identity : public LinearSolverBase< Matrix_type >
{
public:
  typedef LinearSolverBase< Matrix_type > Base_type;

  typedef typename Base_type::SparseVector_type SparseVector_type;
  typedef typename Base_type::SparseVectorView_type SparseVectorView_type;
  typedef typename NonZeroPatternType<Matrix_type>::type SystemNonZeroPattern;

//-----------------------------------------------------------------------------
  Identity( AlgebraicEquationSetBase<Matrix_type>& f, LinearSystemSolve solve = RegularSolve ) :
    Base_type(solve), f_(f) { init(); }
  explicit Identity( const PyDict& d, AlgebraicEquationSetBase<Matrix_type>& f, LinearSystemSolve solve = RegularSolve ) :
    Base_type(solve), f_(f) { init(); }
  virtual ~Identity() { delete A_; }

//-----------------------------------------------------------------------------
  virtual void factorize() override
  {
    // update the matrix
    *A_ = 0;
    f_.jacobian(*A_);
  }

//-----------------------------------------------------------------------------
  virtual LinearSolveStatus backsolve( const SparseVectorView_type& b, SparseVectorView_type& x ) const override
  {
    //The identity predoncidioner does nothing by design. Mostly useful for unit testing.
    x = b;

    return LinearSolveStatus(true);
  }

protected:
  void init()
  {
    SystemNonZeroPattern nz(f_.matrixSize());
    f_.jacobian(nz);

    A_ = new Matrix_type(nz);
  }

  using Base_type::A_;
  AlgebraicEquationSetBase<Matrix_type>& f_;
};

} //namespace SLA
} //namespace numpack 

#endif //IDENTITY_H
