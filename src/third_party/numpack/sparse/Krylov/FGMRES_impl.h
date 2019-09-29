// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

// No include blocker here. This file should only be included once in a cpp file for explicit instantiation

#if !defined(FGMRES_INSTANTIATE) && !defined(SANS_HEADERCOMPILE_CHECK)
#error "This file should only be included in a cpp file for explicit instantiations"
#endif

#include "FGMRES.h"

#include "numpack/sparse/tools/Array.h"
#include "numpack/sparse/tools/UpperTriangMatrix.h"
#include "numpack/sparse/sparse_Add.h"
#include "numpack/sparse/sparse_Mul.h"

#include "numpack/dense/tools/norm.h"
#include "numpack/dense/tools/dot.h"
#include "numpack/sparse/tools/norm.h"
#include "numpack/sparse/tools/dot.h"

#include "tools/minmax.h"

#include <cmath>
#include <iostream>
#include <sstream>


namespace numpack 
{
namespace SLA
{

//-----------------------------------------------------------------------------
template< class Matrix_type >
FGMRES<Matrix_type>::FGMRES( const PyDict& d, Solver_ptr M ) :
  Base_type(M->systemSolve()),
  params(FGMRESParam::params),
  M_(M),
  tol_(d.get(params.tol)),
  tolFinal_(0),
  nInner_(d.get(params.nInner)),
  nOuter_(d.get(params.nOuter)),
  printConv_(d.get(params.PrintCovergence)),
  itTotal_(0)
{
}

//-----------------------------------------------------------------------------
template< class Matrix_type >
FGMRES<Matrix_type>::~FGMRES() {}

//-----------------------------------------------------------------------------
template< class Matrix_type >
void FGMRES<Matrix_type>::factorize()
{
  M_->factorize();

  // the pointer simply references the matrix stored in the preconditioner
  A_ = &const_cast<Matrix_type&>(M_->A());
}

//-----------------------------------------------------------------------------
template< class Matrix_type >
LinearSolveStatus FGMRES<Matrix_type>::backsolve( const SparseVectorView_type& b, SparseVectorView_type& x ) const
{
  int itOuter = 0, mip, i, j;
  static const Real delta = Real(1.0E-03);
  Real av = 0, y1 = 0, y2 = 0, g1 = 0, g2 = 0, mu = 0, htmp = 0, hmip = 0, hmi = 0;
  Real rho = 0;

  //Allocate temporary vectors
  SparseVector_type r(x.size());
  std::vector<SparseVector_type> v(nInner_+1, x.size()), z(nInner_, x.size()); //Need z for Flexible GMRES
  UpperTriangMatrix h(nInner_+1);
  Array<Real> y(nInner_+1), c(nInner_+1), s(nInner_+1), g(nInner_+1);

  const Matrix_type& A = M_->A();

  const_cast<FGMRES*>(this)->itTotal_ = 0;

  y = 0;
  rho = tol_ + Real(1);

  while (rho > tol_ && itOuter < nOuter_)
  {
    r   = b - A*x;    // Compute the residual vector r
    rho = norm(r,2);  // Compute norm of the residual

    if (rho < tol_) break; //Finish so we do not divide by zero.

    v[0] = r/rho;

    g = 0.;
    h = 0.;

    g[0] = rho;

    int mi = 0;
    if (printConv_)
      std::cout << "FGMRES | Outer : " << itOuter << " | Inner : " << mi << " | Resid : " << rho << std::endl;

    //Always want to do at least one iteration
    rho = rho+1;
    while (rho > tol_ && mi < nInner_)
    {
      mip = mi + 1;
      M_->backsolve(v[mi], z[mi]);  // Flexible GMRES
      v[mip] = A*z[mi];             // Flexible GMRES
      //v[mip] = A*v[mi];           // Original GMRES
      av = norm(v[mip],2);

      //
      // Othogonalize the vectors
      //
      for (j = 0; j <= mi; ++j)
      {
        h(j,mi) = dot(v[mip],v[j]);
        v[mip] -= h(j,mi)*v[j];
      }
      hmip = norm(v[mip],2);

      if ( av + delta * hmip == av)
      {
         for (j = 0; j <= mi; ++j)
         {
           htmp     = dot(v[mip],v[j]);
           h(j,mi) += htmp;
           v[mip]  -= htmp * v[j];
         }
         hmip = norm(v[mip],2);
      }

      v[mip] /= hmip;

      //
      //  Multiply by givens rotation
      //
      for (j = 0; j < mi; ++j)
      {
        y1 = c[j] * h(j,mi) - s[j] * h(j+1,mi);
        y2 = s[j] * h(j,mi) + c[j] * h(j+1,mi);

        h(j,mi)   = y1;
        h(j+1,mi) = y2;
        y[j]      = y1;
        y[j+1]    = y2;
      }

      hmi   = h(mi,mi);
      mu    = sqrt(hmi*hmi + hmip*hmip);
      c[mi] =  hmi / mu;
      s[mi] = -hmip / mu;

      h(mi,mi) = c[mi] * hmi - s[mi] * hmip;

      //
      //  Multiply by givens rotation
      //
      g1 = c[mi] * g[mi] - s[mi] * g[mip];
      g2 = s[mi] * g[mi] + c[mi] * g[mip];

      g[mi]  = g1;
      g[mip] = g2;

      rho = fabs(g1);

      mi++;

      if (printConv_)
        std::cout << "FGMRES | Outer : " << itOuter << " | Inner : " << mi << " | Resid : " << rho << std::endl;

    }// while

    //
    //Solve the upper-triangular system
    //
    y[mi] = g[mi]/hmip; //h(mip,mi);

    for (j = mi-1; j >= 0; --j)
    {
      htmp = 0;
      for (i = j+1; i <= mi; ++i)
        htmp += h(j,i)*y[i];

      y[j] = (g[j] - htmp)/h(j,j);
    }

    //
    // Compute the solution vector
    //
    for (i = 0; i < mi; ++i)
      x += y[i]*z[i];   // Flexible GMRES
//    z = M_->backsolve(A)*v[i]; // Original GMRES
//    x += y[i]*z;      // Original GMRES

    itOuter++;
    itTotal_ += nInner_;
  }

  tolFinal_ = rho;

  if (printConv_)
  {
    //std::stringstream print;
    //print << "\rFGMRES converged to " << rho << " in " << ittot << " iterations." << std::endl << std::endl;
    std::cout << "FGMRES converged to " << rho << " in " << itTotal_ << " iterations." << std::endl << std::endl;
    //CommMatrix.Group.erase_print();
    //CommMatrix.Group.print( print.str() );
  }

  return LinearSolveStatus(true);
}


} //namespace SLA
} //namespace numpack 
