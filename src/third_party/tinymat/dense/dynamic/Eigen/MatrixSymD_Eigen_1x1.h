// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "../Eigen.h"

#include "../MatrixSymD.h"
#include "../VectorD.h"

#include "../MatrixD_Det.h"
#include "../MatrixD_Trace.h"

#include "tinymat/types/SurrealD.h"

#include "tools/SANSnumerics.h"

#include <cassert>

namespace tinymat
{
namespace DLA
{

template<class T>
void
EigenSystem_1x1(const MatrixSymD<T>& A, VectorD<T>& L, MatrixD<T>& E )
{
  L[0] = A(0,0);
  E(0,0) = 1;
}

}
}
