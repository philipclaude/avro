// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MDFORDERING_H
#define MDFORDERING_H

#include <petscsys.h>

#include <vector>

#include "LinearAlgebra/SparseLinAlg/ScalarMatrix_CRS.h"

namespace SANS
{
namespace SLA
{

std::vector<PetscInt> computeOrdering_MDF(const SparseMatrix_CRS<PetscScalar>& Cmat);

} //namespace SLA
} //namespace SANS

#endif //MDFORDERING_H
