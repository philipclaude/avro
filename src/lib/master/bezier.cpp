//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "master/basis.h"
#include "master/simplex.h"

namespace avro
{

void
Bezier<Simplex>::eval(const ReferenceElement<Simplex>& ref , const double* x , double* phi )
{
  printf("eval in bezier!\n");
}

void
Bezier<Simplex>::grad( const ReferenceElement<Simplex>& ref , const double* x , double* gphi )
{
  printf("grad in bezier!\n");
}

void
Bezier<Simplex>::hess(const ReferenceElement<Simplex>& ref , const double* x , double* hphi )
{
  printf("hess in bezier!\n");
}

} // avro
