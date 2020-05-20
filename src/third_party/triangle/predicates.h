// avro: Adaptive Voronoi Remesher
// Copyright 2017-2018, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef AVRO_THIRD_PARTY_PREDICATES_H_
#define AVRO_THIRD_PARTY_PREDICATES_H_

// provided in triangle
extern "C" void exactinit();
extern "C" double orient2d(const double*,const double*,const double*);

// provided in tetgen/predicates.cxx
double insphere(double *pa, double *pb, double *pc, double *pd, double *pe);
double incircle(double *pa, double *pb, double *pc, double *pd);

#endif
