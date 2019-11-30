// avro: Adaptive Voronoi Remesher
// Copyright 2017-2018, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef AVRO_THIRD_PARTY_TETGEN_PREDICATES_H_
#define AVRO_THIRD_PARTY_TETGEN_PREDICATES_H_
REAL orient3d( REAL* pa, REAL* pb, REAL* pc, REAL* pd );
REAL orient4d(REAL* pa, REAL* pb, REAL* pc, REAL* pd, REAL* pe,
              REAL aheight, REAL bheight, REAL cheight, REAL dheight,
              REAL eheight);
void exactinit(int verbose, int noexact, int nofilter, REAL maxx, REAL maxy, REAL maxz);
#endif
