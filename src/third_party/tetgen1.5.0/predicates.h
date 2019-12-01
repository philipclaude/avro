#ifndef LUNA_THIRD_PARTY_TETGEN_PREDICATES_H_
#define LUNA_THIRD_PARTY_TETGEN_PREDICATES_H_
REAL orient3d(const REAL* pa, const REAL* pb, const REAL* pc, const REAL* pd);
void exactinit(int verbose, int noexact, int nofilter, REAL maxx, REAL maxy, REAL maxz);
#endif
