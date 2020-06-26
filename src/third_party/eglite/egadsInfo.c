/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Base Object Functions
 *
 *      Copyright 2011-2019, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#if !defined(WIN32) && !defined(__CYGWIN__)
#define LONG long
#include <execinfo.h>
#else
#ifdef _OCC64
#define LONG long long
#else
#define LONG long
#endif
#endif

#include "egadsTypes.h"
#include "egadsInternals.h"


#define ZERO            1.e-5           /* allow for float-like precision */
#define STRING(a)       #a
#define STR(a)          STRING(a)


  static char *EGADSprop[2] = {STR(EGADSPROP),
                               "\nEGADSprop: Copyright 2011-2019 MIT. All Rights Reserved."};


  extern void EG_initOCC( );
  extern void EG_exactInit( );
  extern int  EG_destroyGeometry( egObject *geom );
  extern int  EG_destroyTopology( egObject *topo );
  extern int  EG_copyGeometry( const egObject *geom, /*@null@*/ double *xform,
                               egObject **copy );
  extern int  EG_copyTopology( const egObject *topo, /*@null@*/ double *xform,
                               egObject **copy );
  extern int  EG_flipGeometry( const egObject *geom, egObject **copy );
  extern int  EG_flipTopology( const egObject *topo, egObject **copy );
  extern int  EG_getTopology( const egObject *topo, egObject **geom,
                              int *ocls, int *type, /*@null@*/ double *limits,
                              int *nobjs, egObject ***objs, int **senses );
  extern int  EG_getGeometry( const egObject *geom, int *oclass, int *mtype,
                              egObject **refGeom, int **ivec, double **rvec );
  extern int  EG_getTolerance( const egObject *topo, double *tol );

  extern int  EG_computeTessMap( egTessel *btess, int outLevel );
  extern int  EG_initTessBody( egObject *object, egObject **tess );
  extern int  EG_getTessEdge( const egObject *tess, int eIndex, int *len,
                              const double **xyz, const double **t );
  extern int  EG_setTessEdge( egObject *tess, int eIndex, int len,
                              const double *xyz, const double *t );
  extern int  EG_getTessFace( const egObject *tess, int fIndex, int *len,
                              const double **xyz, const double **uv,
                              const int **ptype, const int **pindex, int *ntri,
                              const int **tris, const int **tric );
  extern int  EG_setTessFace( egObject *tess, int fIndex, int len,
                              const double *xyz, const double *uv, int ntri,
                              const int *tris );



static void
EG_traceback()
{
#if !defined(WIN32) && !defined(__CYGWIN__)
  int    i;
  void   *array[100];
  size_t size;
  char   **strings;

  size = backtrace(array, 100);
  strings = backtrace_symbols (array, size);
  i = size;
  printf ("\nObtained %d stack frames:\n", i);
  for (i = 0; i < size; i++)
    printf ("%s\n", strings[i]);
  free (strings);
  printf("\n");
#endif
}


int
EG_getInfo(const egObject *object, int *oclass, int *mtype, egObject **top,
           egObject **prev, egObject **next)
{
  if (object == NULL)               return EGADS_NULLOBJ;
  if (object->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (object->oclass == EMPTY)      return EGADS_EMPTY;

  *oclass = object->oclass;
  *mtype  = object->mtype;
  *top    = object->topObj;
  *prev   = object->prev;
  *next   = object->next;

  return EGADS_SUCCESS;
}
