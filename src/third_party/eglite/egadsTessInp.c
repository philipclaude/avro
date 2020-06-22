/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Tessellation Input Functions
 *
 *      Copyright 2011-2019, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "egadsTypes.h"
#include "egadsInternals.h"
#include "egadsTris.h"
#include "emp.h"

#include "regQuads.h"


  typedef struct {
    void     *mutex;              /* the mutex or NULL for single thread */
    long     master;              /* master thread ID */
    int      end;                 /* end of loop */
    int      index;               /* current loop index */
    egTessel *ntess;              /* tessellation structure */
    /*@dependent@*/
    bodyQuad *bodydata;           /* the quad storage */
  } EMPquad;


  typedef struct {
    int vert2;                    /* the second vert */
    int nvert;                    /* the new vert */
    int next;                     /* the next in the list for the first vert */
  } midside;


#define MINARC
#define REGULAR
#define NOTFILLED       -1
#define EPS10           1.e-10


  extern void EG_cleanupTess( egTessel *btess );
  extern void EG_cleanupTessMaps( egTessel *btess );
  extern void EG_makeConnect( int k1, int k2, int *tri, int *kedge, int *ntable,
                              connect *etable, int face );
  extern int  EG_fillArea( int ncontours, const int *cntr,
                           const double *vertices, int *triangles, int *n_fig8,
                           int pass, fillArea *fa );

  extern int  EG_getTopology( const egObject *topo, egObject **geom, int *oclas,
                              int *type, /*@null@*/ double *limits, int *nChild,
                              egObject ***children, int **senses );
  extern int  EG_getBodyTopos( const egObject *body, /*@null@*/ egObject *src,
                               int oclass, int *ntopo,
                               /*@null@*/ egObject ***topos );
  extern int  EG_indexBodyTopo( const egObject *body, const egObject *src );
  extern int  EG_evaluate( const egObject *geom, const double *param,
                           double *result );
  extern int  EG_getRange( const egObject *geom, double *range, int *pflag );
  extern int  EG_invEvaluateGuess( const egObject *geom, double *xyz,
                                   double *param, double *result );
  extern int  EG_getEdgeUV( const egObject *face, const egObject *edge,
                            int sense, double t, double *result );
  extern int  EG_getEdgeUVeval( const egObject *face, const egObject *topo,
                                int sense, double t, double *result);
  extern int  EG_getTessEdge( const egObject *tess, int indx, int *len,
                              const double **xyz, const double **t );
  extern int  EG_getTessFace( const egObject *tess, int indx, int *len,
                              const double **xyz, const double **uv,
                              const int **ptype, const int **pindex,
                              int *ntri, const int **tris, const int **tric );
  extern int  EG_attributeRet( const egObject *obj, const char *name,
                               int *type, int *len, /*@null@*/ const int **ints,
                               /*@null@*/ const double **reals,
                               /*@null@*/ const char   **str );
#ifndef LITE
  extern int  EG_attributeDel( egObject *obj, /*@null@*/ const char *name );
  extern int  EG_attributeAdd( egObject *obj, const char *name, int type,
                               int len, /*@null@*/ const int    *ints,
                                        /*@null@*/ const double *reals,
                                        /*@null@*/ const char   *str );
#endif


int
EG_openTessBody(egObject *tess)
{
  egTessel *btess;
  egObject *obj;

  if  (tess == NULL)                 return EGADS_NULLOBJ;
  if  (tess->magicnumber != MAGIC)   return EGADS_NOTOBJ;
  if  (tess->oclass != TESSELLATION) return EGADS_NOTTESS;
  btess = (egTessel *) tess->blind;
  if (btess == NULL) {
    printf(" EGADS Error: NULL Blind Object (EG_openTessBody)!\n");
    return EGADS_NOTFOUND;
  }
  obj = btess->src;
  if (obj == NULL) {
    printf(" EGADS Error: NULL Source Object (EG_openTessBody)!\n");
    return EGADS_NULLOBJ;
  }
  if (obj->magicnumber != MAGIC) {
    printf(" EGADS Error: Source Not an Object (EG_openTessBody)!\n");
    return EGADS_NOTOBJ;
  }
  if (obj->oclass != BODY) {
    printf(" EGADS Error: Source Not Body (EG_openTessBody)!\n");
    return EGADS_NOTBODY;
  }
  if (btess->done != 1) return EGADS_TESSTATE;

  /* set open state and clean up any local/global mappings */
  btess->done = 0;
  EG_cleanupTessMaps(btess);

  return EGADS_SUCCESS;
}


int
EG_initTessBody(egObject *object, egObject **tess)
{
  int      i, j, stat, outLevel, nedge, nface, oclass, mtype, nnode, *senses;
  double   limits[2];
  egTessel *btess;
  egObject *ttess, *context, *geom, **edges, **nodes;

  *tess = NULL;
  if (object == NULL)               return EGADS_NULLOBJ;
  if (object->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (object->oclass != BODY)       return EGADS_NOTBODY;
  outLevel = EG_outLevel(object);
  context  = EG_context(object);
  if (context == NULL)              return EGADS_NULLOBJ;

  stat = EG_getBodyTopos(object, NULL, EDGE, &nedge, &edges);
  if (stat != EGADS_SUCCESS) return stat;
  stat = EG_getBodyTopos(object, NULL, FACE, &nface, NULL);
  if (stat  != EGADS_SUCCESS) return stat;

  btess = (egTessel *) EG_alloc(sizeof(egTessel));
  if (btess == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Blind Malloc (EG_initTessBody)!\n");
    EG_free(edges);
    return EGADS_MALLOC;
  }
  btess->src       = object;
  btess->xyzs      = NULL;
  btess->tess1d    = NULL;
  btess->tess2d    = NULL;
  btess->globals   = NULL;
  btess->nGlobal   = 0;
  btess->nEdge     = nedge;
  btess->nFace     = nface;
  btess->nu        = 0;
  btess->nv        = 0;
  btess->done      = 0;
  btess->params[0] = 0.0;
  btess->params[1] = 0.0;
  btess->params[2] = 0.0;
  for (i = 0; i < MTESSPARAM; i++) btess->tparam[i] = 0.0;

  btess->tess1d = (egTess1D *) EG_alloc(nedge*sizeof(egTess1D));
  if (btess->tess1d == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Alloc %d Edges (EG_initTessBody)!\n", nedge);
    EG_free(edges);
    EG_free(btess);
    return EGADS_MALLOC;
  }
  for (j = 0; j < nedge; j++) {
    btess->tess1d[j].obj            = edges[j];
    btess->tess1d[j].faces[0].index = 0;
    btess->tess1d[j].faces[0].nface = 0;
    btess->tess1d[j].faces[0].faces = NULL;
    btess->tess1d[j].faces[0].tric  = NULL;
    btess->tess1d[j].faces[1].index = 0;
    btess->tess1d[j].faces[1].nface = 0;
    btess->tess1d[j].faces[1].faces = NULL;
    btess->tess1d[j].faces[1].tric  = NULL;
    btess->tess1d[j].nodes[0]       = 0;
    btess->tess1d[j].nodes[1]       = 0;
    btess->tess1d[j].xyz            = NULL;
    btess->tess1d[j].t              = NULL;
    btess->tess1d[j].global         = NULL;
    btess->tess1d[j].npts           = 0;
  }
  EG_free(edges);
  for (j = 0; j < nedge; j++) {
    stat = EG_getTopology(btess->tess1d[j].obj, &geom, &oclass, &mtype, limits,
                          &nnode, &nodes, &senses);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: EG_getTopology = %d (EG_initTessBody)!\n", stat);
      EG_cleanupTess(btess);
      EG_free(btess);
      return EGADS_MALLOC;
    }
    stat = EG_indexBodyTopo(object, nodes[0]);
    if (stat <= EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: EG_indexBodyTopo0 = %d (EG_initTessBody)!\n",
               stat);
      EG_cleanupTess(btess);
      EG_free(btess);
      return EGADS_MALLOC;
    }
    btess->tess1d[j].nodes[0] = btess->tess1d[j].nodes[1] =  stat;
    if (mtype == DEGENERATE)    btess->tess1d[j].nodes[1] = -stat;
    if (nnode == 1) continue;
    stat = EG_indexBodyTopo(object, nodes[1]);
    if (stat < EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: EG_indexBodyTopo1 = %d (EG_initTessBody)!\n",
               stat);
      EG_cleanupTess(btess);
      EG_free(btess);
      return EGADS_MALLOC;
    }
    btess->tess1d[j].nodes[1] = stat;
  }

  if (nface != 0) {
    btess->tess2d = (egTess2D *) EG_alloc(2*nface*sizeof(egTess2D));
    if (btess->tess2d == NULL) {
      printf(" EGADS Error: Alloc %d Faces (EG_initTessBody)!\n", nface);
      EG_cleanupTess(btess);
      EG_free(btess);
      return EGADS_MALLOC;
    }
    for (j = 0; j < 2*nface; j++) {
      btess->tess2d[j].mKnots = NULL;
      btess->tess2d[j].xyz    = NULL;
      btess->tess2d[j].uv     = NULL;
      btess->tess2d[j].global = NULL;
      btess->tess2d[j].ptype  = NULL;
      btess->tess2d[j].pindex = NULL;
      btess->tess2d[j].bary   = NULL;
      btess->tess2d[j].frame  = NULL;
      btess->tess2d[j].frlps  = NULL;
      btess->tess2d[j].tris   = NULL;
      btess->tess2d[j].tric   = NULL;
      btess->tess2d[j].patch  = NULL;
      btess->tess2d[j].npts   = 0;
      btess->tess2d[j].nframe = 0;
      btess->tess2d[j].nfrlps = 0;
      btess->tess2d[j].ntris  = 0;
      btess->tess2d[j].npatch = 0;
      btess->tess2d[j].tfi    = 0;
    }
  }

  stat = EG_makeObject(context, &ttess);
  if (stat != EGADS_SUCCESS) {
    EG_cleanupTess(btess);
    EG_free(btess);
    return stat;
  }
  ttess->oclass = TESSELLATION;
  ttess->blind  = btess;
  EG_referenceObject(ttess,  context);
  EG_referenceTopObj(object, ttess);
  *tess = ttess;

  return EGADS_SUCCESS;
}


int
EG_computeTessMap(egTessel *btess, int outLevel)
{
  int i, j, k, n, npts, pt, pi, nNode, *inode;

  if (btess->nGlobal !=    0) return EGADS_EXISTS;
  if (btess->globals != NULL) return EGADS_EXISTS;

  /* special case of a degenerate WIREBODY -- NodeBody */
  if ((btess->nEdge == 1) && (btess->tess1d[0].obj->mtype == DEGENERATE) &&
      (btess->tess1d[0].nodes[0] == 1)) {
    btess->xyzs = (double *) EG_alloc(3*sizeof(double));
    if (btess->xyzs == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Allocation of 1 Nodes (EG_computeTessMap)!\n");
      return EGADS_MALLOC;
    }
    btess->xyzs[0] = btess->tess1d[0].xyz[0];
    btess->xyzs[1] = btess->tess1d[0].xyz[1];
    btess->xyzs[2] = btess->tess1d[0].xyz[2];

    btess->tess1d[0].global = (int *) EG_alloc(2*sizeof(int));
    if (btess->tess1d[0].global == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Allocation of 1 Global (EG_computeTessMap)!\n");
      return EGADS_MALLOC;
    }
    btess->tess1d[0].global[0] = 1;
    btess->tess1d[0].global[1] = 1;

    btess->globals = (int *) EG_alloc(2*sizeof(int));
    if (btess->globals == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Allocation of 1 Nodes (EG_computeTessMap)!\n");
      return EGADS_MALLOC;
    }
    btess->globals[0] = 0;
    btess->globals[1] = 1;
    btess->nGlobal    = 1;
    return EGADS_SUCCESS;
  }

  /* get Node and Edge sizes */
  for (nNode = npts = i = 0; i < btess->nEdge; i++) {
    if (btess->tess1d[i].obj        ==       NULL) continue;
    if (btess->tess1d[i].obj->mtype == DEGENERATE) continue;
    if (btess->tess1d[i].nodes[1]   <           0) continue;
    npts += btess->tess1d[i].npts;
    if (nNode < btess->tess1d[i].nodes[0]) nNode = btess->tess1d[i].nodes[0];
    if (nNode < btess->tess1d[i].nodes[1]) nNode = btess->tess1d[i].nodes[1];
  }
  inode = (int *) EG_alloc(nNode*sizeof(int));
  if (inode == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Allocation of %d Nodes (EG_computeTessMap)!\n",
             nNode);
    return EGADS_MALLOC;
  }
  btess->xyzs = (double *) EG_alloc(3*nNode*sizeof(double));
  if (btess->xyzs == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Allocation Of %d Nodes (EG_computeTessMap)!\n",
             nNode);
    EG_free(inode);
    return EGADS_MALLOC;
  }

  for (i = 0; i < nNode; i++) {
    inode[i]           = 0;
    btess->xyzs[3*i  ] = 0.0;
    btess->xyzs[3*i+1] = 0.0;
    btess->xyzs[3*i+2] = 0.0;
  }
  for (i = 0; i < btess->nEdge; i++) {
    if (btess->tess1d[i].obj        ==       NULL) continue;
    if (btess->tess1d[i].obj->mtype == DEGENERATE) continue;
    if (btess->tess1d[i].nodes[1]   <           0) continue;
    j = btess->tess1d[i].nodes[0]-1;
    btess->xyzs[3*j  ] = btess->tess1d[i].xyz[0];
    btess->xyzs[3*j+1] = btess->tess1d[i].xyz[1];
    btess->xyzs[3*j+2] = btess->tess1d[i].xyz[2];
    inode[j]++;
    j = btess->tess1d[i].nodes[1]-1;
    k = btess->tess1d[i].npts;
    btess->xyzs[3*j  ] = btess->tess1d[i].xyz[3*k-3];
    btess->xyzs[3*j+1] = btess->tess1d[i].xyz[3*k-2];
    btess->xyzs[3*j+2] = btess->tess1d[i].xyz[3*k-1];
    inode[j]++;
    btess->tess1d[i].global = (int *) EG_alloc(k*sizeof(int));
    if (btess->tess1d[i].global == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: %d Allocation of %d Global (EG_computeTessMap)!\n",
               i+1, k);
      EG_cleanupTessMaps(btess);
      EG_free(inode);
      return EGADS_MALLOC;
    }
    for (j = 0; j < k; j++) btess->tess1d[i].global[j] = 0;
  }

  if (btess->nFace == 0) {

    /* deal with wirebodies */

    for (i = 0; i < nNode; i++)
      if (inode[i] != 0) npts -= inode[i]-1;

    btess->globals = (int *) EG_alloc(2*npts*sizeof(int));
    if (btess->globals == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Allocation of %d Globals (EG_computeTessMap)!\n",
               npts);
      EG_cleanupTessMaps(btess);
      EG_free(inode);
      return EGADS_MALLOC;
    }

    for (k = i = 0; i < btess->nEdge; i++) {
      if (btess->tess1d[i].obj        ==       NULL) continue;
      if (btess->tess1d[i].obj->mtype == DEGENERATE) continue;
      if (btess->tess1d[i].nodes[1]   <           0) continue;
      n = btess->tess1d[i].nodes[0] - 1;
      if (inode[n] > 0) {
        btess->globals[2*k  ]      =   0;
        btess->globals[2*k+1]      = n+1;
        k++;
        btess->tess1d[i].global[0] =   k;
        inode[n] = -k;
      } else {
        btess->tess1d[i].global[0] = -inode[n];
      }

      for (j = 1; j < btess->tess1d[i].npts-1; j++) {
        btess->globals[2*k  ]      = j+1;
        btess->globals[2*k+1]      = i+1;
        k++;
        btess->tess1d[i].global[j] =   k;
      }

      n = btess->tess1d[i].nodes[1] - 1;
      j = btess->tess1d[i].npts     - 1;
      if (inode[n] > 0) {
        btess->globals[2*k  ]      =   0;
        btess->globals[2*k+1]      = n+1;
        k++;
        btess->tess1d[i].global[j] =   k;
        inode[n] = -k;
      } else {
        btess->tess1d[i].global[j] =  -inode[n];
      }

    }

    EG_free(inode);
    btess->nGlobal = npts;
    return EGADS_SUCCESS;
  }

  for (k = i = 0; i < btess->nFace; i++) {
    n = btess->tess2d[i].npts;
    if (n == 0) continue;
    btess->tess2d[i].global = (int *) EG_alloc(n*sizeof(int));
    if (btess->tess2d[i].global == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: %d Allocation of %d Globals (EG_computeTessMap)!\n",
               i+1, n);
      EG_cleanupTessMaps(btess);
      EG_free(inode);
      return EGADS_MALLOC;
    }
    for (j = 0; j < n; j++) btess->tess2d[i].global[j] = 0;
    k += n;
  }
  btess->globals = (int *) EG_alloc(2*k*sizeof(int));
  if (btess->globals == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Allocation of %d Globals (EG_computeTessMap)!\n",
             k);
    EG_cleanupTessMaps(btess);
    EG_free(inode);
    return EGADS_MALLOC;
  }

  for (i = 0; i < nNode; i++) inode[i] = 0;
  for (k = i = 0; i < btess->nFace; i++) {
    n = btess->tess2d[i].npts;
    for (j = 0; j < n; j++) {
      pt = btess->tess2d[i].ptype[j];
      pi = btess->tess2d[i].pindex[j];
      if (pt == 0) {
        if (inode[pi-1] == 0) {
          btess->globals[2*k  ]            = pt;
          btess->globals[2*k+1]            = pi;
          k++;
          btess->tess2d[i].global[j]       = k;
          inode[pi-1]                      = k;
        } else {
          btess->tess2d[i].global[j]       = inode[pi-1];
        }
      } else if (pt > 0) {
        if (btess->tess1d[pi-1].global[pt-1] == 0) {
          btess->globals[2*k  ]            = pt;
          btess->globals[2*k+1]            = pi;
          k++;
          btess->tess2d[i].global[j]       = k;
          btess->tess1d[pi-1].global[pt-1] = k;
        } else {
          btess->tess2d[i].global[j]       = btess->tess1d[pi-1].global[pt-1];
        }
      } else {
        btess->globals[2*k  ]              = -j-1;
        btess->globals[2*k+1]              =  i+1;
        k++;
        btess->tess2d[i].global[j]         =  k;
      }
    }
  }

  /* patch up beginning and end of Edges */
  for (i = 0; i < btess->nEdge; i++) {
    if (btess->tess1d[i].obj        ==       NULL) continue;
    if (btess->tess1d[i].obj->mtype == DEGENERATE) continue;
    if (btess->tess1d[i].nodes[1]   <           0) continue;
    n = btess->tess1d[i].nodes[0] - 1;
    if (inode[n] == 0) {
      btess->globals[2*k  ]      = 0;
      btess->globals[2*k+1]      = n+1;
      k++;
      inode[n]                   = k;
    } else {
      btess->tess1d[i].global[0] = inode[n];
    }
    n = btess->tess1d[i].nodes[1] - 1;
    j = btess->tess1d[i].npts     - 1;
    if (inode[n] == 0) {
      btess->globals[2*k  ]      = 0;
      btess->globals[2*k+1]      = n+1;
      k++;
      inode[n]                   = k;
    } else {
      btess->tess1d[i].global[j] = inode[n];
    }
  }

  EG_free(inode);
  btess->nGlobal = k;

  return EGADS_SUCCESS;
}


int
EG_statusTessBody(egObject *tess, egObject **body, int *state, int *npts)
{
  int          i, j, k, stat, outLevel, atype, alen;
  egTessel     *btess;
  egObject     *obj;
  const int    *ints;
  const double *reals;
  const char   *str;

  *body  = NULL;
  *state = *npts = 0;
  if (tess == NULL)                 return EGADS_NULLOBJ;
  if (tess->magicnumber != MAGIC)   return EGADS_NOTOBJ;
  if (tess->oclass != TESSELLATION) return EGADS_NOTTESS;
  outLevel = EG_outLevel(tess);

  btess = (egTessel *) tess->blind;
  if (btess == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Blind Object (EG_statusTessBody)!\n");
    return EGADS_NOTFOUND;
  }
  obj = btess->src;
  if (obj == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Source Object (EG_statusTessBody)!\n");
    return EGADS_NULLOBJ;
  }
  if (obj->magicnumber != MAGIC) {
    if (outLevel > 0)
      printf(" EGADS Error: Source Not an Object (EG_statusTessBody)!\n");
    return EGADS_NOTOBJ;
  }
  if (obj->oclass != BODY) {
    if (outLevel > 0)
      printf(" EGADS Error: Source Not Body (EG_statusTessBody)!\n");
    return EGADS_NOTBODY;
  }
  if (btess->tess1d == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: No Edge Tessellations (EG_statusTessBody)!\n");
    return EGADS_NODATA;
  }
  if ((btess->tess2d == NULL) && (btess->nFace != 0)) {
    if (outLevel > 0)
      printf(" EGADS Error: No Face Tessellations (EG_statusTessBody)!\n");
    return EGADS_NODATA;
  }
  *body  = obj;
  *state = btess->done;
  *npts  = btess->nGlobal;

  if (btess->done == 0) {

    /* are we done? */
    for (i = 0; i < btess->nEdge; i++) {
      if (btess->tess1d[i].nodes[0] == -btess->tess1d[i].nodes[1]) continue;
      if (btess->tess1d[i].npts == 0) return EGADS_OUTSIDE;
    }
    for (i = 0; i < btess->nFace; i++)
      if (btess->tess2d[i].npts == 0) return EGADS_OUTSIDE;

    stat = EG_attributeRet(tess, ".mixed", &atype, &alen, &ints, &reals, &str);
    if (stat == EGADS_SUCCESS) {
      if ((alen != btess->nFace) || (atype != ATTRINT)) {
#ifndef LITE
        stat = EG_attributeDel(tess, ".mixed");
        if ((stat != EGADS_SUCCESS) && (outLevel > 0))
          printf(" EGADS Error: Deleting Attribute %d (EG_statusTessBody)!\n",
                 stat);
#endif
      } else {
        for (j = i = 0; i < btess->nFace; i++) {
          if (2*ints[i] == btess->tess2d[i].ntris) j++;
          if (2*ints[i] >  btess->tess2d[i].ntris) {
            printf(" EGADS Error: %d 2*nQuads (%d) > nTris (%d) (EG_statusTessBody)!\n",
                   i+1, ints[i], btess->tess2d[i].ntris);
            j = -1;
            break;
          }
          /* check pairs */
          for (k = btess->tess2d[i].ntris-2*ints[i]; k < btess->tess2d[i].ntris;
               k += 2)
            if ((btess->tess2d[i].tris[3*k  ] != btess->tess2d[i].tris[3*k+3]) ||
                (btess->tess2d[i].tris[3*k+2] != btess->tess2d[i].tris[3*k+4])) {
              printf(" EGADS Error: %d %d Bad Pair %d %d %d - %d %d %d!\n", i+1,
                     k+1, btess->tess2d[i].tris[3*k], btess->tess2d[i].tris[3*k+1],
                     btess->tess2d[i].tris[3*k+2], btess->tess2d[i].tris[3*k+3],
                     btess->tess2d[i].tris[3*k+4], btess->tess2d[i].tris[3*k+5]);
              j = -1;
              break;
            }
          if (j == -1) break;
        }
        if (j == -1) {
#ifndef LITE
          stat = EG_attributeDel(tess, ".mixed");
          if ((stat != EGADS_SUCCESS) && (outLevel > 0))
            printf(" EGADS Error: Deleting Attribute %d (EG_statusTessBody)!\n",
                   stat);
#endif
        } else if (j == btess->nFace) {
#ifndef LITE
          stat = EG_attributeAdd(tess, ".tessType", ATTRSTRING, 4, NULL, NULL,
                                 "Quad");
          if (stat != EGADS_SUCCESS)
            if (outLevel > 0)
              printf(" EGADS Warning: EG_attributeAdd = %d (EG_statusTessBody)!\n",
                     stat);
#endif
        } else {
#ifndef LITE
          stat = EG_attributeAdd(tess, ".tessType", ATTRSTRING, 5, NULL, NULL,
                                 "Mixed");
          if (stat != EGADS_SUCCESS)
            if (outLevel > 0)
              printf(" EGADS Warning: EG_attributeAdd = %d (EG_statusTessBody)!\n",
                     stat);
#endif
        }
      }
    }

    *state = btess->done = 1;
  }

  if (btess->globals != NULL) return EGADS_SUCCESS;

  /* compute the mappings and return the number of global vertices */
  stat = EG_computeTessMap(btess, outLevel);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_computeTessMap = %d (EG_statusTessBody)!\n",
              stat);
    return stat;
  }
  *npts = btess->nGlobal;

  return EGADS_SUCCESS;
}


int
EG_setTessEdge(const egObject *tess, int index, int len, const double *xyz,
               const double *t)
{
  int      i, j, k, n, stat, outLevel, oclass, mtype, nnode, *senses;
  double   xyz0[3], xyz1[3], trange[2], *xyzs, *ts;
  egTessel *btess;
  egObject *obj, **nodes, *geom, **objs;

  if  (tess == NULL)                 return EGADS_NULLOBJ;
  if  (tess->magicnumber != MAGIC)   return EGADS_NOTOBJ;
  if  (tess->oclass != TESSELLATION) return EGADS_NOTTESS;
  if  (len <= 1)                     return EGADS_NODATA;
  if ((xyz == NULL) || (t == NULL))  return EGADS_NODATA;
  outLevel = EG_outLevel(tess);

  btess = (egTessel *) tess->blind;
  if (btess == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Blind Object (EG_setTessEdge)!\n");
    return EGADS_NOTFOUND;
  }
  obj = btess->src;
  if (obj == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Source Object (EG_setTessEdge)!\n");
    return EGADS_NULLOBJ;
  }
  if (obj->magicnumber != MAGIC) {
    if (outLevel > 0)
      printf(" EGADS Error: Source Not an Object (EG_setTessEdge)!\n");
    return EGADS_NOTOBJ;
  }
  if (obj->oclass != BODY) {
    if (outLevel > 0)
      printf(" EGADS Error: Source Not Body (EG_setTessEdge)!\n");
    return EGADS_NOTBODY;
  }
  if (btess->tess1d == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: No Edge Tessellations (EG_setTessEdge)!\n");
    return EGADS_NODATA;
  }
  if (btess->done == 1) {
    if (outLevel > 0)
      printf(" EGADS Error: Complete Tessellation (EG_setTessEdge)!\n");
    return EGADS_EXISTS;
  }
  if ((index < 1) || (index > btess->nEdge)) {
    if (outLevel > 0)
      printf(" EGADS Error: Index = %d [1-%d] (EG_setTessEdge)!\n",
             index, btess->nEdge);
    return EGADS_INDEXERR;
  }

  for (i = 1; i < len; i++) {
    if (t[i] > t[i-1]) continue;
    printf(" EGADS Error: ts not in order %d %lf  %d %lf (EG_setTessEdge)!\n",
           i-1, t[i-1], i, t[i]);
    return EGADS_RANGERR;
  }

  /* are any of our Faces already set? */
  stat = EG_getBodyTopos(obj, btess->tess1d[index-1].obj, FACE, &n, &objs);
  if (stat  != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: Edge %d - EG_getBodyTopos = %d (EG_setTessEdge)!\n",
             index, stat);
    return stat;
  }
  if ((n != 0) && (objs != NULL))
    for (i = 0; i < n; i++) {
      j = EG_indexBodyTopo(obj, objs[i]);
      if (j <= EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Edge %d - indexTopoBody %d = %d (EG_setTessEdge)!\n",
                 index, i+1, j);
        EG_free(objs);
        return j;
      }
      if (btess->tess2d[j-1].npts != 0) {
        if (outLevel > 0)
          printf(" EGADS Error: Edge %d - Face %d set (EG_setTessEdge)!\n",
                 index, j);
        EG_free(objs);
        return EGADS_EXISTS;
      }
    }
  if (objs != NULL) EG_free(objs);

  /* get the bounding information */
  stat = EG_getTopology(btess->tess1d[index-1].obj, &geom, &oclass, &mtype,
                        trange, &nnode, &nodes, &senses);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_getTopology = %d for Edge %d (EG_setTessEdge)!\n",
             stat, index);
    return stat;
  }

  stat = EG_getTopology(nodes[0], &geom, &oclass, &mtype, xyz0, &i, &objs,
                        &senses);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_getTopology = %d - Edge %d/0 (EG_setTessEdge)!\n",
             stat, index);
    return stat;
  }
  j = nnode - 1;
  stat = EG_getTopology(nodes[j], &geom, &oclass, &mtype, xyz1, &i, &objs,
                        &senses);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_getTopology = %d - Edge %d/1 (EG_setTessEdge)!\n",
             stat, index);
    return stat;
  }

  /* allocate the data */
  xyzs = (double *) EG_alloc(3*len*sizeof(double));
  if (xyzs == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Allocating %d Coordinates (EG_setTessEdge)!\n",
             len);
    return EGADS_MALLOC;
  }
  ts = (double *) EG_alloc(len*sizeof(double));
  if (ts == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Allocating %d Parameters (EG_setTessEdge)!\n",
             len);
    EG_free(xyzs);
    return EGADS_MALLOC;
  }
  ts[0]   = trange[0];
  xyzs[0] = xyz0[0];
  xyzs[1] = xyz0[1];
  xyzs[2] = xyz0[2];
  for (i = 1; i < len-1; i++) {
    ts[i]       = t[i];
    xyzs[3*i  ] = xyz[3*i  ];
    xyzs[3*i+1] = xyz[3*i+1];
    xyzs[3*i+2] = xyz[3*i+2];
  }
  ts[len-1]     = trange[1];
  xyzs[3*len-3] = xyz1[0];
  xyzs[3*len-2] = xyz1[1];
  xyzs[3*len-1] = xyz1[2];

  /* set the data */
  if (btess->tess1d[index-1].xyz != NULL) EG_free(btess->tess1d[index-1].xyz);
  if (btess->tess1d[index-1].t   != NULL) EG_free(btess->tess1d[index-1].t);

  btess->tess1d[index-1].npts = len;
  btess->tess1d[index-1].xyz  = xyzs;
  btess->tess1d[index-1].t    = ts;

  if (n > 0) {
    if (btess->tess1d[index-1].faces[0].tric == NULL)
      btess->tess1d[index-1].faces[0].tric = (int *)
                                             EG_alloc((n*(len-1))*sizeof(int));
    if (btess->tess1d[index-1].faces[0].tric == NULL) {
      if (outLevel > 0)
        printf(" EGADS Warning: Alloc %d Tric- Edge %d (EG_setTessEdge)!\n",
               len, index);
    } else {
      for (i = 0; i < len-1; i++)
        for (k = 0; k < n; k++)
          btess->tess1d[index-1].faces[0].tric[i*n+k] = 0;
    }

    if (btess->tess1d[index-1].faces[1].tric == NULL)
      btess->tess1d[index-1].faces[1].tric = (int *)
                                             EG_alloc((n*(len-1))*sizeof(int));
    if (btess->tess1d[index-1].faces[1].tric == NULL) {
      if (outLevel > 0)
        printf(" EGADS Warning: Alloc %d Tric+ Edge %d (EG_setTessEdge)!\n",
               len, index);
    } else {
      for (i = 0; i < len-1; i++)
        for (k = 0; k < n; k++)
          btess->tess1d[index-1].faces[1].tric[i*n+k] = 0;
    }
  }

  return EGADS_SUCCESS;
}


static int
findPoint(double *range, double *uvp, int ptype, int pindex,
          int len, int *table, const double *uv, double *uvx)
{
  int    i, j;
  double du, dv;

  for (i = 0; i < len; i++)
    if ((ptype == table[2*i  ]) && (pindex == table[2*i+1])) {
      du = fabs(uv[2*i  ] - uvp[0])/(range[1]-range[0]);
      dv = fabs(uv[2*i+1] - uvp[1])/(range[3]-range[2]);
      if ((du < 0.25) && (dv < 0.25)) {
        uvx[0] = uv[2*i  ];
        uvx[1] = uv[2*i+1];
        return i;
      }
    }

  /* Degenerate Node? -- update and let pass */
  if (ptype == 0) {
    for (j = i = 0; i < len; i++)
      if ((ptype == table[2*i  ]) && (pindex == table[2*i+1])) j++;
    if (j == 1)
      for (i = 0; i < len; i++)
        if ((ptype == table[2*i  ]) && (pindex == table[2*i+1])) {
/*        printf(" EGADS Info: Differing UV @ Node %d  %lf %lf -- %lf %lf!\n",
                 table[2*i+1], uvp[0], uvp[1], uv[2*i  ], uv[2*i+1]);  */
          du     = fabs(uv[2*i  ] - uvp[0])/(range[1]-range[0]);
          dv     = fabs(uv[2*i+1] - uvp[1])/(range[3]-range[2]);
          uvx[0] = uv[2*i  ];
          uvx[1] = uv[2*i+1];
          if (du >= 0.25) uvx[0] = uvp[0];
          if (dv >= 0.25) uvx[1] = uvp[1];
          return i;
        }
  }

  return EGADS_NOTFOUND;
}


static int
makeNeighbors(int f, int nverts, int ntri, int *tris, int *tric,
              int nseg, triSeg *segs)
{
  int     *ntab, nside, j;
  connect *etab;

  ntab = (int *) EG_alloc(nverts*sizeof(int));
  if (ntab == NULL) {
    printf(" EGADS Error: Vert Table Malloc (EG_setTessFace)!\n");
    return EGADS_MALLOC;
  }
  etab = (connect *) EG_alloc(ntri*3*sizeof(connect));
  if (etab == NULL) {
    printf(" EGADS Error: Edge Table Malloc (EG_setTessFace)!\n");
    EG_free(ntab);
    return EGADS_MALLOC;
  }

  nside = -1;
  for (j = 0; j < nverts; j++) ntab[j] = NOTFILLED;
  for (j = 0; j < ntri;  j++) {
    EG_makeConnect(tris[3*j+1], tris[3*j+2], &tric[3*j  ], &nside,ntab,etab, f);
    EG_makeConnect(tris[3*j  ], tris[3*j+2], &tric[3*j+1], &nside,ntab,etab, f);
    EG_makeConnect(tris[3*j  ], tris[3*j+1], &tric[3*j+2], &nside,ntab,etab, f);
  }

  for (j = 0; j < nseg; j++)
    EG_makeConnect(segs[j].indices[0], segs[j].indices[1], &segs[j].neighbor,
                   &nside, ntab, etab, f);

  /* report any unconnected triangle sides */
  for (j = 0; j <= nside; j++) {
    if (etab[j].tri == NULL) continue;
    printf(" EGADS Info: Face %d, Unconnected Side %d %d = %d\n",
           f, etab[j].node1+1, etab[j].node2+1, *etab[j].tri);
    *etab[j].tri = 0;
  }

  EG_free(etab);
  EG_free(ntab);
  return EGADS_SUCCESS;
}


int
EG_setTessFace(const egObject *tess, int index, int len, const double *xyz,
               const double *uv, int ntri, const int *tris)
{
  int      i, j, k, m, n, hit, iedge, outLevel, stat, nedge, *table, *map;
  int      oclass, mtype, nloop, np, sen, or, lor, pt, pi, *senses, *lsenses;
  int      ntot, st, nseg, ntrix, nf8, nd, mm, mp;
  int      *frlps, *frame, *ptype, *pindex, *trix, *tric, *sns;
  double   smallu, smallv;
  double   range[4], trange[2], uvm[2], uvp[2], uvx[2], *uvs, *xyzs, *intEdg;
  triSeg   *segs;
  fillArea fast;
  egTessel *btess;
  egObject *obj, *geom, *face, **faces, **loops, **edges, **nds;
  static double scl[3][2] = {{1.0, 1.0},  {10.0, 1.0},  {0.1, 10.0}};

  if  (tess == NULL)                 return EGADS_NULLOBJ;
  if  (tess->magicnumber != MAGIC)   return EGADS_NOTOBJ;
  if  (tess->oclass != TESSELLATION) return EGADS_NOTTESS;
  if ((len <= 1) || (ntri < 1))      return EGADS_NODATA;
  if ((xyz == NULL) || (uv == NULL)) return EGADS_NODATA;
  if  (tris == NULL)                 return EGADS_NODATA;
  outLevel = EG_outLevel(tess);

  btess = (egTessel *) tess->blind;
  if (btess == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Blind Object (EG_setTessFace)!\n");
    return EGADS_NOTFOUND;
  }
  obj = btess->src;
  if (obj == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Source Object (EG_setTessFace)!\n");
    return EGADS_NULLOBJ;
  }
  if (obj->magicnumber != MAGIC) {
    if (outLevel > 0)
      printf(" EGADS Error: Source Not an Object (EG_setTessFace)!\n");
    return EGADS_NOTOBJ;
  }
  if (obj->oclass != BODY) {
    if (outLevel > 0)
      printf(" EGADS Error: Source Not Body (EG_setTessFace)!\n");
    return EGADS_NOTBODY;
  }
  if (btess->tess2d == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: No Face Tessellations (EG_setTessFace)!\n");
    return EGADS_NODATA;
  }
  if (btess->done == 1) {
    if (outLevel > 0)
      printf(" EGADS Error: Complete Tessellation (EG_setTessFace)!\n");
    return EGADS_EXISTS;
  }
  if ((index < 1) || (index > btess->nFace)) {
    if (outLevel > 0)
      printf(" EGADS Error: Index = %d [1-%d] (EG_setTessFace)!\n",
             index, btess->nFace);
    return EGADS_INDEXERR;
  }

  /* check triangle indices */
  for (i = 0; i < ntri; i++) {
    for (j = 0; j < 3; j++) {
      if ((tris[3*i+j] < 1) || (tris[3*i+j] > len)) {
        printf(" EGADS Error: Face %d - tris %d/%d = %d [1-%d] (EG_setTessFace)!\n",
               index, i+1, j, tris[3*i+j], len);
        return EGADS_INDEXERR;
      }
    }
    if ((tris[3*i  ] == tris[3*i+1]) || (tris[3*i  ] == tris[3*i+2]) ||
        (tris[3*i+1] == tris[3*i+2])) {
      printf(" EGADS Error: Face %d - tris %d is degenerate = %d %d %d (EG_setTessFace)!\n",
             index, i+1, tris[3*i  ], tris[3*i+1], tris[3*i+2]);
      return EGADS_INDEXERR;
    }
  }

  /* get our Face object */
  stat = EG_getBodyTopos(obj, NULL, FACE, &i, &faces);
  if (stat  != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: Face %d - EG_getBodyTopos Faces = %d (EG_setTessFace)!\n",
             index, stat);
    return stat;
  }
  face = faces[index-1];
  EG_free(faces);

  /* make sure we have all of the edge tessellations */
  stat = EG_getBodyTopos(obj, face, EDGE, &nedge, &edges);
  if (stat  != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: Face %d - EG_getBodyTopos Edges = %d (EG_setTessFace)!\n",
             index, stat);
    return stat;
  }
  for (i = 0; i < nedge; i++) {
    iedge = EG_indexBodyTopo(obj, edges[i]);
    if (iedge <= 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Face %d - EG_indexTopoBody Edge %d = %d (EG_setTessFace)!\n",
               index, i+1, iedge);
      return iedge;
    }
    if (edges[i]->mtype == DEGENERATE) continue;
    if (btess->tess1d[iedge-1].nodes[0] == -btess->tess1d[iedge-1].nodes[1])
      continue;
    if (btess->tess1d[iedge-1].npts > 0) continue;
    printf(" EGADS Error: Face %d - No Tessellation for Edge %d (EG_setTessFace)!\n",
           index, iedge);
    EG_free(edges);
    return EGADS_NOTFOUND;
  }

  /* setup ptype/pindex */
  table = (int *) EG_alloc(3*len*sizeof(int));
  if (table == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Face %d - Allocating %d Vert Table (EG_setTessFace)!\n",
             index, len);
    EG_free(edges);
    return EGADS_MALLOC;
  }
  map = &table[2*len];
  for (i = 0; i < len; i++)
    table[2*i  ] = table[2*i+1] = map[i] = -1;
  for (k = 0; k < len; k++) {
    for (hit = i = 0; i < nedge; i++) {
      iedge = EG_indexBodyTopo(obj, edges[i]);
      if (edges[i]->mtype == DEGENERATE) continue;
      if (btess->tess1d[iedge-1].nodes[0] == -btess->tess1d[iedge-1].nodes[1])
        continue;
      for (j = 0; j < btess->tess1d[iedge-1].npts; j++) {
        if (xyz[3*k  ] != btess->tess1d[iedge-1].xyz[3*j  ]) continue;
        if (xyz[3*k+1] != btess->tess1d[iedge-1].xyz[3*j+1]) continue;
        if (xyz[3*k+2] != btess->tess1d[iedge-1].xyz[3*j+2]) continue;
        table[2*k  ] = j+1;
        table[2*k+1] = iedge;
        if (j == 0) {
          table[2*k  ] = 0;
          table[2*k+1] = btess->tess1d[iedge-1].nodes[0];
        }
        if (j == btess->tess1d[iedge-1].npts-1) {
          table[2*k  ] = 0;
          table[2*k+1] = btess->tess1d[iedge-1].nodes[1];
        }
        hit++;
        break;
      }
      if (hit != 0) break;
    }
  }
  EG_free(edges);

  /* reorder based on loops */
  stat = EG_getTopology(face, &geom, &oclass, &or, range, &nloop, &loops,
                        &lsenses);
  if (stat  != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: Face %d - EG_getTopology = %d (EG_setTessFace)!\n",
             index, stat);
    EG_free(table);
    return stat;
  }
  smallu = 0.00005*(range[1] - range[0]);
  smallv = 0.00005*(range[3] - range[2]);

  /* get total number of points in all of the loops */
  for (ntot = i = 0; i < nloop; i++) {
    stat = EG_getTopology(loops[i], &geom, &oclass, &mtype, NULL, &nedge,
                          &edges, &senses);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Face %d - EG_getTopology Loop %d = %d (EG_setTessFace)!\n",
               index, i+1, stat);
      EG_free(table);
      return stat;
    }
    for (j = 0; j < nedge; j++) {
      iedge = EG_indexBodyTopo(obj, edges[j]);
      if (edges[j]->mtype == DEGENERATE) continue;
      if (btess->tess1d[iedge-1].nodes[0] == -btess->tess1d[iedge-1].nodes[1])
        continue;
      ntot += btess->tess1d[iedge-1].npts-1;
    }
  }
  ntrix = ntot-2 + 2*(nloop-1);
  segs  = (triSeg *) EG_alloc(ntot*sizeof(triSeg));
  if (segs == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Face %d - Allocating %d Segs (EG_setTessFace)!\n",
             index, ntot);
    EG_free(table);
    return EGADS_MALLOC;
  }
  uvs = (double *) EG_alloc((2*ntot+2)*sizeof(double));
  if (uvs == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Face %d - Allocating %d uvs (EG_setTessFace)!\n",
             index, ntot);
    EG_free(segs);
    EG_free(table);
    return EGADS_MALLOC;
  }
  uvs[0] = uvs[1] = 0.0;

  /* find the Edge vertices */
  frlps = (int *) EG_alloc(nloop*sizeof(int));
  if (frlps == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Face %d - Allocating %d Loops (EG_setTessFace)!\n",
             index, nloop);
    EG_free(uvs);
    EG_free(segs);
    EG_free(table);
    return EGADS_MALLOC;
  }
  for (np = i = 0; i < nloop; i++) {
    st   = np;
    stat = EG_getTopology(loops[i], &geom, &oclass, &mtype, NULL, &nedge,
                          &edges, &senses);
    if (stat != EGADS_SUCCESS) continue;
    lor = 1;
    if ((lsenses[i] == 2) || (lsenses[i] == -2)) lor = -1;
    n = 0;
    if (or*lor == SREVERSE) n = nedge-1;
    for (j = 0; j < nedge; j++, n += or*lor) {
      iedge = EG_indexBodyTopo(obj, edges[n]);
      if (edges[n]->mtype == DEGENERATE) continue;
      if (btess->tess1d[iedge-1].nodes[0] == -btess->tess1d[iedge-1].nodes[1])
        continue;
      stat = EG_getTopology(edges[n], &geom, &oclass, &mtype, trange, &nd,
                            &nds, &sns);
      if (stat != EGADS_SUCCESS) {
        printf(" EGADS Error: EG_getTopology = %d  for Face = %d, Edge = %d\n",
               stat, index, iedge);
        EG_free(frlps);
        EG_free(uvs);
        EG_free(segs);
        EG_free(table);
        return stat;
      }
      sen = senses[n]*or*lor;

      /* internal Edge? */
      intEdg = NULL;
      for (m = 0; m < nedge; m++) {
        if (m == n) continue;
        if (iedge == EG_indexBodyTopo(obj, edges[m])) {
          uvm[0] = uvm[1] = -1.0;
          uvp[0] = uvp[1] =  1.0;
          EG_getEdgeUV(face, edges[n], -1, 0.5*(trange[0]+trange[1]), uvm);
          EG_getEdgeUV(face, edges[n],  1, 0.5*(trange[0]+trange[1]), uvp);
          if ((uvm[0] == uvp[0]) && (uvm[1] == uvp[1]) && (intEdg == NULL)) {
            printf(" EGADS Info: ");
            printf("Face #%d -> Edge #%d (%d) Internally in Loop %d %d, sen = %d!\n",
                   index, iedge, nedge, n+1, m+1, sen);
            intEdg = (double *)
                     EG_alloc(4*btess->tess1d[iedge-1].npts*sizeof(double));
            if (intEdg == NULL) {
              printf(" EGADS Internal: Cannot Allocate %d intEgdes!\n",
                     btess->tess1d[iedge-1].npts);
              continue;
            }
            for (m = 0; m < btess->tess1d[iedge-1].npts; m++) {
              stat = EG_getEdgeUV(face, edges[n], senses[n]*lor,
                                  btess->tess1d[iedge-1].t[m], &intEdg[4*m]);
              if (stat != EGADS_SUCCESS) {
                printf(" EGADS Error: getEdgeUV! = %d  for Face %d, Edge = %d\n",
                       stat, index, iedge);
                EG_free(intEdg);
                EG_free(frlps);
                EG_free(uvs);
                EG_free(segs);
                EG_free(table);
                return stat;
              }
            }
            for (m = 0; m < btess->tess1d[iedge-1].npts; m++) {
              mm = m - 1;
              mp = m + 1;
              if (mm <  0) mm = 0;
              if (mp >= btess->tess1d[iedge-1].npts)
                mp = btess->tess1d[iedge-1].npts - 1;
              uvm[0] = intEdg[4*mp  ] - intEdg[4*mm  ];
              uvm[1] = intEdg[4*mp+1] - intEdg[4*mm+1];
              uvp[0] = atan2(-uvm[1], uvm[0]);
              intEdg[4*m+2] = sen*smallu*sin(uvp[0]);
              intEdg[4*m+3] = sen*smallv*cos(uvp[0]);
            }
          }
        }
      }

      if (sen == 1) {
        for (m = 0; m < btess->tess1d[iedge-1].npts-1; m++, np++) {
          pt = m+1;
          pi = iedge;
          if (m == 0) {
            pt = 0;
            pi = btess->tess1d[iedge-1].nodes[0];
          }
          uvp[0] = uvp[1] = 0.0;
          stat   = EG_getEdgeUV(face, edges[n], senses[n],
                                btess->tess1d[iedge-1].t[m], uvp);
          if (stat != EGADS_SUCCESS)
            printf(" EGADS Internal: Face %d - EdgeUV+ = %d (EG_setTessFace)!\n",
                   index, stat);
          stat = findPoint(range, uvp, pt, pi, len, table, uv, uvx);
          if (stat < EGADS_SUCCESS) {
            printf(" EGADS Error: Face %d - FindPt+ %d/%d = %d (EG_setTessFace)!\n",
                   index, pt, pi, stat);
            if (intEdg != NULL) EG_free(intEdg);
            EG_free(frlps);
            EG_free(uvs);
            EG_free(segs);
            EG_free(table);
            return stat;
          }
          map[stat]           =  np;
          if (intEdg == NULL) {
            uvs[2*np+2]       =  uvx[0];
            uvs[2*np+3]       =  uvx[1];
          } else {
            uvs[2*np+2]       =  intEdg[4*m  ] + intEdg[4*m+2];
            uvs[2*np+3]       =  intEdg[4*m+1] + intEdg[4*m+3];
          }
          segs[np].indices[0] =  np+1;
          segs[np].indices[1] =  np+2;
          segs[np].neighbor   = -iedge;
          segs[np].edge       =  senses[n]*lor*iedge;
          segs[np].index      =  m+1;
        }
      } else {
        for (m = btess->tess1d[iedge-1].npts-1; m > 0; m--, np++) {
          pt = m+1;
          pi = iedge;
          if (m == btess->tess1d[iedge-1].npts-1) {
            pt = 0;
            pi = btess->tess1d[iedge-1].nodes[1];
          }
          uvm[0] = uvm[1] = 0.0;
          stat   = EG_getEdgeUV(face, edges[n], senses[n],
                                btess->tess1d[iedge-1].t[m], uvm);
          if (stat != EGADS_SUCCESS)
            printf(" EGADS Internal: Face %d - EdgeUV- = %d (EG_setTessFace)!\n",
                   index, stat);
          stat = findPoint(range, uvm, pt, pi, len, table, uv, uvx);
          if (stat < EGADS_SUCCESS) {
            printf(" EGADS Error: Face %d - FindPt- %d/%d = %d (EG_setTessFace)!\n",
                   index, pt, pi, stat);
            if (intEdg != NULL) EG_free(intEdg);
            EG_free(frlps);
            EG_free(uvs);
            EG_free(segs);
            EG_free(table);
            return stat;
          }
          map[stat]           =  np;
          if (intEdg == NULL) {
            uvs[2*np+2]       =  uvx[0];
            uvs[2*np+3]       =  uvx[1];
          } else {
            uvs[2*np+2]       =  intEdg[4*m  ] + intEdg[4*m+2];
            uvs[2*np+3]       =  intEdg[4*m+1] + intEdg[4*m+3];
          }
          segs[np].indices[0] =  np+1;
          segs[np].indices[1] =  np+2;
          segs[np].neighbor   = -iedge;
          segs[np].edge       =  senses[n]*lor*iedge;
          segs[np].index      =  m+1;
        }
      }
      if (intEdg != NULL) EG_free(intEdg);
    }
    if (np > 0) segs[np-1].indices[1] = st+1;
    frlps[i] = np - st;
  }
  nseg = np;

  /* fill up 2D tess structure */

  for (i = 0; i < len; i++)
    if (map[i] == -1) {
      map[i] = np;
      np++;
    }

  /* get the frame */
  fast.pts   = NULL;
  fast.segs  = NULL;
  fast.front = NULL;
  frame = (int *) EG_alloc(3*ntrix*sizeof(int));
  if (frame == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Face %d - Allocating %d frame (EG_setTessFace)!\n",
             index, ntrix);
    EG_free(frlps);
    EG_free(uvs);
    EG_free(segs);
    EG_free(table);
    return EGADS_MALLOC;
  }
  n = EG_fillArea(nloop, frlps, uvs, frame, &nf8, 0, &fast);
  /* adjust for figure 8 configurations */
  if (nf8 != 0) {
    printf(" EGADS Warning: Face %d -> Found %d figure 8's!\n", index, nf8);
    for (i = 0; i < nf8; i++) if (n+2*i == ntrix) ntrix = n;
  }
  if (n != ntrix) {
    range[0] = range[2] = uvs[2];
    range[1] = range[3] = uvs[3];
    for (i = 2; i <= ntot; i++) {
      if (uvs[2*i  ] < range[0]) range[0] = uvs[2*i  ];
      if (uvs[2*i+1] < range[1]) range[1] = uvs[2*i+1];
      if (uvs[2*i  ] > range[2]) range[2] = uvs[2*i  ];
      if (uvs[2*i+1] > range[3]) range[3] = uvs[2*i+1];
    }
    for (i = 1; i <= ntot; i++) {
      uvs[2*i  ] = (uvs[2*i  ]-range[0])/(range[2]-range[0]);
      uvs[2*i+1] = (uvs[2*i+1]-range[1])/(range[3]-range[1]);
    }
    for (j = 0; j < 3; j++) {
      for (i = 1; i <= ntot; i++) {
        uvs[2*i  ] *= scl[j][0];
        uvs[2*i+1] *= scl[j][1];
      }
      n = EG_fillArea(nloop, frlps, uvs, frame, &nf8, 1, &fast);
      printf(" EGADS Internal: Face %d -> Renormalizing %d, ntris = %d (%d)!\n",
             index, j, ntrix, n);
      if (n == ntrix) break;
    }
  }
  if (fast.segs  != NULL) EG_free(fast.segs);
  if (fast.pts   != NULL) EG_free(fast.pts);
  if (fast.front != NULL) EG_free(fast.front);
  EG_free(uvs);
  if (n != ntrix) {
    printf(" EGADS Error: Face %d - Can't Triangulate Frame (EG_setTessFace)!\n",
           index);
    EG_free(frame);
    EG_free(frlps);
    EG_free(segs);
    EG_free(table);
    return EGADS_DEGEN;
  }

  /* set the triangle data */
  trix = (int *) EG_alloc(3*ntri*sizeof(int));
  tric = (int *) EG_alloc(3*ntri*sizeof(int));
  if ((trix == NULL) || (tric == NULL)) {
    if (trix != NULL) EG_free(trix);
    if (tric != NULL) EG_free(tric);
    if (outLevel > 0)
      printf(" EGADS Error: Face %d - Allocating %d tris (EG_setTessFace)!\n",
             index, ntri);
    EG_free(frame);
    EG_free(frlps);
    EG_free(segs);
    EG_free(table);
    return EGADS_MALLOC;
  }
  for (i = 0; i < ntri; i++) {
    trix[3*i  ] = map[tris[3*i  ]-1] + 1;
    trix[3*i+1] = map[tris[3*i+1]-1] + 1;
    trix[3*i+2] = map[tris[3*i+2]-1] + 1;
    tric[3*i  ] = i+1;
    tric[3*i+1] = i+1;
    tric[3*i+2] = i+1;
  }
  stat = makeNeighbors(index, len, ntri, trix, tric, nseg, segs);
  EG_free(segs);
  if (stat != EGADS_SUCCESS) {
    EG_free(trix);
    EG_free(tric);
    EG_free(frame);
    EG_free(frlps);
    EG_free(segs);
    EG_free(table);
    return stat;
  }

  /* set the reordered vertices */
  ptype  = (int *)    EG_alloc(  len*sizeof(int));
  pindex = (int *)    EG_alloc(  len*sizeof(int));
  uvs    = (double *) EG_alloc(2*len*sizeof(double));
  xyzs   = (double *) EG_alloc(3*len*sizeof(double));
  if ((ptype == NULL) || (pindex == NULL) || (xyzs == NULL) || (uvs == NULL)) {
    if (ptype  != NULL) EG_free(ptype);
    if (pindex != NULL) EG_free(pindex);
    if (uvs    != NULL) EG_free(uvs);
    if (xyzs   != NULL) EG_free(xyzs);
    if (outLevel > 0)
      printf(" EGADS Error: Face %d - Allocating %d verts (EG_setTessFace)!\n",
             index, len);
    EG_free(trix);
    EG_free(tric);
    EG_free(frame);
    EG_free(frlps);
    EG_free(table);
    return EGADS_MALLOC;
  }
  for (j = 0; j < len; j++) {
    i           = map[j];
    ptype[i]    = table[2*j  ];
    pindex[i]   = table[2*j+1];
    uvs[2*i  ]  = uv[2*j  ];
    uvs[2*i+1]  = uv[2*j+1];
    xyzs[3*i  ] = xyz[3*j  ];
    xyzs[3*i+1] = xyz[3*j+1];
    xyzs[3*i+2] = xyz[3*j+2];
  }
  EG_free(table);

  /* update the Face pointers */
  if (btess->tess2d[index-1].xyz    != NULL)
    EG_free(btess->tess2d[index-1].xyz);
  if (btess->tess2d[index-1].uv     != NULL)
    EG_free(btess->tess2d[index-1].uv);
  if (btess->tess2d[index-1].ptype  != NULL)
    EG_free(btess->tess2d[index-1].ptype);
  if (btess->tess2d[index-1].pindex != NULL)
    EG_free(btess->tess2d[index-1].pindex);
  if (btess->tess2d[index-1].bary   != NULL)
    EG_free(btess->tess2d[index-1].bary);
  if (btess->tess2d[index-1].frame  != NULL)
    EG_free(btess->tess2d[index-1].frame);
  if (btess->tess2d[index-1].frlps  != NULL)
    EG_free(btess->tess2d[index-1].frlps);
  if (btess->tess2d[index-1].tris   != NULL)
    EG_free(btess->tess2d[index-1].tris);
  if (btess->tess2d[index-1].tric   != NULL)
    EG_free(btess->tess2d[index-1].tric);
  btess->tess2d[index-1].npts   = len;
  btess->tess2d[index-1].xyz    = xyzs;
  btess->tess2d[index-1].uv     = uvs;
  btess->tess2d[index-1].ptype  = ptype;
  btess->tess2d[index-1].pindex = pindex;
  btess->tess2d[index-1].ntris  = ntri;
  btess->tess2d[index-1].tris   = trix;
  btess->tess2d[index-1].tric   = tric;
  btess->tess2d[index-1].bary   = NULL;
  btess->tess2d[index-1].nframe = ntrix;
  btess->tess2d[index-1].frame  = frame;
  btess->tess2d[index-1].frlps  = frlps;
  btess->tess2d[index-1].nfrlps = nloop;

  for (i = 1; i < nloop; i++)
    btess->tess2d[index-1].frlps[i] += btess->tess2d[index-1].frlps[i-1];

  return EGADS_SUCCESS;
}


int
EG_localToGlobal(const egObject *tess, int index, int local, int *global)
{
  int      stat;
  egTessel *btess;

  if  (tess == NULL)                          return EGADS_NULLOBJ;
  if  (tess->magicnumber != MAGIC)            return EGADS_NOTOBJ;
  if  (tess->oclass != TESSELLATION)          return EGADS_NOTTESS;
  if  (index == 0)                            return EGADS_INDEXERR;
  if  (local <  1)                            return EGADS_RANGERR;
  btess = (egTessel *) tess->blind;
  if (btess == NULL)                          return EGADS_NOTFOUND;
  if (btess->done == 0)                       return EGADS_TESSTATE;
  if ((index < 0) && (-index > btess->nEdge)) return EGADS_INDEXERR;
  if ((index > 0) && ( index > btess->nFace)) return EGADS_INDEXERR;

  if (btess->globals == NULL) {
    stat = EG_computeTessMap(btess, EG_outLevel(tess));
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: EG_computeTessMap = %d (EG_localToGlobal)!\n", stat);
      return stat;
    }
  }

  if (index < 0) {
    if (btess->tess1d[-index-1].global == NULL) return EGADS_DEGEN;
    if (local > btess->tess1d[-index-1].npts) return EGADS_RANGERR;
    *global = btess->tess1d[-index-1].global[local-1];
  } else {
    if (btess->tess2d[ index-1].global == NULL) return EGADS_DEGEN;
    if (local > btess->tess2d[ index-1].npts) return EGADS_RANGERR;
    *global = btess->tess2d[ index-1].global[local-1];
  }

  return EGADS_SUCCESS;
}


int
EG_getGlobal(const egObject *tess, int global, int *ptype, int *pindex,
             /*@null@*/ double *xyz)
{
  int      i, j, stat;
  egTessel *btess;

  if (tess == NULL)                 return EGADS_NULLOBJ;
  if (tess->magicnumber != MAGIC)   return EGADS_NOTOBJ;
  if (tess->oclass != TESSELLATION) return EGADS_NOTTESS;
  btess = (egTessel *) tess->blind;
  if (btess == NULL)                return EGADS_NOTFOUND;
  if (btess->done == 0)             return EGADS_TESSTATE;

  if (btess->globals == NULL) {
    stat = EG_computeTessMap(btess, EG_outLevel(tess));
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: EG_computeTessMap = %d (EG_getGlobal)!\n", stat);
      return stat;
    }
  }
  if ((global < 1) || (global > btess->nGlobal)) return EGADS_INDEXERR;

  *ptype  = i = btess->globals[2*global-2];
  *pindex = j = btess->globals[2*global-1];
  if (xyz == NULL) return EGADS_SUCCESS;

  if (i == 0) {
    xyz[0] = btess->xyzs[3*j-3];
    xyz[1] = btess->xyzs[3*j-2];
    xyz[2] = btess->xyzs[3*j-1];
  } else if (i > 0) {
    xyz[0] = btess->tess1d[j-1].xyz[3*i-3];
    xyz[1] = btess->tess1d[j-1].xyz[3*i-2];
    xyz[2] = btess->tess1d[j-1].xyz[3*i-1];
  } else {
    xyz[0] = btess->tess2d[j-1].xyz[-3*i-3];
    xyz[1] = btess->tess2d[j-1].xyz[-3*i-2];
    xyz[2] = btess->tess2d[j-1].xyz[-3*i-1];
  }

  return EGADS_SUCCESS;
}


static int
EG_findMidSide(int i1, int i2, int *table, midside *mid)
{
  int index, last;

  if (table[i1] == NOTFILLED) return 0;

  index = last = table[i1];
  while (index != NOTFILLED) {
    if (mid[index].vert2 == i2) return mid[index].nvert;
    last  = index;
    index = mid[index].next;
  }

  return -last-1;
}


void
EG_getInterior(const ego face, double *xyz, double *uv)
{
  int    count, stat, diverge = 0;
  double a00, a10, a11, b0, b1, det, dist, ldist, frac = 1.0;
  double dx[3], du[2], range[4], result[18];

  stat = EG_getRange(face, range, &count);
  if (stat != EGADS_SUCCESS) {
    printf(" EGADS Warning: EG_getRange = %d (EG_getInterior)!\n", stat);
    return;
  }
  if ((uv[0] < range[0]) || (uv[0] > range[1]) || (uv[1] < range[2]) ||
      (uv[1] > range[3])) {
    printf(" EGADS Warning: Out of Range (EG_getInterior)!\n");
    return;
  }
  
  /* newton iteration */
  ldist = 0.0;
  for (count = 0; count < 15; count++) {
    stat = EG_evaluate(face, uv, result);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Warning: EG_evaluate = %d (EG_getInterior)!\n", stat);
      return;
    }
    dx[0] = result[0] - xyz[0];
    dx[1] = result[1] - xyz[1];
    dx[2] = result[2] - xyz[2];
    dist  = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
    if (dist < EPS10) return;

    b0  =    -dx[0]*result[ 3] -     dx[1]*result[ 4] -     dx[2]*result[ 5];
    b1  =    -dx[0]*result[ 6] -     dx[1]*result[ 7] -     dx[2]*result[ 8];
    a00 = result[3]*result[ 3] + result[4]*result[ 4] + result[5]*result[ 5] +
              dx[0]*result[ 9] +     dx[1]*result[10] +     dx[2]*result[11];
    a10 = result[3]*result[ 6] + result[4]*result[ 7] + result[5]*result[ 8] +
              dx[0]*result[12] +     dx[1]*result[13] +     dx[2]*result[14];
    a11 = result[6]*result[ 6] + result[7]*result[ 7] + result[8]*result[ 8] +
              dx[0]*result[15] +     dx[1]*result[16] +     dx[2]*result[17];

    det = a00*a11 - a10*a10;
    if (det == 0.0) {
      printf(" EGADS Warning: det = zero (EG_getInterior)!\n");
      return;
    }
    det   = 1.0/det;
    du[0] = det*(b0*a11 - b1*a10);
    du[1] = det*(b1*a00 - b0*a10);
    if ((fabs(du[0]) < EPS10) && (fabs(du[1]) < EPS10)) return;
    if (count != 0)
      if (ldist < dist) {
        printf(" EG_getInterior -- diverge %d  %le %le!\n",
               count+1, du[0], du[1]);
/*      return;  */
        diverge++;
      }
    ldist  = dist;
    
    /* partial update? */
    if (uv[0]+frac*du[0] < range[0]) frac = (uv[0]-range[0])/fabs(du[0]);
    if (uv[0]+frac*du[0] > range[1]) frac = (range[1]-uv[0])/fabs(du[0]);
    if (uv[1]+frac*du[1] < range[2]) frac = (uv[1]-range[2])/fabs(du[1]);
    if (uv[1]+frac*du[1] > range[3]) frac = (range[3]-uv[1])/fabs(du[1]);
    
    uv[0] += frac*du[0];
    uv[1] += frac*du[1];
    if (frac != 1.0) {
      if (diverge != 0) printf(" EG_getInterior: boundary exit %d\n", diverge);
      return;
    }
  }
  
  printf(" EG_getInterior: not converged %d!\n", diverge);
}


#ifdef MINARC

// Find the point P = S(u,v) midpoint of uv0 and uv1 that minimizes the distance between
// min 1/2 ( l_0 ^ 2 + l_1 ^ 2) + lambda * (l_1 - l_0) = L(u, v, lambda)
// grad (L) = (L_1, L_2, L_3) = (0, 0, 0) --> solution using Newton x_n+1 = x_n + delta_n (3x3 system)
void
EG_getSidepoint(const ego face, double fact, const double *uvm,
                const double *uvp, double *uv)
{
  int    i, it, nT = 50;
  double b, dlu0, dlu1, dlv0, dlv1, ddl0, ddl1, detJ, x2, ctt;
  double l[2], pM[18], pP[18], pIT[18], J[3][3], ATJ[3][3], r0[3], r1[3];
  double delta[3], L[3], uvIT[3], xyz[3], range[4];
#ifdef DEBUG
  double e1 = 0.0, e2 = 0.0, x0, x1;
#endif

  i = EG_getRange(face, range, &it);
  if (i  != EGADS_SUCCESS) {
    printf(" EG_getSidepoint: EG_egetRange %d !!\n ", i);
    return;
  }
  /* Initial guess uv  */
  uvIT[0] = fact*uvp[0]  + (1.0-fact)*uvm[0];
  uvIT[1] = fact*uvp[1]  + (1.0-fact)*uvm[1];
  uv[0]   = uvIT[0];
  uv[1]   = uvIT[1];
  uvIT[2] = 0.0;
  i       = EG_evaluate(face, uvm , pM);
  i      += EG_evaluate(face, uvp , pP);
  i      += EG_evaluate(face, uvIT, pIT);
  if (i  != EGADS_SUCCESS) {
    printf(" EG_getSidepoint: EG_evaluate %d !!\n ", i);
    return;
  }

#ifdef DEBUG
  r0[0] = pIT[0] - pM[0];
  r0[1] = pIT[1] - pM[1];
  r0[2] = pIT[2] - pM[2];

  r1[0] = pIT[0] - pP[0];
  r1[1] = pIT[1] - pP[1];
  r1[2] = pIT[2] - pP[2];

  l[0]  = DOT(r0, r0);
  l[1]  = DOT(r1, r1);
  b     = sqrt(l[0]) + sqrt(l[1]);
  printf(" \n\n --------------------------------------------------- \n");
  printf(" INITIAL GUESS: TOTAL ARC %lf l0 = ||pIT - pM||^2 = %lf l1 = ||pIT - pP||^2 = %lf\n ratios %lf %lf (EXPECTED %lf %lf)\n",
         b, l[0], l[1], l[0]/(b*b), l[1]/(b*b), fact*fact, (1.0-fact)*(1.0-fact));
  l[0]  = sqrt(l[0]);
  l[1]  = sqrt(l[1]);
  printf(" As distances sqrt(l0) = ||pIT - pM|| = %lf sqrt (l1) = ||pIT - pP|| = %lf\n ratios %lf %lf (EXPECTED %lf %lf)\n",
         l[0], l[1], l[0]/b, l[1]/b, fact, (1.0-fact));
  printf(" --------------------------------------------------- \n");
#endif

  ctt     = (1.0 - fact) * (1.0 - fact) / (fact * fact);
  for (it = 0; it < nT; it++) {
    r0[0]   = pIT[0] - pM[0];
    r0[1]   = pIT[1] - pM[1];
    r0[2]   = pIT[2] - pM[2];

    r1[0]   = pIT[0] - pP[0];
    r1[1]   = pIT[1] - pP[1];
    r1[2]   = pIT[2] - pP[2];

    l[0]    = DOT(r0, r0);
    l[1]    = DOT(r1, r1);
    dlu0    = 2.0 * (r0[0] * pIT[3] + r0[1] * pIT[4] + r0[2] * pIT[5]);
    dlu1    = 2.0 * (r1[0] * pIT[3] + r1[1] * pIT[4] + r1[2] * pIT[5]);
    dlv0    = 2.0 * (r0[0] * pIT[6] + r0[1] * pIT[7] + r0[2] * pIT[8]);
    dlv1    = 2.0 * (r1[0] * pIT[6] + r1[1] * pIT[7] + r1[2] * pIT[8]);

    /* UPDATE FUNCTION L */
    L[0]    = dlu0 * (l[0] - ctt * uvIT[2]) + dlu1 * (l[1] + uvIT[2]);
    L[1]    = dlv0 * (l[0] - ctt * uvIT[2]) + dlv1 * (l[1] + uvIT[2]);
    L[2]    = l[1] - ctt * l[0];

    /* duu */
    ddl0    = 2.0 * (pIT[9 ] * r0[0] + pIT[3] * pIT[3] +
                     pIT[10] * r0[1] + pIT[4] * pIT[4] +
                     pIT[11] * r0[2] + pIT[5] * pIT[5]);
    ddl1    = 2.0 * (pIT[9 ] * r1[0] + pIT[3] * pIT[3] +
                     pIT[10] * r1[1] + pIT[4] * pIT[4] +
                     pIT[11] * r1[2] + pIT[5] * pIT[5]);
    b       = ddl0 * (l[0] - ctt * uvIT[2]) + ddl1* (l[1] + uvIT[2]);
    J[0][0] = dlu0 * dlu0 + dlu1 * dlu1 + b;

    /* duv */
    ddl0    = 2.0 * (pIT[12] * r0[0] + pIT[3] * pIT[6] +
                     pIT[13] * r0[1] + pIT[4] * pIT[7] +
                     pIT[14] * r0[2] + pIT[5] * pIT[8]);

    ddl1    = 2.0 * (pIT[12] * r1[0] + pIT[3] * pIT[6] +
                     pIT[13] * r1[1] + pIT[4] * pIT[7] +
                     pIT[14] * r1[2] + pIT[5] * pIT[8]);
    b       = ddl0 * (l[0] - ctt * uvIT[2]) + ddl1 * (l[1] + uvIT[2]);
    J[0][1] = dlu0 * dlv0 + dlu1 * dlv1 + b;

    /* dvv */
    ddl0    = 2.0 * (pIT[15] * r0[0] + pIT[6] * pIT[6] +
                     pIT[16] * r0[1] + pIT[7] * pIT[7] +
                     pIT[17] * r0[2] + pIT[8] * pIT[8]);

    ddl1    = 2.0 * (pIT[15] * r1[0] + pIT[6] * pIT[6] +
                     pIT[16] * r1[1] + pIT[7] * pIT[7] +
                     pIT[17] * r1[2] + pIT[8] * pIT[8]);
    b       = ddl0 * (l[0] - ctt * uvIT[2]) + ddl1 * (l[1] + uvIT[2]);
    J[1][1] = dlv0 * dlv0 + dlv1 * dlv1 + b;

    J[0][2] = dlu1 - ctt * dlu0;
    J[1][2] = dlv1 - ctt * dlv0;
    J[1][0] = J[0][1];
    J[2][0] = J[0][2];
    J[2][1] = J[1][2];
    J[2][2] = 0.0;

    /* Solve Linear System: J * delta = - L
       For now: Invert Jacobian directly --> delta = J^-1 L  */
    detJ    = J[0][0] * J[1][1] * J[2][2] + J[0][1] * J[1][2] * J[2][0] +
    J[0][2] * J[2][1] * J[1][0] - J[0][2] * J[1][1] * J[2][0] -
    J[1][2] * J[2][1] * J[0][0] - J[2][2] * J[1][0] * J[0][1];
    if (fabs(detJ) < EPS10) break;

    ATJ[0][0] =   J[1][1] * J[2][2] - J[2][1] * J[1][2];
    ATJ[0][1] = -(J[0][1] * J[2][2] - J[2][1] * J[0][2]);
    ATJ[0][2] =   J[0][1] * J[1][2] - J[1][1] * J[0][2];
    ATJ[1][0] = -(J[1][0] * J[2][2] - J[2][0] * J[1][2]);
    ATJ[1][1] =   J[0][0] * J[2][2] - J[2][0] * J[0][2];
    ATJ[1][2] = -(J[0][0] * J[1][2] - J[1][0] * J[0][2]);
    ATJ[2][0] =   J[1][0] * J[2][1] - J[2][0] * J[1][1];
    ATJ[2][1] = -(J[0][0] * J[2][1] - J[2][0] * J[0][1]);
    ATJ[2][2] =   J[0][0] * J[1][1] - J[1][0] * J[0][1];

    detJ      = 1.0 / detJ;
    delta[0]  = -detJ * (ATJ[0][0] * L[0] + ATJ[0][1] * L[1] + ATJ[0][2] * L[2]);
    delta[1]  = -detJ * (ATJ[1][0] * L[0] + ATJ[1][1] * L[1] + ATJ[1][2] * L[2]);
    delta[2]  = -detJ * (ATJ[2][0] * L[0] + ATJ[2][1] * L[1] + ATJ[2][2] * L[2]);
    uvIT[0]  += delta[0];
    uvIT[1]  += delta[1];
    uvIT[2]  += delta[2];
    if (uvIT[0] < range[0] || uvIT[0] > range[1] ||
        uvIT[1] < range[2] || uvIT[1] > range[3]) {
#ifdef DEBUG
      printf(" UV %lf %lf OUT OF RANGE: %lf %lf %lf %lf\n", uvIT[0], uvIT[1],
             range[0], range[1], range[2], range[3]);
#endif
      pIT[0] = (1.0-fact)*pM[0] + fact*pP[0];
      pIT[1] = (1.0-fact)*pM[1] + fact*pP[1];
      pIT[2] = (1.0-fact)*pM[2] + fact*pP[2];
      i      = EG_invEvaluateGuess(face, pIT, uv, xyz);
      if (i != EGADS_SUCCESS || uv[0] < range[0] || uv[0] > range[1] ||
          uv[1] < range[2] || uv[1] > range[3]) {
        uv[0] = fact*uvp[0]  + (1.0-fact)*uvm[0];
        uv[1] = fact*uvp[1]  + (1.0-fact)*uvm[1];
      }
      return;
    }
    x2 = sqrt(DOT(delta, delta));
#ifdef DEBUG
    if      (it == 0) x0 = x2;
    else if (it == 1) x1 = x2;
    else {
      e1 = fabs(x1 / x0);
      e2 = fabs(x2 / x1);
      x0 = x1;
      x1 = x2;
    }
#endif
    i = EG_evaluate(face, uvIT, pIT);
    if (i != EGADS_SUCCESS || x2 < EPS10) break;
  }
#ifdef DEBUG
  if (i != EGADS_SUCCESS) printf("EG_evaluate %d !!\n", i);
  printf("IT %d DELTA SIZE %1.2e < %1.2e L [%lf  %lf  %lf] \n",
         it, x2, EPS10, L[0], L[1], L[2]);
  if (e1 > EPS10 && e2 > EPS10) printf("CONVERGENCE RATE %lf \n", log(e2)/log(e1));
  b    = sqrt(l[0]) + sqrt(l[1]);
  printf(" --------------------------------------------------- \n");
  printf(" NEW TOTAL ARC %lf with l0 = ||pIT - pM||^2 = %lf l1 = ||PIT - pP||^2 = %lf\n ratios %lf %lf (EXPECTED %lf %lf)\n",
         b, l[0], l[1], l[0]/(b*b), l[1]/(b*b), fact*fact, (1.0-fact)*(1.0-fact));
  l[0] = sqrt(l[0]);
  l[1] = sqrt(l[1]);
  printf(" As Actual chords ||pIT - pM|| = %lf ||pIT - pP|| = %lf\n ratios %lf %lf (EXPECTED %lf %lf)\n",
         l[0], l[1], l[0]/b, l[1]/b, fact, (1.0 - fact));
  printf(" --------------------------------------------------- \n");
#endif

  if (it == nT) printf(" EG_getSidepoint: not converged!\n");
  uv[0] = uvIT[0];
  uv[1] = uvIT[1];
}

#else

void
EG_getSidepoint(const ego face, double fac, const double *uvm,
                const double *uvp, double *uv)
{
  int    stat;
  double result[18], xyzm[3], xyzp[3], xyz[3];

  stat    = EG_evaluate(face, uvm, result);
  if (stat != EGADS_SUCCESS) {
    printf(" EGADS Warning: EG_evaluate M = %d (EG_getSidepoint)!\n", stat);
    return;
  }
  xyzm[0] = result[0];
  xyzm[1] = result[1];
  xyzm[2] = result[2];
  stat    = EG_evaluate(face, uvp, result);
  if (stat != EGADS_SUCCESS) {
    printf(" EGADS Warning: EG_evaluate P = %d (EG_getSidepoint)!\n", stat);
    return;
  }
  xyzp[0] = result[0];
  xyzp[1] = result[1];
  xyzp[2] = result[2];

  uv[0]  = fac*uvp[0]  + (1.0-fac)*uvm[0];
  uv[1]  = fac*uvp[1]  + (1.0-fac)*uvm[1];
  xyz[0] = fac*xyzp[0] + (1.0-fac)*xyzm[0];
  xyz[1] = fac*xyzp[1] + (1.0-fac)*xyzm[1];
  xyz[2] = fac*xyzp[2] + (1.0-fac)*xyzm[2];

  EG_getInterior(face, xyz, uv);
}

#endif


static void
EG_quadThread(void *struc)
{
  int     index, stat;
  long    ID;
  EMPquad *qthread;

  qthread = (EMPquad *) struc;

  /* get our identifier */
  ID = EMP_ThreadID();

  /* look for work */
  for (;;) {

    /* only one thread at a time here -- controlled by a mutex! */
    if (qthread->mutex != NULL) EMP_LockSet(qthread->mutex);
    for (index = qthread->index; index < qthread->end; index++) {
      if (qthread->ntess->tess2d[index].tfi ==    1) continue;
      if (qthread->bodydata->qm[index]      == NULL) continue;
      if (qthread->bodydata->qm[index]->fID ==    0) continue;
      break;
    }
    qthread->index = index+1;
    if (qthread->mutex != NULL) EMP_LockRelease(qthread->mutex);
    if (index >= qthread->end) break;

    /* do the work */
    stat = EG_meshRegularization(qthread->bodydata->qm[index]);
    if (stat != EGADS_SUCCESS)
      printf(" EGADS Warning: EG_fullMeshRegularization %d = %d (EG_quadTess)!\n",
             index+1, stat);
  }

  /* exhausted all work -- exit */
  if (ID != qthread->master) EMP_ThreadExit();
}


int
EG_quadTess(const egObject *tess, egObject **quadTess)
{
  int          i, j, k, m, n, nedges, nfaces, stat, npts, ntris, nside, alen;
  int          outLevel, iv, np, nt, i1, i2, is, ie, ien, oclass, mtype, atype;
  int          sum[2], side[4], degens[2], iuv[2], *table, *tris, *senses;
  const int    *ptype, *pindex, *trs, *trc, *ints;
  double       result[18], xyz[3], uv[2], uvm[2], uvp[2], trange[2], t;
  double       *coords, *parms;
  long         start;
  const double *xyzs, *ts, *uvs, *reals;
  const char   *str;
  egTessel     *btess, *ntess;
  egObject     *obj, *geom, *newTess, **nodes, **edges, **faces, **objs;
  midside      *mid;
  bodyQuad     bodydata;
  EMPquad      qthread;
  void         **threads = NULL;
#ifdef TRIOUT
  FILE         *fp;
  char         filename[100];
#endif
  static int   sides[3][2] = {{1,2}, {2,0}, {0,1}       };
  static int   sideq[4][2] = {{1,2}, {2,5}, {5,0}, {0,1}};
  static int   neigq[4]    = { 0,     3,     4,     2   };

  *quadTess = NULL;
  if (tess == NULL)                 return EGADS_NULLOBJ;
  if (tess->magicnumber != MAGIC)   return EGADS_NOTOBJ;
  if (tess->oclass != TESSELLATION) return EGADS_NOTTESS;
  outLevel = EG_outLevel(tess);

  btess = (egTessel *) tess->blind;
  if (btess == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Blind Object (EG_quadTess)!\n");
    return EGADS_NOTFOUND;
  }
  if (btess->done == 0) {
    if (outLevel > 0)
      printf(" EGADS Error: Tessellation is open (EG_quadTess)!\n");
    return EGADS_TESSTATE;
  }
  obj = btess->src;
  if (obj == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Source Object (EG_quadTess)!\n");
    return EGADS_NULLOBJ;
  }
  if (obj->magicnumber != MAGIC) {
    if (outLevel > 0)
      printf(" EGADS Error: Source Not an Object (EG_quadTess)!\n");
    return EGADS_NOTOBJ;
  }
  if (obj->oclass != BODY) {
    if (outLevel > 0)
      printf(" EGADS Error: Source Not Body (EG_quadTess)!\n");
    return EGADS_NOTBODY;
  }
  if (btess->tess1d == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: No Edge Tessellations (EG_quadTess)!\n");
    return EGADS_NODATA;
  }
  if ((btess->tess2d == NULL) && (btess->nFace != 0)) {
    if (outLevel > 0)
      printf(" EGADS Error: No Face Tessellations (EG_quadTess)!\n");
    return EGADS_NODATA;
  }

  /* initialize the new tessellation object */
  stat = EG_initTessBody(obj, &newTess);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_initTessBody = %d (EG_quadTess)!\n", stat);
    return stat;
  }

  stat = EG_getBodyTopos(obj, NULL, EDGE, &nedges, &edges);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_getBodyTopos E = %d (EG_quadTess)!\n", stat);
    EG_deleteObject(newTess);
    return stat;
  }

  /* rebuild the Edges */

  for (j = i = 0; i < nedges; i++) {
    if (edges[i]->mtype == DEGENERATE) continue;
    stat = EG_getTessEdge(tess, i+1, &npts, &xyzs, &ts);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: EG_getTessEdge %d = %d (EG_quadTess)!\n",
               i+1, stat);
      EG_free(edges);
      EG_deleteObject(newTess);
      return stat;
    }
    if (npts == 0) {
      if (outLevel > 0)
        printf(" EGADS Error: EG_getTessEdge %d -- no points (EG_quadTess)!\n",
               i+1);
      EG_free(edges);
      EG_deleteObject(newTess);
      return EGADS_INDEXERR;
    }
    if (npts > j) j = npts;
  }
  /* allocate to the maximum length */
  coords = (double *) EG_alloc(4*(2*j-1)*sizeof(double));
  if (coords == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Malloc for %d points (EG_quadTess)!\n", j);
    EG_free(edges);
    EG_deleteObject(newTess);
    return EGADS_MALLOC;
  }
  parms = &coords[3*(2*j-1)];

  for (i = 0; i < nedges; i++) {
    if (edges[i]->mtype == DEGENERATE) continue;
    stat = EG_getTessEdge(tess, i+1, &npts, &xyzs, &ts);
    if (stat != EGADS_SUCCESS) continue;

    for (j = 0; j < npts-1; j++) {
      parms[2*j  ] =      ts[j];
      parms[2*j+1] = 0.5*(ts[j] + ts[j+1]);
      stat = EG_evaluate(edges[i], &parms[2*j+1], result);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: EG_evaluate Edge %d/%d = %d (EG_quadTess)!\n",
                 i+1, j+1, stat);
        EG_free(coords);
        EG_free(edges);
        EG_deleteObject(newTess);
        return stat;
      }
      coords[6*j  ] = xyzs[3*j  ];
      coords[6*j+1] = xyzs[3*j+1];
      coords[6*j+2] = xyzs[3*j+2];
      coords[6*j+3] = result[0];
      coords[6*j+4] = result[1];
      coords[6*j+5] = result[2];
    }
    j = npts-1;
    parms[2*j  ]  = ts[j];
    coords[6*j  ] = xyzs[3*j  ];
    coords[6*j+1] = xyzs[3*j+1];
    coords[6*j+2] = xyzs[3*j+2];

    stat = EG_setTessEdge(newTess, i+1, 2*npts-1, coords, parms);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: EG_setTessEdge %d = %d (EG_quadTess)!\n",
               i+1, stat);
      EG_free(coords);
      EG_free(edges);
      EG_deleteObject(newTess);
      return stat;
    }
  }
  EG_free(coords);
  ntess = (egTessel *) newTess->blind;

  /* size and allocate temporary arrays */

  stat = EG_getBodyTopos(obj, NULL, FACE, &nfaces, &faces);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_getBodyTopos F = %d (EG_quadTess)!\n", stat);
    EG_free(edges);
    EG_deleteObject(newTess);
    return stat;
  }
  for (nside = npts = ntris = i = 0; i < nfaces; i++) {
    stat = EG_getTessFace(tess, i+1, &np, &xyzs, &uvs, &ptype, &pindex,
                          &nt, &trs, &trc);
    if ((stat != EGADS_SUCCESS) || (nt == 0)) {
      if (outLevel > 0)
        if (stat != EGADS_SUCCESS) {
          printf(" EGADS Error: EG_getTessFace %d = %d (EG_quadTess)!\n",
                 i+1, stat);
        } else {
          printf(" EGADS Error: Face %d has no tessellation (EG_quadTess)!\n",
                 i+1);
        }
      EG_free(faces);
      EG_free(edges);
      EG_deleteObject(newTess);
      return stat;
    }
    for (sum[0] = sum[1] = j = 0; j < nt; j++)
      for (k = 0; k < 3; k++)
        if (trc[3*j+k] > 0) {
          sum[0]++;
        } else {
          sum[1]++;
        }
    k = sum[0]/2 + sum[1];
    if (k  > nside) nside = k;
    if (np > npts)  npts  = np;
    if (nt > ntris) ntris = nt;
  }
  k      = npts + ntris + nside;
  coords = (double *) EG_alloc(5*k*sizeof(double));
  if (coords == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Malloc for %d Face points (EG_quadTess)!\n", k);
    EG_free(faces);
    EG_free(edges);
    EG_deleteObject(newTess);
    return EGADS_MALLOC;
  }
  parms = &coords[3*k];
  tris  = (int *) EG_alloc(6*3*ntris*sizeof(int));
  if (tris == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Malloc for %d Face tris (EG_quadTess)!\n", 3*ntris);
    EG_free(coords);
    EG_free(faces);
    EG_free(edges);
    EG_deleteObject(newTess);
    return EGADS_MALLOC;
  }
  mid = (midside *) EG_alloc(nside*sizeof(midside));
  if (mid == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Malloc for %d Face mids (EG_quadTess)!\n", nside);
    EG_free(tris);
    EG_free(coords);
    EG_free(faces);
    EG_free(edges);
    EG_deleteObject(newTess);
    return EGADS_MALLOC;
  }
  table = (int *) EG_alloc(npts*sizeof(int));
  if (table == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Malloc for %d Face table (EG_quadTess)!\n", npts);
    EG_free(mid);
    EG_free(tris);
    EG_free(coords);
    EG_free(faces);
    EG_free(edges);
    EG_deleteObject(newTess);
    return EGADS_MALLOC;
  }

  /* fill in the quads (as pairs of tris) a Face at a time */

  for (i = 0; i < nfaces; i++) {
    if (btess->tess2d[i].tfi == 1) continue;
    stat = EG_getTessFace(tess, i+1, &np, &xyzs, &uvs, &ptype, &pindex,
                          &nt, &trs, &trc);
    if (stat != EGADS_SUCCESS) continue;

    /* find degenerate nodes (if any) */
    degens[0] = degens[1] = 0;
    iuv[0]    = iuv[1]    = 0;
    stat      = EG_getBodyTopos(obj, faces[i], EDGE, &k, &objs);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Internal: EG_getBodyTopos on Face %d = %d\n", i+1, stat);
    } else {
      for (j = 0; j < k; j++) {
        stat = EG_getTopology(objs[j], &geom, &oclass, &mtype,
                              trange, &n, &nodes, &senses);
        if (stat != EGADS_SUCCESS) {
          printf(" EGADS Internal: EG_getTopology on Edge = %d\n", stat);
          continue;
        }
        if (mtype != DEGENERATE) continue;
        stat = EG_getEdgeUVeval(faces[i], objs[j], 0, trange[0], result);
        if (stat != EGADS_SUCCESS) {
          printf(" EGADS Internal: EG_getEdgeUVeval = %d\n", stat);
          continue;
        }
        n = EG_indexBodyTopo(obj, nodes[0]);
        if (n > 0) {
          if (degens[0] == 0) {
            degens[0] = n;
            if (result[3] != 0.0) iuv[0] = 1;
          } else if (degens[1] == 0) {
            degens[1] = n;
            if (result[3] != 0.0) iuv[1] = 1;
          } else {
            printf(" EGADS Info: More than 2 Degen Nodes in Face %d!\n", i+1);
          }
        }
      }
      EG_free(objs);
    }

    /* make the vertices */
    for (j = 0; j < np; j++) {
      coords[3*j  ] = xyzs[3*j  ];
      coords[3*j+1] = xyzs[3*j+1];
      coords[3*j+2] = xyzs[3*j+2];
      parms[2*j  ]  = uvs[2*j  ];
      parms[2*j+1]  = uvs[2*j+1];
      table[j]      = NOTFILLED;
    }
    iv = np;

    /* get triangle side midpoint insertions */
    for (is = j = 0; j < nt; j++)
      for (k = 0; k < 3; k++) {
        i1 = trs[3*j+sides[k][0]] - 1;
        i2 = trs[3*j+sides[k][1]] - 1;
        if (i2 < i1) {
          stat = i1;
          i1   = i2;
          i2   = stat;
        }
        m = EG_findMidSide(i1, i2, table, mid);
        if (m > 0) continue;
        uvm[0] = uvs[2*i1  ];
        uvm[1] = uvs[2*i1+1];
        uvp[0] = uvs[2*i2  ];
        uvp[1] = uvs[2*i2+1];
        if (trc[3*j+k] > 0) {
          /* Interior side */
          if ((ptype[i1] == 0) && (pindex[i1] == degens[0])) {
            if (iuv[0] == 0) {
              uvm[0] = uvs[2*i2  ];
            } else {
              uvm[1] = uvs[2*i2+1];
            }
          } else if ((ptype[i1] == 0) && (pindex[i1] == degens[1])) {
            if (iuv[1] == 0) {
              uvm[0] = uvs[2*i2  ];
            } else {
              uvm[1] = uvs[2*i2+1];
            }
          } else if ((ptype[i2] == 0) && (pindex[i2] == degens[0])) {
            if (iuv[0] == 0) {
              uvp[0] = uvs[2*i1  ];
            } else {
              uvp[1] = uvs[2*i1+1];
            }
          } else if ((ptype[i2] == 0) && (pindex[i2] == degens[1])) {
            if (iuv[1] == 0) {
              uvp[0] = uvs[2*i1  ];
            } else {
              uvp[1] = uvs[2*i1+1];
            }
          }
          EG_getSidepoint(faces[i], 0.5, uvm, uvp, uv);
          stat = EG_evaluate(faces[i], uv, result);
        } else {
          /* Edge side */
          ie   = -trc[3*j+k] - 1;
          stat = EG_getTopology(edges[ie], &geom, &oclass, &mtype,
                                trange, &n, &objs, &senses);
          if (stat != EGADS_SUCCESS) {
            if (outLevel > 0)
              printf(" EGADS Error: EG_getTopology %d/%d = %d (EG_quadTess)!\n",
                     i+1, ie+1, stat);
            EG_free(table);
            EG_free(mid);
            EG_free(tris);
            EG_free(coords);
            EG_free(faces);
            EG_free(edges);
            EG_deleteObject(newTess);
            return stat;
          }
          n = EG_indexBodyTopo(obj, objs[0]);
          if ((ptype[i1] == 0) && (ptype[i2] == 0)) {
            ien = 1;
          } else if (ptype[i1] == 0) {
            if (mtype == ONENODE) {
              if ((xyzs[3*i2  ] == ntess->tess1d[ie].xyz[6]) &&
                  (xyzs[3*i2+1] == ntess->tess1d[ie].xyz[7]) &&
                  (xyzs[3*i2+2] == ntess->tess1d[ie].xyz[8])) {
                ien = 1;
              } else {
                ien = ntess->tess1d[ie].npts - 2;
              }
            } else {
              if (pindex[i1] == n) {
                ien = 1;
              } else {
                ien = 2*ptype[i2] - 1;
              }
            }
            if (ie+1 != pindex[i2])
              printf(" EGADS Info: Edge mismatch %d %d\n", ie+1, pindex[i2]);
          } else if (ptype[i2] == 0) {
            if (mtype == ONENODE) {
              if ((xyzs[3*i1  ] == ntess->tess1d[ie].xyz[6]) &&
                  (xyzs[3*i1+1] == ntess->tess1d[ie].xyz[7]) &&
                  (xyzs[3*i1+2] == ntess->tess1d[ie].xyz[8])) {
                ien = 1;
              } else {
                ien = ntess->tess1d[ie].npts - 2;
              }
            } else {
              if (pindex[i2] == n) {
                ien = 1;
              } else {
                ien = 2*ptype[i1] - 1;
              }
            }
            if (ie+1 != pindex[i1])
              printf(" EGADS Info: Edge mismatch %d %d\n", ie+1, pindex[i1]);
          } else {
            ien = 2*ptype[i1] - 1;
            if (ien > 2*ptype[i2]-1) ien = 2*ptype[i2] - 1;
            if (ien < 0) {
              printf(" EGADS Internal: EG_quadTess ien = %d!\n", ien);
              continue;
            }
            if ((ie+1 != pindex[i1]) || (ie+1 != pindex[i2]))
              printf(" EGADS Info: Edge mismatch %d %d %d\n", ie+1, pindex[i1],
                     pindex[i2]);
          }
          t    = ntess->tess1d[ie].t[ien];
          stat = EG_getEdgeUV(faces[i], edges[ie], 0, t, uv);
          if (stat == EGADS_TOPOERR) {
            /* sense in Face twice! */
            stat = EG_getEdgeUV(faces[i], edges[ie], -1, t, uvm);
            if (stat != EGADS_SUCCESS) {
              printf(" EGADS Info:  EG_getEdgeUV -1\n");
            } else {
              stat = EG_getEdgeUV(faces[i], edges[ie],  1, t, uvp);
              if (stat != EGADS_SUCCESS) {
                printf(" EGADS Info:  EG_getEdgeUV +1\n");
              } else {
                result[0] = sqrt((uvm[0]-uv[0])*(uvm[0]-uv[0]) +
                                 (uvm[1]-uv[1])*(uvm[1]-uv[1]));
                result[1] = sqrt((uvp[0]-uv[0])*(uvp[0]-uv[0]) +
                                 (uvp[1]-uv[1])*(uvp[1]-uv[1]));
                if (result[0] < result[1]) {
                  uv[0] = uvm[0];
                  uv[1] = uvm[1];
                } else {
                  uv[0] = uvp[0];
                  uv[1] = uvp[1];
                }
              }
            }
          }
          if (stat != EGADS_SUCCESS) {
            if (outLevel > 0)
              printf(" EGADS Error: EG_getEdgeUV %d/%d = %d (EG_quadTess)!\n",
                     i+1, ie+1, stat);
            EG_free(table);
            EG_free(mid);
            EG_free(tris);
            EG_free(coords);
            EG_free(faces);
            EG_free(edges);
            EG_deleteObject(newTess);
            return stat;
          }
          result[0] = ntess->tess1d[ie].xyz[3*ien  ];
          result[1] = ntess->tess1d[ie].xyz[3*ien+1];
          result[2] = ntess->tess1d[ie].xyz[3*ien+2];
        }
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: EG_evaluate %d %d/%d = %d (EG_quadTess)!\n",
                   i+1, j+1, k+1, stat);
          EG_free(table);
          EG_free(mid);
          EG_free(tris);
          EG_free(coords);
          EG_free(faces);
          EG_free(edges);
          EG_deleteObject(newTess);
          return stat;
        }
        if (m == 0) {
          table[i1]      = is;
        } else {
          mid[-m-1].next = is;
        }
        mid[is].vert2 = i2;
        mid[is].nvert = iv + 1;
        mid[is].next  = NOTFILLED;
        is++;
        coords[3*iv  ] = result[0];
        coords[3*iv+1] = result[1];
        coords[3*iv+2] = result[2];
        parms[2*iv  ]  = uv[0];
        parms[2*iv+1]  = uv[1];
        iv++;
      }

    /* fill in the new triangulation based on Quads */
    for (j = 0; j < nt; j++) {
      for (k = 0; k < 3; k++) {
        i1 = trs[3*j+sides[k][0]] - 1;
        i2 = trs[3*j+sides[k][1]] - 1;
        if (i2 < i1) {
          stat = i1;
          i1   = i2;
          i2   = stat;
        }
        side[k] = EG_findMidSide(i1, i2, table, mid);
        if (side[k] <= 0) {
          if (outLevel > 0)
            printf(" EGADS Error: findMidSide %d = %d (EG_quadTess)!\n",
                   i+1, side[k]);
          EG_free(table);
          EG_free(mid);
          EG_free(tris);
          EG_free(coords);
          EG_free(faces);
          EG_free(edges);
          EG_deleteObject(newTess);
          return stat;
        }
      }
      /* get middle vertex based on mid-sides */
      xyz[0]  = (coords[3*(side[0]-1)  ] + coords[3*(side[1]-1)  ] +
                 coords[3*(side[2]-1)  ])/3.0;
      xyz[1]  = (coords[3*(side[0]-1)+1] + coords[3*(side[1]-1)+1] +
                 coords[3*(side[2]-1)+1])/3.0;
      xyz[2]  = (coords[3*(side[0]-1)+2] + coords[3*(side[1]-1)+2] +
                 coords[3*(side[2]-1)+2])/3.0;
      uv[0]   = (parms[2*(side[0]-1)  ]  + parms[2*(side[1]-1)  ] +
                 parms[2*(side[2]-1)  ])/3.0;
      uv[1]   = (parms[2*(side[0]-1)+1]  + parms[2*(side[1]-1)+1] +
                 parms[2*(side[2]-1)+1])/3.0;
      EG_getInterior(faces[i], xyz, uv);
      stat  = EG_evaluate(faces[i], uv, result);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: EG_evaluate %d %d = %d (EG_quadTess)!\n",
                 i+1, j+1, stat);
        EG_free(table);
        EG_free(mid);
        EG_free(tris);
        EG_free(coords);
        EG_free(faces);
        EG_free(edges);
        EG_deleteObject(newTess);
        return stat;
      }
      coords[3*iv  ] = result[0];
      coords[3*iv+1] = result[1];
      coords[3*iv+2] = result[2];
      parms[2*iv  ]  = uv[0];
      parms[2*iv+1]  = uv[1];
      iv++;

      tris[18*j   ] = trs[3*j  ];
      tris[18*j+ 1] = side[2];
      tris[18*j+ 2] = iv;
      tris[18*j+ 3] = trs[3*j  ];
      tris[18*j+ 4] = iv;
      tris[18*j+ 5] = side[1];

      tris[18*j+ 6] = trs[3*j+1];
      tris[18*j+ 7] = side[0];
      tris[18*j+ 8] = iv;
      tris[18*j+ 9] = trs[3*j+1];
      tris[18*j+10] = iv;
      tris[18*j+11] = side[2];

      tris[18*j+12] = trs[3*j+2];
      tris[18*j+13] = side[1];
      tris[18*j+14] = iv;
      tris[18*j+15] = trs[3*j+2];
      tris[18*j+16] = iv;
      tris[18*j+17] = side[0];
    }

    stat = EG_setTessFace(newTess, i+1, iv, coords, parms, 6*nt, tris);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: EG_setTessFace %d = %d (EG_quadTess)!\n",
               i+1, stat);
      EG_free(table);
      EG_free(mid);
      EG_free(tris);
      EG_free(coords);
      EG_free(faces);
      EG_free(edges);
      EG_deleteObject(newTess);
      return stat;
    }
#ifdef TRIOUT
    sprintf(filename, "Components.%d.i.tri", i+1);
    fp = fopen(filename, "w");
    fprintf(fp," %d %d\n", iv, 6*nt);
    for (j = 0; j < iv; j++)
      fprintf(fp, " %lf %lf %lf\n",coords[3*j  ],coords[3*j+1],coords[3*j+2]);
    for (j = 0; j < 6*nt; j++)
      fprintf(fp, " %d %d %d\n", tris[3*j  ], tris[3*j+1], tris[3*j+2]);
    for (j = 0; j < 6*nt; j++)
      fprintf(fp, " 1\n");
    fclose(fp);
#endif
  }

  /* fill in the quads from TFI Faces */

  for (i = 0; i < nfaces; i++) {
    if (btess->tess2d[i].tfi != 1) continue;
    stat = EG_getTessFace(tess, i+1, &np, &xyzs, &uvs, &ptype, &pindex,
                          &nt, &trs, &trc);
    if (stat != EGADS_SUCCESS) continue;

    /* find degenerate nodes (if any) */
    degens[0] = degens[1] = 0;
    iuv[0]    = iuv[1]    = 0;
    stat = EG_getBodyTopos(obj, faces[i], EDGE, &k, &objs);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Internal: EG_getBodyTopos on Face %d = %d\n", i+1, stat);
    } else {
      for (j = 0; j < k; j++) {
        stat = EG_getTopology(objs[j], &geom, &oclass, &mtype,
                              trange, &n, &nodes, &senses);
        if (stat != EGADS_SUCCESS) {
          printf(" EGADS Internal: EG_getTopology on Edge = %d\n", stat);
          continue;
        }
        if (mtype != DEGENERATE) continue;
        stat = EG_getEdgeUVeval(faces[i], objs[j], 0, trange[0], result);
        if (stat != EGADS_SUCCESS) {
          printf(" EGADS Internal: EG_getEdgeUVeval = %d\n", stat);
          continue;
        }
        n = EG_indexBodyTopo(obj, nodes[0]);
        if (n > 0) {
          if (degens[0] == 0) {
            degens[0] = n;
            if (result[3] != 0.0) iuv[0] = 1;
          } else if (degens[1] == 0) {
            degens[1] = n;
            if (result[3] != 0.0) iuv[1] = 1;
          } else {
            printf(" EGADS Info: More than 2 Degen Nodes in Face %d!\n",
                   i+1);
          }
        }
      }
      EG_free(objs);
    }

    /* make the vertices */
    for (j = 0; j < np; j++) {
      coords[3*j  ] = xyzs[3*j  ];
      coords[3*j+1] = xyzs[3*j+1];
      coords[3*j+2] = xyzs[3*j+2];
      parms[2*j  ]  = uvs[2*j  ];
      parms[2*j+1]  = uvs[2*j+1];
      table[j]      = NOTFILLED;
    }
    iv = np;

    /* get quad side midpoint insertions */
    for (is = j = 0; j < nt/2; j++)
      for (k = 0; k < 4; k++) {
        i1 = trs[6*j+sideq[k][0]] - 1;
        i2 = trs[6*j+sideq[k][1]] - 1;
        if (i2 < i1) {
          stat = i1;
          i1   = i2;
          i2   = stat;
        }
        m = EG_findMidSide(i1, i2, table, mid);
        if (m > 0) continue;
        uv[0] = 0.5*(uvs[2*i1  ] + uvs[2*i2  ]);
        uv[1] = 0.5*(uvs[2*i1+1] + uvs[2*i2+1]);
        ie = trc[6*j+neigq[k]];
        if (ie > 0) {
          /* Interior side */
          if ((ptype[i1] == 0) && (pindex[i1] == degens[0])) {
            if (iuv[0] == 0) {
              uv[0] = uvs[2*i2  ];
            } else {
              uv[1] = uvs[2*i2+1];
            }
          } else if ((ptype[i1] == 0) && (pindex[i1] == degens[1])) {
            if (iuv[1] == 0) {
              uv[0] = uvs[2*i2  ];
            } else {
              uv[1] = uvs[2*i2+1];
            }
          } else if ((ptype[i2] == 0) && (pindex[i2] == degens[0])) {
            if (iuv[0] == 0) {
              uv[0] = uvs[2*i1  ];
            } else {
              uv[1] = uvs[2*i1+1];
            }
          } else if ((ptype[i2] == 0) && (pindex[i2] == degens[1])) {
            if (iuv[1] == 0) {
              uv[0] = uvs[2*i1  ];
            } else {
              uv[1] = uvs[2*i1+1];
            }
          }
          stat = EG_evaluate(faces[i], uv, result);
        } else {
          /* Edge side */
          ie   = -ie - 1;
          stat = EG_getTopology(edges[ie], &geom, &oclass, &mtype,
                                trange, &n, &objs, &senses);
          if (stat != EGADS_SUCCESS) {
            if (outLevel > 0)
              printf(" EGADS Error: EG_getTopology %d/%d = %d (EG_quadTess)!\n",
                     i+1, ie+1, stat);
            EG_free(table);
            EG_free(mid);
            EG_free(tris);
            EG_free(coords);
            EG_free(faces);
            EG_free(edges);
            EG_deleteObject(newTess);
            return stat;
          }
          n = EG_indexBodyTopo(obj, objs[0]);
          if ((ptype[i1] == 0) && (ptype[i2] == 0)) {
            ien = 1;
          } else if (ptype[i1] == 0) {
            if (pindex[i1] == n) {
              ien = 1;
            } else {
              ien = 2*ptype[i2] - 1;
            }
            if (ie+1 != pindex[i2])
              printf(" EGADS Info: Edge mismatch %d %d\n", ie+1, pindex[i2]);
          } else if (ptype[i2] == 0) {
            if (pindex[i2] == n) {
              ien = 1;
            } else {
              ien = 2*ptype[i1] - 1;
            }
            if (ie+1 != pindex[i1])
              printf(" EGADS Info: Edge mismatch %d %d\n", ie+1, pindex[i1]);
          } else {
            ien = 2*ptype[i1] - 1;
            if (ien > 2*ptype[i2]-1) ien = 2*ptype[i2] - 1;
            if (ien < 0) {
              printf(" EGADS Internal: EG_quadTess ien = %d!\n", ien);
              continue;
            }
            if ((ie+1 != pindex[i1]) || (ie+1 != pindex[i2]))
              printf(" EGADS Info: Edge mismatch %d %d %d\n", ie+1, pindex[i1],
                     pindex[i2]);
          }
          t    = ntess->tess1d[ie].t[ien];
          stat = EG_getEdgeUV(faces[i], edges[ie], 0, t, uv);
          if (stat == EGADS_TOPOERR) {
            /* sense in Face twice! */
            stat = EG_getEdgeUV(faces[i], edges[ie], -1, t, uvm);
            if (stat != EGADS_SUCCESS) {
              printf(" EGADS Info:  EG_getEdgeUV -1\n");
            } else {
              stat = EG_getEdgeUV(faces[i], edges[ie],  1, t, uvp);
              if (stat != EGADS_SUCCESS) {
                printf(" EGADS Info:  EG_getEdgeUV +1\n");
              } else {
                result[0] = sqrt((uvm[0]-uv[0])*(uvm[0]-uv[0]) +
                                 (uvm[1]-uv[1])*(uvm[1]-uv[1]));
                result[1] = sqrt((uvp[0]-uv[0])*(uvp[0]-uv[0]) +
                                 (uvp[1]-uv[1])*(uvp[1]-uv[1]));
                if (result[0] < result[1]) {
                  uv[0] = uvm[0];
                  uv[1] = uvm[1];
                } else {
                  uv[0] = uvp[0];
                  uv[1] = uvp[1];
                }
              }
            }
          }
          if (stat != EGADS_SUCCESS) {
            if (outLevel > 0)
              printf(" EGADS Error: EG_getEdgeUV %d/%d = %d (EG_quadTess)!\n",
                     i+1, ie+1, stat);
            EG_free(table);
            EG_free(mid);
            EG_free(tris);
            EG_free(coords);
            EG_free(faces);
            EG_free(edges);
            EG_deleteObject(newTess);
            return stat;
          }
          result[0] = ntess->tess1d[ie].xyz[3*ien  ];
          result[1] = ntess->tess1d[ie].xyz[3*ien+1];
          result[2] = ntess->tess1d[ie].xyz[3*ien+2];
        }
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: EG_evaluate %d %d/%d = %d (EG_quadTess)!\n",
                   i+1, j+1, k+1, stat);
          EG_free(table);
          EG_free(mid);
          EG_free(tris);
          EG_free(coords);
          EG_free(faces);
          EG_free(edges);
          EG_deleteObject(newTess);
          return stat;
        }
        if (m == 0) {
          table[i1]      = is;
        } else {
          mid[-m-1].next = is;
        }
        mid[is].vert2 = i2;
        mid[is].nvert = iv + 1;
        mid[is].next  = NOTFILLED;
        is++;
        coords[3*iv  ] = result[0];
        coords[3*iv+1] = result[1];
        coords[3*iv+2] = result[2];
        parms[2*iv  ]  = uv[0];
        parms[2*iv+1]  = uv[1];
        iv++;
      }

    /* fill in the new triangulation based on Quads */
    for (j = 0; j < nt/2; j++) {
      for (k = 0; k < 4; k++) {
        i1 = trs[6*j+sideq[k][0]] - 1;
        i2 = trs[6*j+sideq[k][1]] - 1;
        if (i2 < i1) {
          stat = i1;
          i1   = i2;
          i2   = stat;
        }
        side[k] = EG_findMidSide(i1, i2, table, mid);
        if (side[k] <= 0) {
          if (outLevel > 0)
            printf(" EGADS Error: findMidSide %d = %d (EG_quadTess)!\n",
                   i+1, side[k]);
          EG_free(table);
          EG_free(mid);
          EG_free(tris);
          EG_free(coords);
          EG_free(faces);
          EG_free(edges);
          EG_deleteObject(newTess);
          return stat;
        }
      }
      /* get middle vertex based on mid-sides */
      uv[0] = 0.25*(parms[2*(side[0]-1)  ] + parms[2*(side[1]-1)  ] +
                    parms[2*(side[2]-1)  ] + parms[2*(side[3]-1)  ]);
      uv[1] = 0.25*(parms[2*(side[0]-1)+1] + parms[2*(side[1]-1)+1] +
                    parms[2*(side[2]-1)+1] + parms[2*(side[3]-1)+1]);
      stat  = EG_evaluate(faces[i], uv, result);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: EG_evaluate %d %d = %d (EG_quadTess)!\n",
                 i+1, j+1, stat);
        EG_free(table);
        EG_free(mid);
        EG_free(tris);
        EG_free(coords);
        EG_free(faces);
        EG_free(edges);
        EG_deleteObject(newTess);
        return stat;
      }
      coords[3*iv  ] = result[0];
      coords[3*iv+1] = result[1];
      coords[3*iv+2] = result[2];
      parms[2*iv  ]  = uv[0];
      parms[2*iv+1]  = uv[1];
      iv++;

      tris[24*j   ] = trs[6*j  ];
      tris[24*j+ 1] = side[3];
      tris[24*j+ 2] = iv;
      tris[24*j+ 3] = trs[6*j  ];
      tris[24*j+ 4] = iv;
      tris[24*j+ 5] = side[2];

      tris[24*j+ 6] = trs[6*j+1];
      tris[24*j+ 7] = side[0];
      tris[24*j+ 8] = iv;
      tris[24*j+ 9] = trs[6*j+1];
      tris[24*j+10] = iv;
      tris[24*j+11] = side[3];

      tris[24*j+12] = trs[6*j+2];
      tris[24*j+13] = side[1];
      tris[24*j+14] = iv;
      tris[24*j+15] = trs[6*j+2];
      tris[24*j+16] = iv;
      tris[24*j+17] = side[0];

      tris[24*j+18] = trs[6*j+5];
      tris[24*j+19] = side[2];
      tris[24*j+20] = iv;
      tris[24*j+21] = trs[6*j+5];
      tris[24*j+22] = iv;
      tris[24*j+23] = side[1];
    }

    stat = EG_setTessFace(newTess, i+1, iv, coords, parms, 4*nt, tris);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: EG_setTessFace %d = %d (EG_quadTess)!\n",
               i+1, stat);
      EG_free(table);
      EG_free(mid);
      EG_free(tris);
      EG_free(coords);
      EG_free(faces);
      EG_free(edges);
      EG_deleteObject(newTess);
      return stat;
    }
    /* set tfi flag in new tessellation */
    ntess->tess2d[i].tfi = 1;
#ifdef TRIOUT
    sprintf(filename, "Components.%d.i.tri", i+1);
    fp = fopen(filename, "w");
    fprintf(fp," %d %d\n", iv, 4*nt);
    for (j = 0; j < iv; j++)
      fprintf(fp, " %lf %lf %lf\n",coords[3*j  ],coords[3*j+1],coords[3*j+2]);
    for (j = 0; j < 4*nt; j++)
      fprintf(fp, " %d %d %d\n", tris[3*j  ], tris[3*j+1], tris[3*j+2]);
    for (j = 0; j < 4*nt; j++)
      fprintf(fp, " 1\n");
    fclose(fp);
#endif
  }

  /* clean up temp storage */
  EG_free(table);
  EG_free(mid);
  EG_free(tris);
  EG_free(coords);
  EG_free(faces);
  EG_free(edges);

  /* close up the open tessellation */
  stat = EG_statusTessBody(newTess, &obj, &i, &npts);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_statusTessBody = %d (EG_quadTess)!\n", stat);
    EG_deleteObject(newTess);
    return stat;
  }
  if (i != 1) {
    if (outLevel > 0)
      printf(" EGADS Error: New Tessellation Object is Open (EG_quadTess)!\n");
    EG_deleteObject(newTess);
    return EGADS_TESSTATE;
  }
#ifndef LITE
  table = (int *) EG_alloc(nfaces*sizeof(int));
  if (table != NULL) {
    for (i = 0; i < nfaces; i++) table[i] = ntess->tess2d[i].ntris/2;
    stat = EG_attributeAdd(newTess, ".mixed", ATTRINT, nfaces,
                           table, NULL, NULL);
    if (stat != EGADS_SUCCESS)
      printf(" EGADS Warning: EG_attributeAdd m = %d !\n", stat);
    EG_free(table);
  }
  stat = EG_attributeAdd(newTess, ".tessType", ATTRSTRING, 4,
                         NULL, NULL, "Quad");
  if (stat != EGADS_SUCCESS)
    if (outLevel > 0)
      printf(" EGADS Warning: EG_attributeAdd Q = %d (EG_quadTess)!\n", stat);
#endif

  /* do we regularize? */
  stat = EG_attributeRet(tess, ".qRegular", &atype, &alen, &ints, &reals, &str);
  if (stat == EGADS_SUCCESS)
    if ((atype == ATTRSTRING) && (str != NULL))
      if ((strcmp(str, "off") == 0) || (strcmp(str, "Off") == 0) ||
          (strcmp(str, "OFF") == 0)) {
        *quadTess = newTess;
        return EGADS_SUCCESS;
      }
#ifndef REGULAR
  if (stat == EGADS_NOTFOUND) {
    *quadTess = newTess;
    return EGADS_SUCCESS;
  }
#endif

  /* yes! */
  bodydata.tess = newTess;
  stat = EG_getBodyTopos(btess->src, NULL, EDGE, &bodydata.nedges, NULL);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Warning: EG_getBodyTopos E = %d (EG_quadTess)!\n", stat);
    EG_deleteObject(newTess);
    return stat;
  }
  stat = EG_getBodyTopos(btess->src, NULL, FACE, &bodydata.nfaces,
                         &bodydata.faces);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Warning: EG_getBodyTopos F = %d (EG_quadTess)!\n", stat);
    EG_deleteObject(newTess);
    return stat;
  }
  stat = EG_createMeshMap(&bodydata);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Warning: EG_createMeshMap = %d (EG_quadTess)!\n", stat);
    EG_free(bodydata.faces);
    *quadTess = newTess;
    return EGADS_SUCCESS;
  }

  /* set the thread storage */
  qthread.mutex    = NULL;
  qthread.master   = EMP_ThreadID();
  qthread.end      = bodydata.nfaces;
  qthread.index    = 0;
  qthread.ntess    = ntess;
  qthread.bodydata = &bodydata;

  np = EMP_Init(&start);
  if (outLevel > 1) printf(" EMP NumProcs = %d!\n", np);

  if (np > 1) {
    /* create the mutex to handle list synchronization */
    qthread.mutex = EMP_LockCreate();
    if (qthread.mutex == NULL) {
      printf(" EMP Error: mutex creation = NULL!\n");
      np = 1;
    } else {
      /* get storage for our extra threads */
      threads = (void **) malloc((np-1)*sizeof(void *));
      if (threads == NULL) {
        EMP_LockDestroy(qthread.mutex);
        np = 1;
      }
    }
  }

  /* create the threads and get going! */
  if (threads != NULL)
    for (i = 0; i < np-1; i++) {
      threads[i] = EMP_ThreadCreate(EG_quadThread, &qthread);
      if (threads[i] == NULL)
        printf(" EMP Error Creating Thread #%d!\n", i+1);
    }
  /* now run the thread block from the original thread */
  EG_quadThread(&qthread);

  /* wait for all others to return */
  if (threads != NULL)
    for (i = 0; i < np-1; i++)
      if (threads[i] != NULL) EMP_ThreadWait(threads[i]);

  /* thread cleanup */
  if (threads != NULL)
    for (i = 0; i < np-1; i++)
      if (threads[i] != NULL) EMP_ThreadDestroy(threads[i]);
  if (qthread.mutex != NULL) EMP_LockDestroy(qthread.mutex);
  if (threads != NULL) free(threads);
  if (outLevel > 1)
    printf(" EMP Number of Seconds on Quad Thread Block = %ld\n",
           EMP_Done(&start));

  /* collect all of the Face quads and make the tessellation object */
  stat = EG_makeQuadTess(bodydata, &obj);
  EG_deleteObject(newTess);
  EG_destroyMeshMap(&bodydata);
  EG_free(bodydata.faces);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Warning: EG_makeQuadTess = %d (EG_quadTess)!\n", stat);
    return stat;
  }

  *quadTess = obj;
  return EGADS_SUCCESS;
}
