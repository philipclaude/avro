// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef SCALARVECTOR_H
#define SCALARVECTOR_H

#include "SparseVector.h"
#include "numpack/block/block_Type.h"

namespace numpack 
{
namespace SLA
{

struct ScalarVector
{
  template<class VectorS> ScalarVector( const SparseVector<VectorS>& x );
  template<class VectorS> void setTo( SparseVector<VectorS>& x);

  template<class VectorS> ScalarVector( const SparseVector< DLA::VectorD<VectorS> >& x );
  template<class VectorS> void setTo( SparseVector< DLA::VectorD<VectorS> >& x);

  template<class VectorS> ScalarVector( const DLA::MatrixDView< SparseVector<VectorS> >& x );
  template<class VectorS> void setTo( DLA::MatrixDView< SparseVector<VectorS> >& x);

  template<class VectorS> ScalarVector( const DLA::MatrixDView< SparseVector< DLA::VectorD<VectorS> > >& x );
  template<class VectorS> void setTo( DLA::MatrixDView< SparseVector< DLA::VectorD<VectorS> > >& x);

  template<class Vector0, class Vector1>
  ScalarVector( const BLA::VectorBlock_2< DLA::VectorD<SparseVector<Vector0> >,
                                          DLA::VectorD<SparseVector<Vector1> > >& x );
  template<class Vector0, class Vector1>
  void setTo( BLA::VectorBlock_2< DLA::VectorD<SparseVector<Vector0> >,
                                  DLA::VectorD<SparseVector<Vector1> > >& x );

  template<class Vector0, class Vector1, class Vector2>
  ScalarVector( const BLA::VectorBlock_3< DLA::VectorD<SparseVector<Vector0> >,
                                          DLA::VectorD<SparseVector<Vector1> >,
                                          DLA::VectorD<SparseVector<Vector2> > >& x );
  template<class Vector0, class Vector1, class Vector2>
  void setTo( BLA::VectorBlock_3< DLA::VectorD<SparseVector<Vector0> >,
                                  DLA::VectorD<SparseVector<Vector1> >,
                                  DLA::VectorD<SparseVector<Vector2> > >& x );

  template<class Vector0, class Vector1, class Vector2, class Vector3>
  ScalarVector( const BLA::VectorBlock_4< DLA::VectorD<SparseVector<Vector0> >,
                                          DLA::VectorD<SparseVector<Vector1> >,
                                          DLA::VectorD<SparseVector<Vector2> >,
                                          DLA::VectorD<SparseVector<Vector3> > >& x );
  template<class Vector0, class Vector1, class Vector2, class Vector3>
  void setTo( BLA::VectorBlock_4< DLA::VectorD<SparseVector<Vector0> >,
                                  DLA::VectorD<SparseVector<Vector1> >,
                                  DLA::VectorD<SparseVector<Vector2> >,
                                  DLA::VectorD<SparseVector<Vector3> > >& x );



  ScalarVector(const ScalarVector&) = delete;

  operator double*() { return v; }

  ~ScalarVector() { delete [] v; }

  int m;
  double *v = NULL;
};

} // namespace SLA
} // namespace numpack 

#endif //SCALARVECTOR_H
