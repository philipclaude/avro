// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

// No include blocker here. This file should only be included once in a cpp file for explicit instantiation

#if !defined(SCALARVECTOR_INSTANTIATE) && !defined(SANS_HEADERCOMPILE_CHECK)
#error "This file should only be included in a cpp file for explicit instantiations"
#endif

#include "ScalarVector.h"
#include "LinearAlgebra/DenseLinAlg/StaticSize/VectorS.h"
#include "LinearAlgebra/DenseLinAlg/tools/index.h"
#include "LinearAlgebra/DenseLinAlg/tools/VectorSize.h"

namespace SANS
{
namespace SLA
{

//---------------------------------------------------------------------------//
template<class VectorS>
ScalarVector::ScalarVector( const SparseVector<VectorS>& x )
{
  // Deduce the sizes of VectorS type
  const int M = DLA::VectorSize<VectorS>::M;

  m = x.m()*M;
  v = new double[m];

  int n = 0;
  //Loop over the rows of blocks
  for ( int bi = 0; bi < x.m(); bi++ )
    //Loop over the rows in a fixed block
    for ( int fi = 0; fi < M; fi++ )
      v[n++] = DLA::index(x[bi],fi);
}

template<class VectorS>
void
ScalarVector::setTo( SparseVector<VectorS>& x)
{
  // Deduce the sizes of VectorS type
  const int M = DLA::VectorSize<VectorS>::M;

  int n = 0;
  //Loop over the rows of blocks
  for ( int bi = 0; bi < x.m(); bi++ )
    //Loop over the rows in a fixed block
    for ( int fi = 0; fi < M; fi++ )
      DLA::index(x[bi],fi) = v[n++];
}

//---------------------------------------------------------------------------//
template<class VectorS>
ScalarVector::ScalarVector( const SparseVector< DLA::VectorD<VectorS> >& x )
{
  // Deduce the sizes of VectorS type
  const int M = DLA::VectorSize<VectorS>::M;

  m = x.value_size()*M;
  v = new double[m];

  int n = 0;
  //Loop over the rows of blocks
  for ( int bi = 0; bi < x.m(); bi++ )
    //Loop over the rows in a block
    for ( int i = 0; i < x.block_m(bi); i++ )
      //Loop over the rows in a fixed block
      for ( int fi = 0; fi < M; fi++ )
        v[n++] = DLA::index(x[bi][i],fi);

}

template<class VectorS>
void
ScalarVector::setTo( SparseVector< DLA::VectorD<VectorS> >& x)
{
  // Deduce the sizes of VectorS type
  const int M = DLA::VectorSize<VectorS>::M;

  int n = 0;
  //Loop over the rows of blocks
  for ( int bi = 0; bi < x.m(); bi++ )
    //Loop over the rows in a block
    for ( int i = 0; i < x.block_m(bi); i++ )
      //Loop over the rows in a fixed block
      for ( int fi = 0; fi < M; fi++ )
        DLA::index(x[bi][i],fi) = v[n++];
}

//---------------------------------------------------------------------------//
template<class VectorS>
ScalarVector::ScalarVector( const DLA::MatrixDView< SparseVector<VectorS> >& x )
{
  // Deduce the sizes of VectorS type
  const int M = DLA::VectorSize<VectorS>::M;

  SANS_ASSERT( x.n() == 1 ); // assert x is a vector
  m = 0;
  for ( int i = 0; i < x.m(); i++ )
    m += x(i,0).m()*M;

  SANS_ASSERT_MSG(m>0, "m>0 is required to avoid initializing v as a zero-length array");
  v = new double[m];

  int n = 0;
  //Loop over the dynamic rows
  for ( int i = 0; i < x.m(); i++ )
    //Loop over the sparse rows
    for ( int si = 0; si < x(i,0).m(); si++ )
      //Loop over the static block rows
      for ( int bi = 0; bi < M; bi++ )
        v[n++] = DLA::index(x(i,0)[si],bi);
}

template<class VectorS>
void
ScalarVector::setTo( DLA::MatrixDView< SparseVector<VectorS> >& x)
{
  // Deduce the sizes of VectorS type
  const int M = DLA::VectorSize<VectorS>::M;

  SANS_ASSERT( x.n() == 1 );
  int n = 0;
  //Loop over the dynamic rows
  for ( int i = 0; i < x.m(); i++ )
    //Loop over the sparse rows
    for ( int si = 0; si < x(i,0).m(); si++ )
      //Loop over the static block rows
      for ( int bi = 0; bi < M; bi++ )
        DLA::index(x(i,0)[si],bi) = v[n++];
}

//---------------------------------------------------------------------------//
template<class VectorS>
ScalarVector::ScalarVector( const DLA::MatrixDView< SparseVector< DLA::VectorD<VectorS> > >& x )
{
  // Deduce the sizes of VectorS type
  const int M = DLA::VectorSize<VectorS>::M;

  SANS_ASSERT( x.n() == 1 );
  m = 0;
  for ( int i = 0; i < x.m(); i++ )
    for ( int si = 0; si < x(i,0).m(); si++ )
      m += x(i,0)[si].m()*M;

  v = new double[m];

  int n = 0;
  //Loop over the dynamic rows
  for ( int i = 0; i < x.m(); i++ )
    //Loop over the sparse rows
    for ( int si = 0; si < x(i,0).m(); si++ )
      //Loop over the dynamic block rows
      for ( int di = 0; di < x(i,0)[si].m(); di++ )
        //Loop over the static block rows
        for ( int bi = 0; bi < M; bi++ )
          v[n++] = DLA::index(x(i,0)[si][di],bi);
}

template<class VectorS>
void
ScalarVector::setTo( DLA::MatrixDView< SparseVector< DLA::VectorD<VectorS> > >& x)
{
  // Deduce the sizes of VectorS type
  const int M = DLA::VectorSize<VectorS>::M;

  SANS_ASSERT( x.n() == 1 );
  int n = 0;
  //Loop over the dynamic rows
  for ( int i = 0; i < x.m(); i++ )
    //Loop over the sparse rows
    for ( int si = 0; si < x(i,0).m(); si++ )
      //Loop over the dynamic block rows
      for ( int di = 0; di < x(i,0)[si].m(); di++ )
        //Loop over the static block rows
        for ( int bi = 0; bi < M; bi++ )
          DLA::index(x(i,0)[si][di],bi) = v[n++];
}

//---------------------------------------------------------------------------//
template<class Vector0, class Vector1>
ScalarVector::ScalarVector( const BLA::VectorBlock_2< DLA::VectorD<SparseVector<Vector0> >,
                                                      DLA::VectorD<SparseVector<Vector1> > >& x )
{
  // Deduce the sizes of VectorS types
  const int M0 = DLA::VectorSize<Vector0>::M;
  const int M1 = DLA::VectorSize<Vector1>::M;

  m = 0;
  for ( int i = 0; i < x.v0.m(); i++ )
    m += x.v0[i].m()*M0;

  for ( int i = 0; i < x.v1.m(); i++ )
    m += x.v1[i].m()*M1;

  v = new double[m];


  int n = 0;
  //Loop over the dynamic rows
  for ( int i = 0; i < x.v0.m(); i++ )
    //Loop over the sparse rows
    for ( int si = 0; si < x.v0[i].m(); si++ )
      //Loop over the static block rows
      for ( int bi = 0; bi < M0; bi++ )
        v[n++] = DLA::index(x.v0[i][si],bi);

  //Loop over the dynamic rows
  for ( int i = 0; i < x.v1.m(); i++ )
    //Loop over the sparse rows
    for ( int si = 0; si < x.v1[i].m(); si++ )
      //Loop over the static block rows
      for ( int bi = 0; bi < M1; bi++ )
        v[n++] = DLA::index(x.v1[i][si],bi);
}

template<class Vector0, class Vector1>
void
ScalarVector::setTo( BLA::VectorBlock_2< DLA::VectorD<SparseVector<Vector0> >,
                                         DLA::VectorD<SparseVector<Vector1> > >& x )
{
  // Deduce the sizes of VectorS types
  const int M0 = DLA::VectorSize<Vector0>::M;
  const int M1 = DLA::VectorSize<Vector1>::M;

  int n = 0;
  //Loop over the dynamic rows
  for ( int i = 0; i < x.v0.m(); i++ )
    //Loop over the sparse rows
    for ( int si = 0; si < x.v0[i].m(); si++ )
      //Loop over the static block rows
      for ( int bi = 0; bi < M0; bi++ )
        DLA::index(x.v0[i][si],bi) = v[n++];

  //Loop over the dynamic rows
  for ( int i = 0; i < x.v1.m(); i++ )
    //Loop over the sparse rows
    for ( int si = 0; si < x.v1[i].m(); si++ )
      //Loop over the static block rows
      for ( int bi = 0; bi < M1; bi++ )
        DLA::index(x.v1[i][si],bi) = v[n++];
}

//---------------------------------------------------------------------------//
template<class Vector0, class Vector1, class Vector2>
ScalarVector::ScalarVector( const BLA::VectorBlock_3< DLA::VectorD<SparseVector<Vector0> >,
                                                      DLA::VectorD<SparseVector<Vector1> >,
                                                      DLA::VectorD<SparseVector<Vector2> > >& x )
{
  // Deduce the sizes of VectorS types
  const int M0 = DLA::VectorSize<Vector0>::M;
  const int M1 = DLA::VectorSize<Vector1>::M;
  const int M2 = DLA::VectorSize<Vector2>::M;

  m = 0;
  for ( int i = 0; i < x.v0.m(); i++ )
    m += x.v0[i].m()*M0;

  for ( int i = 0; i < x.v1.m(); i++ )
    m += x.v1[i].m()*M1;

  for ( int i = 0; i < x.v2.m(); i++ )
    m += x.v2[i].m()*M2;

  v = new double[m];


  int n = 0;
  //Loop over the dynamic rows
  for ( int i = 0; i < x.v0.m(); i++ )
    //Loop over the sparse rows
    for ( int si = 0; si < x.v0[i].m(); si++ )
      //Loop over the static block rows
      for ( int bi = 0; bi < M0; bi++ )
        v[n++] = DLA::index(x.v0[i][si],bi);

  //Loop over the dynamic rows
  for ( int i = 0; i < x.v1.m(); i++ )
    //Loop over the sparse rows
    for ( int si = 0; si < x.v1[i].m(); si++ )
      //Loop over the static block rows
      for ( int bi = 0; bi < M1; bi++ )
        v[n++] = DLA::index(x.v1[i][si],bi);

  //Loop over the dynamic rows
  for ( int i = 0; i < x.v2.m(); i++ )
    //Loop over the sparse rows
    for ( int si = 0; si < x.v2[i].m(); si++ )
      //Loop over the static block rows
      for ( int bi = 0; bi < M2; bi++ )
        v[n++] = DLA::index(x.v2[i][si],bi);

}

template<class Vector0, class Vector1, class Vector2>
void
ScalarVector::setTo( BLA::VectorBlock_3< DLA::VectorD<SparseVector<Vector0> >,
                                         DLA::VectorD<SparseVector<Vector1> >,
                                         DLA::VectorD<SparseVector<Vector2> > >& x )
{
  // Deduce the sizes of VectorS types
  const int M0 = DLA::VectorSize<Vector0>::M;
  const int M1 = DLA::VectorSize<Vector1>::M;
  const int M2 = DLA::VectorSize<Vector2>::M;

  int n = 0;
  //Loop over the dynamic rows
  for ( int i = 0; i < x.v0.m(); i++ )
    //Loop over the sparse rows
    for ( int si = 0; si < x.v0[i].m(); si++ )
      //Loop over the static block rows
      for ( int bi = 0; bi < M0; bi++ )
        DLA::index(x.v0[i][si],bi) = v[n++];

  //Loop over the dynamic rows
  for ( int i = 0; i < x.v1.m(); i++ )
    //Loop over the sparse rows
    for ( int si = 0; si < x.v1[i].m(); si++ )
      //Loop over the static block rows
      for ( int bi = 0; bi < M1; bi++ )
        DLA::index(x.v1[i][si],bi) = v[n++];

  //Loop over the dynamic rows
  for ( int i = 0; i < x.v2.m(); i++ )
    //Loop over the sparse rows
    for ( int si = 0; si < x.v2[i].m(); si++ )
      //Loop over the static block rows
      for ( int bi = 0; bi < M2; bi++ )
        DLA::index(x.v2[i][si],bi) = v[n++];
}

//---------------------------------------------------------------------------//
template<class Vector0, class Vector1, class Vector2, class Vector3>
ScalarVector::
ScalarVector( const BLA::VectorBlock_4< DLA::VectorD<SparseVector<Vector0> >,
                                        DLA::VectorD<SparseVector<Vector1> >,
                                        DLA::VectorD<SparseVector<Vector2> >,
                                        DLA::VectorD<SparseVector<Vector3> > >& x )
{
  // Deduce the sizes of VectorS types
  const int M0 = DLA::VectorSize<Vector0>::M;
  const int M1 = DLA::VectorSize<Vector1>::M;
  const int M2 = DLA::VectorSize<Vector2>::M;
  const int M3 = DLA::VectorSize<Vector3>::M;

  m = 0;
  for ( int i = 0; i < x.v0.m(); i++ )
    m += x.v0[i].m()*M0;

  for ( int i = 0; i < x.v1.m(); i++ )
    m += x.v1[i].m()*M1;

  for ( int i = 0; i < x.v2.m(); i++ )
    m += x.v2[i].m()*M2;

  for ( int i = 0; i < x.v3.m(); i++ )
    m += x.v3[i].m()*M3;

  v = new double[m];

  int n = 0;
  //Loop over the dynamic rows
  for ( int i = 0; i < x.v0.m(); i++ )
    //Loop over the sparse rows
    for ( int si = 0; si < x.v0[i].m(); si++ )
      //Loop over the static block rows
      for ( int bi = 0; bi < M0; bi++ )
        v[n++] = DLA::index(x.v0[i][si],bi);

  //Loop over the dynamic rows
  for ( int i = 0; i < x.v1.m(); i++ )
    //Loop over the sparse rows
    for ( int si = 0; si < x.v1[i].m(); si++ )
      //Loop over the static block rows
      for ( int bi = 0; bi < M1; bi++ )
        v[n++] = DLA::index(x.v1[i][si],bi);

  //Loop over the dynamic rows
  for ( int i = 0; i < x.v2.m(); i++ )
    //Loop over the sparse rows
    for ( int si = 0; si < x.v2[i].m(); si++ )
      //Loop over the static block rows
      for ( int bi = 0; bi < M2; bi++ )
        v[n++] = DLA::index(x.v2[i][si],bi);

  //Loop over the dynamic rows
  for ( int i = 0; i < x.v3.m(); i++ )
    //Loop over the sparse rows
    for ( int si = 0; si < x.v3[i].m(); si++ )
      //Loop over the static block rows
      for ( int bi = 0; bi < M3; bi++ )
        v[n++] = DLA::index(x.v3[i][si],bi);

}

template<class Vector0, class Vector1, class Vector2, class Vector3>
void
ScalarVector::
setTo( BLA::VectorBlock_4< DLA::VectorD<SparseVector<Vector0> >,
                           DLA::VectorD<SparseVector<Vector1> >,
                           DLA::VectorD<SparseVector<Vector2> >,
                           DLA::VectorD<SparseVector<Vector3> > >& x )
{
  // Deduce the sizes of VectorS types
  const int M0 = DLA::VectorSize<Vector0>::M;
  const int M1 = DLA::VectorSize<Vector1>::M;
  const int M2 = DLA::VectorSize<Vector2>::M;
  const int M3 = DLA::VectorSize<Vector3>::M;

  int n = 0;
  //Loop over the dynamic rows
  for ( int i = 0; i < x.v0.m(); i++ )
    //Loop over the sparse rows
    for ( int si = 0; si < x.v0[i].m(); si++ )
      //Loop over the static block rows
      for ( int bi = 0; bi < M0; bi++ )
        DLA::index(x.v0[i][si],bi) = v[n++];

  //Loop over the dynamic rows
  for ( int i = 0; i < x.v1.m(); i++ )
    //Loop over the sparse rows
    for ( int si = 0; si < x.v1[i].m(); si++ )
      //Loop over the static block rows
      for ( int bi = 0; bi < M1; bi++ )
        DLA::index(x.v1[i][si],bi) = v[n++];

  //Loop over the dynamic rows
  for ( int i = 0; i < x.v2.m(); i++ )
    //Loop over the sparse rows
    for ( int si = 0; si < x.v2[i].m(); si++ )
      //Loop over the static block rows
      for ( int bi = 0; bi < M2; bi++ )
        DLA::index(x.v2[i][si],bi) = v[n++];

  //Loop over the dynamic rows
  for ( int i = 0; i < x.v3.m(); i++ )
    //Loop over the sparse rows
    for ( int si = 0; si < x.v3[i].m(); si++ )
      //Loop over the static block rows
      for ( int bi = 0; bi < M3; bi++ )
        DLA::index(x.v3[i][si],bi) = v[n++];
}

} // namespace SLA
} // namespace SANS
