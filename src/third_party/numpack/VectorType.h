// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef VECTORTYPE_H
#define VECTORTYPE_H

#include "tools/SANSnumerics.h"

#include "DenseLinAlg/DynamicSize/MatrixD_Type.h"
#include "DenseLinAlg/StaticSize/MatrixS_Type.h"
#include "BlockLinAlg/BlockLinAlg_Type.h"
#include "SparseLinAlg/SparseLinAlg_Type.h"

namespace numpack 
{

//=============================================================================
//A template metafunction for choosing the appropriate vector for a given matrix type
template< class Matrix_type >
struct VectorType;

// specializations for different matrix types as follows
//
template< >
struct VectorType< DLA::MatrixD<Real> >
{
  typedef DLA::VectorD    <Real> type;
  typedef DLA::VectorDView<Real> Viewtype;
};

template< >
struct VectorType< DLA::MatrixD<DLA::MatrixD<Real>> >
{
  typedef DLA::VectorD    <DLA::VectorD<Real>> type;
  typedef DLA::VectorDView<DLA::VectorD<Real>> Viewtype;
};

template< int M, int N >
struct VectorType< DLA::MatrixD<DLA::MatrixD<DLA::MatrixS<M,N,Real>>> >
{
  typedef DLA::VectorD    <DLA::VectorD<DLA::VectorS<N,Real>>> type;
  typedef DLA::VectorDView<DLA::VectorD<DLA::VectorS<N,Real>>> Viewtype;
};

template< >
struct VectorType< SLA::SparseMatrix_CRS<Real> >
{
  typedef SLA::SparseVector<Real> type;
  typedef SLA::SparseVector<Real> Viewtype;
};

template< >
struct VectorType< SLA::SparseMatrix_CRS< DLA::MatrixD<Real> > >
{
  typedef SLA::SparseVector< DLA::VectorD<Real> > type;
  typedef SLA::SparseVector< DLA::VectorD<Real> > Viewtype;
};

template< int M >
struct VectorType< SLA::SparseMatrix_CRS< DLA::MatrixS<M, 1, Real> > >
{
  typedef SLA::SparseVector< Real > type;
  typedef SLA::SparseVector< Real > Viewtype;
};

template< int M, int N >
struct VectorType< SLA::SparseMatrix_CRS< DLA::MatrixS<M, N, Real> > >
{
  typedef SLA::SparseVector< DLA::VectorS<N,Real> > type;
  typedef SLA::SparseVector< DLA::VectorS<N,Real> > Viewtype;
};

template< int M, int N >
struct VectorType< SLA::SparseMatrix_CRS< DLA::MatrixD<DLA::MatrixS<M,N,Real> > > >
{
  typedef SLA::SparseVector< DLA::VectorD< DLA::VectorS<N,Real> > > type;
  typedef SLA::SparseVector< DLA::VectorD< DLA::VectorS<N,Real> > > Viewtype;
};

template<>
struct VectorType< DLA::MatrixD< SLA::SparseMatrix_CRS<Real> > >
{
  typedef DLA::VectorD    < SLA::SparseVector<Real> > type;
  typedef DLA::VectorDView< SLA::SparseVector<Real> > Viewtype;
};

template< >
struct VectorType< DLA::MatrixD< SLA::SparseMatrix_CRS< DLA::MatrixD<Real> > > >
{
  typedef DLA::VectorD    < SLA::SparseVector< DLA::VectorD<Real> > > type;
  typedef DLA::VectorDView< SLA::SparseVector< DLA::VectorD<Real> > > Viewtype;
};

template< int M >
struct VectorType< DLA::MatrixD< SLA::SparseMatrix_CRS< DLA::MatrixS<M,1,Real> > > >
{
  typedef DLA::VectorD    < SLA::SparseVector< Real > > type;
  typedef DLA::VectorDView< SLA::SparseVector< Real > > Viewtype;
};

template< int M, int N >
struct VectorType< DLA::MatrixD< SLA::SparseMatrix_CRS< DLA::MatrixS<M,N,Real> > > >
{
  typedef DLA::VectorD    < SLA::SparseVector< DLA::VectorS<N,Real> > > type;
  typedef DLA::VectorDView< SLA::SparseVector< DLA::VectorS<N,Real> > > Viewtype;
};

template< int M, int N >
struct VectorType< DLA::MatrixD< SLA::SparseMatrix_CRS< DLA::MatrixD<DLA::MatrixS<M,N,Real> > > > >
{
  typedef DLA::VectorD    < SLA::SparseVector< DLA::VectorD< DLA::VectorS<N,Real> > > > type;
  typedef DLA::VectorDView< SLA::SparseVector< DLA::VectorD< DLA::VectorS<N,Real> > > > Viewtype;
};

template<  >
struct VectorType<
         BLA::MatrixBlock_2x2<
           DLA::MatrixD<DLA::MatrixD<Real>>, DLA::MatrixD<DLA::MatrixD<Real>>,
           DLA::MatrixD<DLA::MatrixD<Real>>, DLA::MatrixD<DLA::MatrixD<Real>> > >
{
  typedef BLA::VectorBlock_2< DLA::VectorD<DLA::VectorD<Real>>,
                              DLA::VectorD<DLA::VectorD<Real>> > type;
  typedef type Viewtype;
};

template< int M, int N >
struct VectorType<
         BLA::MatrixBlock_2x2<
           DLA::MatrixD<DLA::MatrixD<DLA::MatrixS<M,N,Real>>>, DLA::MatrixD<DLA::MatrixD<DLA::MatrixS<M,1,Real>>>,
           DLA::MatrixD<DLA::MatrixD<DLA::MatrixS<1,N,Real>>>, DLA::MatrixD<DLA::MatrixD<Real>>                  > >
{
  typedef BLA::VectorBlock_2< DLA::VectorD<DLA::VectorD<DLA::VectorS<N,Real>>>, DLA::VectorD<DLA::VectorD<Real>> > type;
  typedef type Viewtype;
};

template<  >
struct VectorType<
         BLA::MatrixBlock_2x2<
           DLA::MatrixD<SLA::SparseMatrix_CRS<Real>>, DLA::MatrixD<SLA::SparseMatrix_CRS<Real>>,
           DLA::MatrixD<SLA::SparseMatrix_CRS<Real>>, DLA::MatrixD<SLA::SparseMatrix_CRS<Real>> > >
{
  typedef BLA::VectorBlock_2< DLA::VectorD<SLA::SparseVector<Real>>,
                              DLA::VectorD<SLA::SparseVector<Real>> > type;
  typedef type Viewtype;
};

template< int M, int N >
struct VectorType<
         BLA::MatrixBlock_2x2<
           DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<M,N,Real>>>, DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<M,1,Real>>>,
           DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<1,N,Real>>>, DLA::MatrixD<SLA::SparseMatrix_CRS<Real>>                  > >
{
  typedef BLA::VectorBlock_2< DLA::VectorD<SLA::SparseVector<DLA::VectorS<N,Real>>>, DLA::VectorD<SLA::SparseVector<Real>> > type;
  typedef type Viewtype;
};

template< int M0, int M1 >
struct VectorType<
         BLA::MatrixBlock_2x2<
           DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<M0,M0,Real>>>, DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<M0,M1,Real>>>,
           DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<M1,M0,Real>>>, DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<M1,M1,Real>>>
                 >
                         >
{
  typedef BLA::VectorBlock_2< DLA::VectorD<SLA::SparseVector<DLA::VectorS<M0,Real>>>,
                              DLA::VectorD<SLA::SparseVector<DLA::VectorS<M1,Real>>> > type;
  typedef type Viewtype;
};

template<int L, int M, int N, int J, int K>
struct VectorType<
         BLA::MatrixBlock_3x3< DLA::MatrixD<SLA::SparseMatrix_CRS<Real                  >>,
                               DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<1,M,Real>>>,
                               DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<1,N,Real>>>,

                               DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<J,L,Real>>>,
                               DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<J,M,Real>>>,
                               DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<J,N,Real>>>,

                               DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<K,L,Real>>>,
                               DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<K,M,Real>>>,
                               DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<K,N,Real>>> > >
{
  typedef BLA::VectorBlock_3< DLA::VectorD<SLA::SparseVector<Real>>,
                              DLA::VectorD<SLA::SparseVector<DLA::VectorS<M,Real>>>,
                              DLA::VectorD<SLA::SparseVector<DLA::VectorS<N,Real>>> > type;
  typedef type Viewtype;
};

template<int M0, int M1, int M2, int M3>
struct VectorType<
         BLA::MatrixBlock_4x4<DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<M0,M0,Real> > >,
                              DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<M0,M1,Real> > >,
                              DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<M0,M2,Real> > >,
                              DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<M0,M3,Real> > >,
                              //
                              DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<M1,M0,Real> > >,
                              DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<M1,M1,Real> > >,
                              DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<M1,M2,Real> > >,
                              DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<M1,M3,Real> > >,
                              //
                              DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<M2,M0,Real> > >,
                              DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<M2,M1,Real> > >,
                              DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<M2,M2,Real> > >,
                              DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<M2,M3,Real> > >,
                              //
                              DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<M3,M0,Real> > >,
                              DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<M3,M1,Real> > >,
                              DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<M3,M2,Real> > >,
                              DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<M3,M3,Real> > > >
                 >
{
  typedef BLA::VectorBlock_4< DLA::VectorD<SLA::SparseVector<DLA::VectorS<M0,Real> > >,
                              DLA::VectorD<SLA::SparseVector<DLA::VectorS<M1,Real> > >,
                              DLA::VectorD<SLA::SparseVector<DLA::VectorS<M2,Real> > >,
                              DLA::VectorD<SLA::SparseVector<DLA::VectorS<M3,Real> > > > type;
  typedef type Viewtype;
};


// partial specialization for IBL2D/panel coupling
template<int M0, int M2>
struct VectorType<
         BLA::MatrixBlock_4x4<DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<M0,M0,Real> > >,
                              DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<M0,1,Real> > >,
                              DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<M0,M2,Real> > >,
                              DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<M0,1,Real> > >,
                              //
                              DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<1,M0,Real> > >,
                              DLA::MatrixD<SLA::SparseMatrix_CRS<Real> >,
                              DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<1,M2,Real> > >,
                              DLA::MatrixD<SLA::SparseMatrix_CRS<Real> >,
                              //
                              DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<M2,M0,Real> > >,
                              DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<M2,1,Real> > >,
                              DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<M2,M2,Real> > >,
                              DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<M2,1,Real> > >,
                              //
                              DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<1,M0,Real> > >,
                              DLA::MatrixD<SLA::SparseMatrix_CRS<Real> >,
                              DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<1,M2,Real> > >,
                              DLA::MatrixD<SLA::SparseMatrix_CRS<Real> > >
                 >
{
  typedef BLA::VectorBlock_4< DLA::VectorD<SLA::SparseVector<DLA::VectorS<M0,Real> > >,
                              DLA::VectorD<SLA::SparseVector<Real> >,
                              DLA::VectorD<SLA::SparseVector<DLA::VectorS<M2,Real> > >,
                              DLA::VectorD<SLA::SparseVector<Real> > > type;
  typedef type Viewtype;
};


} //namespace numpack 


#endif //VECTORTYPE_H
