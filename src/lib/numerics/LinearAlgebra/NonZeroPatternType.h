// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef NONZEROPATTERNTYPE_H
#define NONZEROPATTERNTYPE_H

#include "tools/SANSnumerics.h"

#include "DenseLinAlg/DynamicSize/MatrixD_Type.h"
#include "DenseLinAlg/StaticSize/MatrixS_Type.h"
#include "BlockLinAlg/BlockLinAlg_Type.h"
#include "SparseLinAlg/SparseLinAlg_Type.h"

namespace SANS
{

//=============================================================================
//A template metafunction for choosing the appropriate non-zero pattern for a given matrix type
template< class Matrix_type >
struct NonZeroPatternType;

template< >
struct NonZeroPatternType< DLA::MatrixD<Real> >
{
  typedef DLA::DenseNonZeroPattern<Real> type;
  typedef DLA::DenseNonZeroPattern<Real> Viewtype;
};

template< >
struct NonZeroPatternType< DLA::MatrixD<DLA::MatrixD<Real>> >
{
  typedef DLA::MatrixD    <DLA::DenseNonZeroPattern<Real>> type;
  typedef DLA::MatrixDView<DLA::DenseNonZeroPattern<Real>> Viewtype;
};

template< int M, int N >
struct NonZeroPatternType< DLA::MatrixD<DLA::MatrixD<DLA::MatrixS<M,N,Real>>> >
{
  typedef DLA::MatrixD    <DLA::DenseNonZeroPattern<DLA::MatrixS<M,N,Real>>> type;
  typedef DLA::MatrixDView<DLA::DenseNonZeroPattern<DLA::MatrixS<M,N,Real>>> Viewtype;
};

template< >
struct NonZeroPatternType< SLA::SparseMatrix_CRS<Real> >
{
  typedef SLA::SparseNonZeroPattern<Real> type;
  typedef SLA::SparseNonZeroPattern<Real> Viewtype;
};

template< >
struct NonZeroPatternType< SLA::SparseMatrix_CRS< DLA::MatrixD<Real> > >
{
  typedef SLA::SparseNonZeroPattern< DLA::MatrixD<Real> > type;
  typedef SLA::SparseNonZeroPattern< DLA::MatrixD<Real> > Viewtype;
};

template< int M, int N >
struct NonZeroPatternType< SLA::SparseMatrix_CRS< DLA::MatrixS<M, N, Real> > >
{
  typedef SLA::SparseNonZeroPattern< DLA::MatrixS<M,N,Real> > type;
  typedef SLA::SparseNonZeroPattern< DLA::MatrixS<M,N,Real> > Viewtype;
};

template< int M, int N >
struct NonZeroPatternType< SLA::SparseMatrix_CRS< DLA::MatrixD<DLA::MatrixS<M,N,Real> > > >
{
  typedef SLA::SparseNonZeroPattern< DLA::MatrixD< DLA::MatrixS<M,N,Real> > > type;
  typedef SLA::SparseNonZeroPattern< DLA::MatrixD< DLA::MatrixS<M,N,Real> > > Viewtype;
};

template<>
struct NonZeroPatternType< DLA::MatrixD< SLA::SparseMatrix_CRS<Real> > >
{
  typedef DLA::MatrixD    < SLA::SparseNonZeroPattern<Real> > type;
  typedef DLA::MatrixDView< SLA::SparseNonZeroPattern<Real> > Viewtype;
};

template< >
struct NonZeroPatternType< DLA::MatrixD< SLA::SparseMatrix_CRS< DLA::MatrixD<Real> > > >
{
  typedef DLA::MatrixD    < SLA::SparseNonZeroPattern< DLA::MatrixD<Real> > > type;
  typedef DLA::MatrixDView< SLA::SparseNonZeroPattern< DLA::MatrixD<Real> > > Viewtype;
};

template< int M, int N >
struct NonZeroPatternType< DLA::MatrixD< SLA::SparseMatrix_CRS< DLA::MatrixS<M,N,Real> > > >
{
  typedef DLA::MatrixD    < SLA::SparseNonZeroPattern< DLA::MatrixS<M,N,Real> > > type;
  typedef DLA::MatrixDView< SLA::SparseNonZeroPattern< DLA::MatrixS<M,N,Real> > > Viewtype;
};

template< int M, int N >
struct NonZeroPatternType< DLA::MatrixD< SLA::SparseMatrix_CRS< DLA::MatrixD<DLA::MatrixS<M,N,Real> > > > >
{
  typedef DLA::MatrixD    < SLA::SparseNonZeroPattern< DLA::MatrixD< DLA::MatrixS<M,N,Real> > > > type;
  typedef DLA::MatrixDView< SLA::SparseNonZeroPattern< DLA::MatrixD< DLA::MatrixS<M,N,Real> > > > Viewtype;
};

template<  >
struct NonZeroPatternType<
         BLA::MatrixBlock_2x2<
           DLA::MatrixD<DLA::MatrixD<Real>>, DLA::MatrixD<DLA::MatrixD<Real>>,
           DLA::MatrixD<DLA::MatrixD<Real>>, DLA::MatrixD<DLA::MatrixD<Real>> > >
{
  typedef BLA::MatrixBlock_2x2<
            DLA::MatrixD<DLA::DenseNonZeroPattern<Real>>, DLA::MatrixD<DLA::DenseNonZeroPattern<Real>>,
            DLA::MatrixD<DLA::DenseNonZeroPattern<Real>>, DLA::MatrixD<DLA::DenseNonZeroPattern<Real>> > type;
  typedef type Viewtype;
};

template< int M, int N >
struct NonZeroPatternType<
         BLA::MatrixBlock_2x2<
           DLA::MatrixD<DLA::MatrixD<DLA::MatrixS<M,N,Real>>>, DLA::MatrixD<DLA::MatrixD<DLA::MatrixS<M,1,Real>>>,
           DLA::MatrixD<DLA::MatrixD<DLA::MatrixS<1,N,Real>>>, DLA::MatrixD<DLA::MatrixD<Real>>                  > >
{
  typedef BLA::MatrixBlock_2x2<
            DLA::MatrixD<DLA::DenseNonZeroPattern<DLA::MatrixS<M,N,Real>>>, DLA::MatrixD<DLA::DenseNonZeroPattern<DLA::MatrixS<M,1,Real>>>,
            DLA::MatrixD<DLA::DenseNonZeroPattern<DLA::MatrixS<1,N,Real>>>, DLA::MatrixD<DLA::DenseNonZeroPattern<Real>>                  > type;
  typedef type Viewtype;
};

template<  >
struct NonZeroPatternType<
         BLA::MatrixBlock_2x2<
           DLA::MatrixD<SLA::SparseMatrix_CRS<Real>>, DLA::MatrixD<SLA::SparseMatrix_CRS<Real>>,
           DLA::MatrixD<SLA::SparseMatrix_CRS<Real>>, DLA::MatrixD<SLA::SparseMatrix_CRS<Real>> > >
{
  typedef BLA::MatrixBlock_2x2<
            DLA::MatrixD<SLA::SparseNonZeroPattern<Real>>, DLA::MatrixD<SLA::SparseNonZeroPattern<Real>>,
            DLA::MatrixD<SLA::SparseNonZeroPattern<Real>>, DLA::MatrixD<SLA::SparseNonZeroPattern<Real>> > type;
  typedef type Viewtype;
};

template< int M, int N >
struct NonZeroPatternType<
         BLA::MatrixBlock_2x2<
           DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<M,N,Real>>>, DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<M,1,Real>>>,
           DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<1,N,Real>>>, DLA::MatrixD<SLA::SparseMatrix_CRS<Real>>                  > >
{
  typedef BLA::MatrixBlock_2x2<
            DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<M,N,Real>>>, DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<M,1,Real>>>,
            DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<1,N,Real>>>, DLA::MatrixD<SLA::SparseNonZeroPattern<Real>>                  > type;
  typedef type Viewtype;
};

template< int M0, int M1 >
struct NonZeroPatternType<
         BLA::MatrixBlock_2x2<
           DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<M0,M0,Real>>>, DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<M0,M1,Real>>>,
           DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<M1,M0,Real>>>, DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<M1,M1,Real>>> > >
{
  typedef BLA::MatrixBlock_2x2<
            DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<M0,M0,Real>>>, DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<M0,M1,Real>>>,
            DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<M1,M0,Real>>>, DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<M1,M1,Real>>>
                              > type;
  typedef type Viewtype;
};

template<int L, int M, int N, int J, int K>
struct NonZeroPatternType<
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
  typedef BLA::MatrixBlock_3x3<
            DLA::MatrixD<SLA::SparseNonZeroPattern<Real                  >>,
            DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<1,M,Real>>>,
            DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<1,N,Real>>>,

            DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<J,L,Real>>>,
            DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<J,M,Real>>>,
            DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<J,N,Real>>>,

            DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<K,L,Real>>>,
            DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<K,M,Real>>>,
            DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<K,N,Real>>> > type;
  typedef type Viewtype;
};

template<int M0, int M1, int M2, int M3>
struct NonZeroPatternType<
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
                              DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<M3,M3,Real> > >
                             >
                 >
{
  typedef BLA::MatrixBlock_4x4< DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<M0,M0,Real> > >,
                                DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<M0,M1,Real> > >,
                                DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<M0,M2,Real> > >,
                                DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<M0,M3,Real> > >,
                                //
                                DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<M1,M0,Real> > >,
                                DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<M1,M1,Real> > >,
                                DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<M1,M2,Real> > >,
                                DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<M1,M3,Real> > >,
                                //
                                DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<M2,M0,Real> > >,
                                DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<M2,M1,Real> > >,
                                DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<M2,M2,Real> > >,
                                DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<M2,M3,Real> > >,
                                //
                                DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<M3,M0,Real> > >,
                                DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<M3,M1,Real> > >,
                                DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<M3,M2,Real> > >,
                                DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<M3,M3,Real> > >
                              > type;
  typedef type Viewtype;
};

// partial specialization for IBL2D/panel coupling
template<int M0, int M2>
struct NonZeroPatternType<
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
                              DLA::MatrixD<SLA::SparseMatrix_CRS<Real> > > >
{
  typedef BLA::MatrixBlock_4x4< DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<M0,M0,Real> > >,
                                DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<M0,1,Real> > >,
                                DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<M0,M2,Real> > >,
                                DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<M0,1,Real> > >,
                                //
                                DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<1,M0,Real> > >,
                                DLA::MatrixD<SLA::SparseNonZeroPattern<Real> >,
                                DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<1,M2,Real> > >,
                                DLA::MatrixD<SLA::SparseNonZeroPattern<Real> >,
                                //
                                DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<M2,M0,Real> > >,
                                DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<M2,1,Real> > >,
                                DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<M2,M2,Real> > >,
                                DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<M2,1,Real> > >,
                                //
                                DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<1,M0,Real> > >,
                                DLA::MatrixD<SLA::SparseNonZeroPattern<Real> >,
                                DLA::MatrixD<SLA::SparseNonZeroPattern<DLA::MatrixS<1,M2,Real> > >,
                                DLA::MatrixD<SLA::SparseNonZeroPattern<Real> >
                              > type;
  typedef type Viewtype;
};

} //namespace SANS


#endif //NONZEROPATTERNTYPE_H
