// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef VECTORSIZETYPE_H
#define VECTORSIZETYPE_H

#include "tools/SANSnumerics.h"

#include "dense/dynamic/MatrixD_Type.h"
#include "dense/static/MatrixS_Type.h"
#include "block/block_Type.h"
#include "sparse/sparse_Type.h"

namespace tinymat 
{

//=============================================================================
//A template metafunction for getting the appropriate vector size type
template< class Vector_type >
struct VectorSizeType;

template< class TV >
struct VectorSizeType< SLA::SparseVector<TV> >
{
  typedef typename SLA::SparseVector<TV>::size_type type;
};

template< class TV >
struct VectorSizeType< DLA::VectorD< SLA::SparseVector<TV> > >
{
  typedef DLA::VectorD< typename SLA::SparseVector<TV>::size_type > type;
};

template< class TV >
struct VectorSizeType< DLA::VectorD<TV> >
{
  typedef DLA::DenseVectorSize type;
};

template< class TV >
struct VectorSizeType< DLA::VectorD< DLA::VectorD<TV> > >
{
  typedef DLA::VectorD< DLA::DenseVectorSize > type;
};

//---------------------------------------------------------------------------//
template< class TV0, class TV1 >
struct VectorSizeType<BLA::VectorBlock_2<DLA::VectorD<SLA::SparseVector<TV0> >,
                                         DLA::VectorD<SLA::SparseVector<TV1> > > >
{
  typedef BLA::VectorBlock_2< DLA::VectorD< typename SLA::SparseVector<TV0>::size_type >,
                              DLA::VectorD< typename SLA::SparseVector<TV1>::size_type > > type;
};

template< class TV0, class TV1 >
struct VectorSizeType<BLA::VectorBlock_2<DLA::VectorD<DLA::VectorD<TV0> >,
                                         DLA::VectorD<DLA::VectorD<TV1> > > >
{
  typedef BLA::VectorBlock_2< DLA::VectorD< DLA::DenseVectorSize >,
                              DLA::VectorD< DLA::DenseVectorSize > > type;
};

//---------------------------------------------------------------------------//
template< class TV0, class TV1, class TV2 >
struct VectorSizeType<BLA::VectorBlock_3<DLA::VectorD<SLA::SparseVector<TV0> >,
                                         DLA::VectorD<SLA::SparseVector<TV1> >,
                                         DLA::VectorD<SLA::SparseVector<TV2> > > >
{
  typedef BLA::VectorBlock_3< DLA::VectorD< typename SLA::SparseVector<TV0>::size_type >,
                              DLA::VectorD< typename SLA::SparseVector<TV1>::size_type >,
                              DLA::VectorD< typename SLA::SparseVector<TV2>::size_type > > type;
};

template< class TV0, class TV1, class TV2 >
struct VectorSizeType<BLA::VectorBlock_3<DLA::VectorD<DLA::VectorD<TV0> >,
                                         DLA::VectorD<DLA::VectorD<TV1> >,
                                         DLA::VectorD<DLA::VectorD<TV2> > > >
{
  typedef BLA::VectorBlock_3< DLA::VectorD< DLA::DenseVectorSize >,
                              DLA::VectorD< DLA::DenseVectorSize >,
                              DLA::VectorD< DLA::DenseVectorSize > > type;
};

//---------------------------------------------------------------------------//
template< class TV0, class TV1, class TV2, class TV3 >
struct VectorSizeType<BLA::VectorBlock_4<DLA::VectorD<SLA::SparseVector<TV0> >,
                                         DLA::VectorD<SLA::SparseVector<TV1> >,
                                         DLA::VectorD<SLA::SparseVector<TV2> >,
                                         DLA::VectorD<SLA::SparseVector<TV3> > > >
{
  typedef BLA::VectorBlock_4< DLA::VectorD< typename SLA::SparseVector<TV0>::size_type >,
                              DLA::VectorD< typename SLA::SparseVector<TV1>::size_type >,
                              DLA::VectorD< typename SLA::SparseVector<TV2>::size_type >,
                              DLA::VectorD< typename SLA::SparseVector<TV3>::size_type > > type;
};

template< class TV0, class TV1, class TV2, class TV3 >
struct VectorSizeType<BLA::VectorBlock_4<DLA::VectorD<DLA::VectorD<TV0> >,
                                         DLA::VectorD<DLA::VectorD<TV1> >,
                                         DLA::VectorD<DLA::VectorD<TV2> >,
                                         DLA::VectorD<DLA::VectorD<TV3> > > >
{
  typedef BLA::VectorBlock_4< DLA::VectorD< DLA::DenseVectorSize >,
                              DLA::VectorD< DLA::DenseVectorSize >,
                              DLA::VectorD< DLA::DenseVectorSize >,
                              DLA::VectorD< DLA::DenseVectorSize > > type;
};

} //namespace tinymat 


#endif //VECTORSIZETYPE_H
