// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXTYPE_H
#define MATRIXTYPE_H

#include "tools/SANSnumerics.h"

#include "dense/dynamic/MatrixD_Type.h"
#include "dense/static/MatrixS_Type.h"
#include "block/block_Type.h"
#include "sparse/sparse_Type.h"

namespace tinymat 
{

//=============================================================================
//A template metafunction for choosing the appropriate matrix type for a given vector type
template< class Vector_type >
struct MatrixType;

template< >
struct MatrixType< SLA::SparseVector<Real> >
{
  typedef SLA::SparseMatrix_CRS<Real> type;
};

template< >
struct MatrixType< SLA::SparseVector< DLA::VectorD<Real> > >
{
  typedef SLA::SparseMatrix_CRS< DLA::MatrixD<Real> > type;
};

template< int M >
struct MatrixType< SLA::SparseVector< DLA::VectorS<M,Real> > >
{
  typedef SLA::SparseMatrix_CRS< DLA::MatrixS<M,M,Real> > type;
};

template< int M >
struct MatrixType< SLA::SparseVector< DLA::VectorD< DLA::VectorS<M,Real> > > >
{
  typedef SLA::SparseMatrix_CRS< DLA::MatrixD<DLA::MatrixS<M,M,Real> > > type;
};

template<>
struct MatrixType< DLA::VectorD< SLA::SparseVector<Real> > >
{
  typedef DLA::MatrixD< SLA::SparseMatrix_CRS<Real> > type;
};

template< >
struct MatrixType< DLA::VectorD< SLA::SparseVector< DLA::VectorD<Real> > > >
{
  typedef DLA::MatrixD< SLA::SparseMatrix_CRS< DLA::MatrixD<Real> > > type;
};

template< int M >
struct MatrixType< DLA::VectorD< SLA::SparseVector< DLA::VectorS<M,Real> > > >
{
  typedef DLA::MatrixD< SLA::SparseMatrix_CRS< DLA::MatrixS<M,M,Real> > > type;
};

template< int M >
struct MatrixType< DLA::VectorD< SLA::SparseVector< DLA::VectorD< DLA::VectorS<M,Real> > > > >
{
  typedef DLA::MatrixD< SLA::SparseMatrix_CRS< DLA::MatrixD<DLA::MatrixS<M,M,Real> > > > type;
};

} //namespace tinymat 


#endif //MATRIXTYPE_H