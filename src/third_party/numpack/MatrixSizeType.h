// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXSIZETYPE_H
#define MATRIXSIZETYPE_H

#include "tools/SANSnumerics.h"

#include "dense/dynamic/MatrixD_Type.h"
#include "dense/static/MatrixS_Type.h"
#include "block/block_Type.h"
#include "sparse/sparse_Type.h"

namespace numpack 
{

//=============================================================================
//A template metafunction for choosing the appropriate matrix size type for a given matrix type
template< class Matrix_type >
struct MatrixSizeType;

template< class TM >
struct MatrixSizeType< SLA::SparseMatrix_CRS<TM> >
{
  typedef SLA::SparseMatrixSize type;
};

template< class TM >
struct MatrixSizeType< DLA::MatrixD< SLA::SparseMatrix_CRS<TM> > >
{
  typedef DLA::MatrixD< SLA::SparseMatrixSize > type;
};

template< class TM >
struct MatrixSizeType< DLA::MatrixD<TM> >
{
  typedef DLA::DenseMatrixSize type;
};

template< class TM >
struct MatrixSizeType< DLA::MatrixD< DLA::MatrixD<TM> > >
{
  typedef DLA::MatrixD< DLA::DenseMatrixSize > type;
};

template< class TM00, class TM01, class TM10, class TM11 >
struct MatrixSizeType<BLA::MatrixBlock_2x2<DLA::MatrixD<SLA::SparseMatrix_CRS<TM00> >,
                                           DLA::MatrixD<SLA::SparseMatrix_CRS<TM01> >,
                                           DLA::MatrixD<SLA::SparseMatrix_CRS<TM10> >,
                                           DLA::MatrixD<SLA::SparseMatrix_CRS<TM11> > > >
{
  typedef BLA::MatrixBlock_2x2< DLA::MatrixD< SLA::SparseMatrixSize >,
                                DLA::MatrixD< SLA::SparseMatrixSize >,
                                DLA::MatrixD< SLA::SparseMatrixSize >,
                                DLA::MatrixD< SLA::SparseMatrixSize > > type;
};

template< class TM00, class TM01, class TM10, class TM11 >
struct MatrixSizeType<BLA::MatrixBlock_2x2<DLA::MatrixD<DLA::MatrixD<TM00> >,
                                           DLA::MatrixD<DLA::MatrixD<TM01> >,
                                           DLA::MatrixD<DLA::MatrixD<TM10> >,
                                           DLA::MatrixD<DLA::MatrixD<TM11> > > >
{
  typedef BLA::MatrixBlock_2x2< DLA::MatrixD< DLA::DenseMatrixSize >,
                                DLA::MatrixD< DLA::DenseMatrixSize >,
                                DLA::MatrixD< DLA::DenseMatrixSize >,
                                DLA::MatrixD< DLA::DenseMatrixSize > > type;
};

template< class TM00, class TM01, class TM02,
          class TM10, class TM11, class TM12,
          class TM20, class TM21, class TM22>
struct MatrixSizeType<BLA::MatrixBlock_3x3<DLA::MatrixD<SLA::SparseMatrix_CRS<TM00> >,
                                           DLA::MatrixD<SLA::SparseMatrix_CRS<TM01> >,
                                           DLA::MatrixD<SLA::SparseMatrix_CRS<TM02> >,

                                           DLA::MatrixD<SLA::SparseMatrix_CRS<TM10> >,
                                           DLA::MatrixD<SLA::SparseMatrix_CRS<TM11> >,
                                           DLA::MatrixD<SLA::SparseMatrix_CRS<TM12> >,

                                           DLA::MatrixD<SLA::SparseMatrix_CRS<TM20> >,
                                           DLA::MatrixD<SLA::SparseMatrix_CRS<TM21> >,
                                           DLA::MatrixD<SLA::SparseMatrix_CRS<TM22> > > >
{
  typedef BLA::MatrixBlock_3x3< DLA::MatrixD< SLA::SparseMatrixSize >,
                                DLA::MatrixD< SLA::SparseMatrixSize >,
                                DLA::MatrixD< SLA::SparseMatrixSize >,

                                DLA::MatrixD< SLA::SparseMatrixSize >,
                                DLA::MatrixD< SLA::SparseMatrixSize >,
                                DLA::MatrixD< SLA::SparseMatrixSize >,

                                DLA::MatrixD< SLA::SparseMatrixSize >,
                                DLA::MatrixD< SLA::SparseMatrixSize >,
                                DLA::MatrixD< SLA::SparseMatrixSize > > type;
};

// block 4x4 matrix
template<class SM00, class SM01, class SM02, class SM03,
         class SM10, class SM11, class SM12, class SM13,
         class SM20, class SM21, class SM22, class SM23,
         class SM30, class SM31, class SM32, class SM33>
struct MatrixSizeType<
         BLA::MatrixBlock_4x4<DLA::MatrixD<SLA::SparseMatrix_CRS<SM00> >,
                              DLA::MatrixD<SLA::SparseMatrix_CRS<SM01> >,
                              DLA::MatrixD<SLA::SparseMatrix_CRS<SM02> >,
                              DLA::MatrixD<SLA::SparseMatrix_CRS<SM03> >,
                              //
                              DLA::MatrixD<SLA::SparseMatrix_CRS<SM10> >,
                              DLA::MatrixD<SLA::SparseMatrix_CRS<SM11> >,
                              DLA::MatrixD<SLA::SparseMatrix_CRS<SM12> >,
                              DLA::MatrixD<SLA::SparseMatrix_CRS<SM13> >,
                              //
                              DLA::MatrixD<SLA::SparseMatrix_CRS<SM20> >,
                              DLA::MatrixD<SLA::SparseMatrix_CRS<SM21> >,
                              DLA::MatrixD<SLA::SparseMatrix_CRS<SM22> >,
                              DLA::MatrixD<SLA::SparseMatrix_CRS<SM23> >,
                              //
                              DLA::MatrixD<SLA::SparseMatrix_CRS<SM30> >,
                              DLA::MatrixD<SLA::SparseMatrix_CRS<SM31> >,
                              DLA::MatrixD<SLA::SparseMatrix_CRS<SM32> >,
                              DLA::MatrixD<SLA::SparseMatrix_CRS<SM33> >
                             >
                     >
{
  typedef BLA::MatrixBlock_4x4< DLA::MatrixD< SLA::SparseMatrixSize >,
                                DLA::MatrixD< SLA::SparseMatrixSize >,
                                DLA::MatrixD< SLA::SparseMatrixSize >,
                                DLA::MatrixD< SLA::SparseMatrixSize >,
                                //
                                DLA::MatrixD< SLA::SparseMatrixSize >,
                                DLA::MatrixD< SLA::SparseMatrixSize >,
                                DLA::MatrixD< SLA::SparseMatrixSize >,
                                DLA::MatrixD< SLA::SparseMatrixSize >,
                                //
                                DLA::MatrixD< SLA::SparseMatrixSize >,
                                DLA::MatrixD< SLA::SparseMatrixSize >,
                                DLA::MatrixD< SLA::SparseMatrixSize >,
                                DLA::MatrixD< SLA::SparseMatrixSize >,
                                //
                                DLA::MatrixD< SLA::SparseMatrixSize >,
                                DLA::MatrixD< SLA::SparseMatrixSize >,
                                DLA::MatrixD< SLA::SparseMatrixSize >,
                                DLA::MatrixD< SLA::SparseMatrixSize > > type;
};

// block 4x4 matrix
template<class SM00, class SM01, class SM02, class SM03,
         class SM10, class SM11, class SM12, class SM13,
         class SM20, class SM21, class SM22, class SM23,
         class SM30, class SM31, class SM32, class SM33>
struct MatrixSizeType<
         BLA::MatrixBlock_4x4<DLA::MatrixD<DLA::MatrixD<SM00> >,
                              DLA::MatrixD<DLA::MatrixD<SM01> >,
                              DLA::MatrixD<DLA::MatrixD<SM02> >,
                              DLA::MatrixD<DLA::MatrixD<SM03> >,
                              //
                              DLA::MatrixD<DLA::MatrixD<SM10> >,
                              DLA::MatrixD<DLA::MatrixD<SM11> >,
                              DLA::MatrixD<DLA::MatrixD<SM12> >,
                              DLA::MatrixD<DLA::MatrixD<SM13> >,
                              //
                              DLA::MatrixD<DLA::MatrixD<SM20> >,
                              DLA::MatrixD<DLA::MatrixD<SM21> >,
                              DLA::MatrixD<DLA::MatrixD<SM22> >,
                              DLA::MatrixD<DLA::MatrixD<SM23> >,
                              //
                              DLA::MatrixD<DLA::MatrixD<SM30> >,
                              DLA::MatrixD<DLA::MatrixD<SM31> >,
                              DLA::MatrixD<DLA::MatrixD<SM32> >,
                              DLA::MatrixD<DLA::MatrixD<SM33> >
                             >
                     >
{
  typedef BLA::MatrixBlock_4x4< DLA::MatrixD< DLA::DenseMatrixSize >,
                                DLA::MatrixD< DLA::DenseMatrixSize >,
                                DLA::MatrixD< DLA::DenseMatrixSize >,
                                DLA::MatrixD< DLA::DenseMatrixSize >,
                                //
                                DLA::MatrixD< DLA::DenseMatrixSize >,
                                DLA::MatrixD< DLA::DenseMatrixSize >,
                                DLA::MatrixD< DLA::DenseMatrixSize >,
                                DLA::MatrixD< DLA::DenseMatrixSize >,
                                //
                                DLA::MatrixD< DLA::DenseMatrixSize >,
                                DLA::MatrixD< DLA::DenseMatrixSize >,
                                DLA::MatrixD< DLA::DenseMatrixSize >,
                                DLA::MatrixD< DLA::DenseMatrixSize >,
                                //
                                DLA::MatrixD< DLA::DenseMatrixSize >,
                                DLA::MatrixD< DLA::DenseMatrixSize >,
                                DLA::MatrixD< DLA::DenseMatrixSize >,
                                DLA::MatrixD< DLA::DenseMatrixSize >> type;
};

} //namespace numpack 


#endif //MATRIXSIZETYPE_H
