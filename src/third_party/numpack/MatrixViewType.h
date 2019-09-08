// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXVIEWTYPE_H
#define MATRIXVIEWTYPE_H

#include "tools/SANSnumerics.h"

#include "DenseLinAlg/DynamicSize/MatrixD_Type.h"
#include "DenseLinAlg/StaticSize/MatrixS_Type.h"
#include "BlockLinAlg/BlockLinAlg_Type.h"
#include "SparseLinAlg/SparseLinAlg_Type.h"

namespace numpack 
{

//=============================================================================
//A template metafunction for choosing the appropriate matrix view type for a given matrix type
template< class Matrix_type >
struct MatrixViewType;

template< class T >
struct MatrixViewType< DLA::MatrixD<T> >
{
  typedef DLA::MatrixDView<T> type;
};

template< class T >
struct MatrixViewType< SLA::SparseMatrix_CRS<T> >
{
  typedef SLA::SparseMatrix_CRS<T> type;
};

template<class SM00, class SM01,
         class SM10, class SM11>
struct MatrixViewType< BLA::MatrixBlock_2x2<SM00, SM01,
                                            SM10, SM11> >
{
  typedef BLA::MatrixBlock_2x2<SM00, SM01,
                               SM10, SM11> type;
};

template<class SM00, class SM01, class SM02,
         class SM10, class SM11, class SM12,
         class SM20, class SM21, class SM22>
struct MatrixViewType< BLA::MatrixBlock_3x3< SM00, SM01, SM02,
                                             SM10, SM11, SM12,
                                             SM20, SM21, SM22 > >
{
  typedef BLA::MatrixBlock_3x3< SM00, SM01, SM02,
                                SM10, SM11, SM12,
                                SM20, SM21, SM22 > type;
};

template<class SM00, class SM01, class SM02, class SM03,
         class SM10, class SM11, class SM12, class SM13,
         class SM20, class SM21, class SM22, class SM23,
         class SM30, class SM31, class SM32, class SM33>
struct MatrixViewType< BLA::MatrixBlock_4x4<SM00, SM01, SM02, SM03,
                                            SM10, SM11, SM12, SM13,
                                            SM20, SM21, SM22, SM23,
                                            SM30, SM31, SM32, SM33> >
{
  typedef BLA::MatrixBlock_4x4<SM00, SM01, SM02, SM03,
                               SM10, SM11, SM12, SM13,
                               SM20, SM21, SM22, SM23,
                               SM30, SM31, SM32, SM33> type;
};

} //namespace numpack 


#endif //MATRIXTYPE_H
