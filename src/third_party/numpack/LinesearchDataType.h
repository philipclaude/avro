// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef LINESEARCHVECTORTYPE_H
#define LINESEARCHVECTORTYPE_H

#include <memory> //shared_ptr
#include <array>

#include "tools/SANSnumerics.h"

#include "DenseLinAlg/DynamicSize/MatrixD_Type.h"
#include "BlockLinAlg/BlockLinAlg_Type.h"
#include "SparseLinAlg/SparseLinAlg_Type.h"

namespace numpack 
{

//=============================================================================
#if 1
template< class SystemVector >
struct LinesearchDataType
{
  typedef DLA::VectorD<SLA::SparseVector<DLA::VectorD<Real>>> VType;
  typedef std::shared_ptr<VType> type;
};
#else
//TODO: this should be the way to go, but can't figure out why it fails now --> block 2x2 not compiling for some cases
//Metafunction for choosing the appropriate line-search data type for a given SystemVector type
template<class SystemVector>
struct LinesearchDataType;

// specializations for different system vector types
//

template< >
struct LinesearchDataType< DLA::VectorD<DLA::VectorD<Real> > >
{
  typedef DLA::VectorD<DLA::VectorD<Real> > VType;
  typedef std::shared_ptr<VType> type;
};

template< >
struct LinesearchDataType< DLA::VectorD<SLA::SparseVector<Real> > >
{
  typedef DLA::VectorD<SLA::SparseVector<Real> > VType;
  typedef std::shared_ptr<VType> type;
};

template< >
struct LinesearchDataType< DLA::VectorD<SLA::SparseVector<DLA::VectorD<Real> > > >
{
  typedef DLA::VectorD<SLA::SparseVector<DLA::VectorD<Real> > > VType;
  typedef std::shared_ptr<VType> type;
};
#endif

//Specialization for 2-block vector
template< class T0, class T1 >
struct LinesearchDataType< BLA::VectorBlock_2<T0, T1> >
{
  typedef DLA::VectorD<SLA::SparseVector<DLA::VectorD<Real>>> VType;
  typedef std::array<std::shared_ptr<VType>, 2> type;
};

//Specialization for 4-block vector
template< class T0, class T1, class T2, class T3 >
struct LinesearchDataType< BLA::VectorBlock_4<T0,T1,T2,T3> >
{
  typedef DLA::VectorD<SLA::SparseVector<DLA::VectorD<Real>>> VType;
  typedef std::array<std::shared_ptr<VType>, 4> type;
};

}

#endif //LINESEARCHVECTORTYPE_H
