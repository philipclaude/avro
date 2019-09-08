// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define SCALARMATRIX_CRS_INSTANTIATE
#include "ScalarMatrix_CRS_impl.h"

namespace SANS
{
namespace SLA
{
template class ScalarMatrix_CRS<int>;

template ScalarMatrix_CRS<int>::ScalarMatrix_CRS(const SparseMatrix_CRS<Real>&);
template ScalarMatrix_CRS<int>::ScalarMatrix_CRS(const DLA::MatrixDView<SparseMatrix_CRS<Real> >&);

template ScalarMatrix_CRS<int>::ScalarMatrix_CRS(const SparseMatrix_CRS< DLA::MatrixS<2, 2, Real> >&);
template ScalarMatrix_CRS<int>::ScalarMatrix_CRS(const SparseMatrix_CRS< DLA::MatrixS<3, 3, Real> >&);
template ScalarMatrix_CRS<int>::ScalarMatrix_CRS(const SparseMatrix_CRS< DLA::MatrixS<8, 8, Real> >&);

template ScalarMatrix_CRS<int>::ScalarMatrix_CRS(const DLA::MatrixDView< SparseMatrix_CRS< DLA::MatrixS<2, 2, Real> > >&);
template ScalarMatrix_CRS<int>::ScalarMatrix_CRS(const DLA::MatrixDView< SparseMatrix_CRS< DLA::MatrixS<3, 3, Real> > >&);
template ScalarMatrix_CRS<int>::ScalarMatrix_CRS(const DLA::MatrixDView< SparseMatrix_CRS< DLA::MatrixS<4, 4, Real> > >&);
template ScalarMatrix_CRS<int>::ScalarMatrix_CRS(const DLA::MatrixDView< SparseMatrix_CRS< DLA::MatrixS<5, 5, Real> > >&);
template ScalarMatrix_CRS<int>::ScalarMatrix_CRS(const DLA::MatrixDView< SparseMatrix_CRS< DLA::MatrixS<6, 6, Real> > >&);
template ScalarMatrix_CRS<int>::ScalarMatrix_CRS(const DLA::MatrixDView< SparseMatrix_CRS< DLA::MatrixS<7, 7, Real> > >&);
template ScalarMatrix_CRS<int>::ScalarMatrix_CRS(const DLA::MatrixDView< SparseMatrix_CRS< DLA::MatrixS<8, 8, Real> > >&);

template ScalarMatrix_CRS<int>::ScalarMatrix_CRS(const DLA::MatrixDView< SparseMatrix_CRS< DLA::MatrixS<1, 2, Real> > >&);
template ScalarMatrix_CRS<int>::ScalarMatrix_CRS(const DLA::MatrixDView< SparseMatrix_CRS< DLA::MatrixS<1, 8, Real> > >&);
template ScalarMatrix_CRS<int>::ScalarMatrix_CRS(const DLA::MatrixDView< SparseMatrix_CRS< DLA::MatrixS<2, 1, Real> > >&);
template ScalarMatrix_CRS<int>::ScalarMatrix_CRS(const DLA::MatrixDView< SparseMatrix_CRS< DLA::MatrixS<2, 8, Real> > >&);
template ScalarMatrix_CRS<int>::ScalarMatrix_CRS(const DLA::MatrixDView< SparseMatrix_CRS< DLA::MatrixS<8, 1, Real> > >&);
template ScalarMatrix_CRS<int>::ScalarMatrix_CRS(const DLA::MatrixDView< SparseMatrix_CRS< DLA::MatrixS<8, 2, Real> > >&);

template ScalarMatrix_CRS<int>::ScalarMatrix_CRS(
    const BLA::MatrixBlock_2x2< DLA::MatrixD<SparseMatrix_CRS<Real> >, DLA::MatrixD<SparseMatrix_CRS<Real> >,
                                DLA::MatrixD<SparseMatrix_CRS<Real> >, DLA::MatrixD<SparseMatrix_CRS<Real> > >&);

template ScalarMatrix_CRS<int>::ScalarMatrix_CRS(
    const BLA::MatrixBlock_2x2< DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<2,2,Real>>>, DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<2,1,Real>>>,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<1,2,Real>>>, DLA::MatrixD<SparseMatrix_CRS<Real>>                  >&);

template ScalarMatrix_CRS<int>::ScalarMatrix_CRS(
    const BLA::MatrixBlock_2x2< DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<2,2,Real>>>, DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<2,2,Real>>>,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<2,2,Real>>>, DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<2,2,Real>>> >&);

template ScalarMatrix_CRS<int>::ScalarMatrix_CRS(
    const BLA::MatrixBlock_2x2< DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<3,3,Real>>>, DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<3,2,Real>>>,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<2,3,Real>>>, DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<2,2,Real>>> >&);

template ScalarMatrix_CRS<int>::ScalarMatrix_CRS(
    const BLA::MatrixBlock_2x2< DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<8,8,Real>>>, DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<8,2,Real>>>,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<2,8,Real>>>, DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<2,2,Real>>> >&);

template ScalarMatrix_CRS<int>::ScalarMatrix_CRS(
    const BLA::MatrixBlock_2x2< DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<3,3,Real>>>, DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<3,1,Real>>>,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<1,3,Real>>>, DLA::MatrixD<SparseMatrix_CRS<Real>>                  >&);

template ScalarMatrix_CRS<int>::ScalarMatrix_CRS(
    const BLA::MatrixBlock_2x2< DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<4,4,Real>>>, DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<4,1,Real>>>,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<1,4,Real>>>, DLA::MatrixD<SparseMatrix_CRS<Real>>                  >&);

template ScalarMatrix_CRS<int>::ScalarMatrix_CRS(
    const BLA::MatrixBlock_3x3< DLA::MatrixD<SparseMatrix_CRS<Real>                  >,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<1,3,Real>>>,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<1,2,Real>>>,

                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<3,1,Real>>>,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<3,3,Real>>>,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<3,2,Real>>>,

                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<2,1,Real>>>,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<2,3,Real>>>,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<2,2,Real>>> >& );

template ScalarMatrix_CRS<int>::ScalarMatrix_CRS(
    const BLA::MatrixBlock_4x4< DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<4,4,Real>>>,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<4,3,Real>>>,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<4,2,Real>>>,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<4,1,Real>>>,

                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<3,4,Real>>>,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<3,3,Real>>>,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<3,2,Real>>>,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<3,1,Real>>>,

                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<2,4,Real>>>,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<2,3,Real>>>,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<2,2,Real>>>,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<2,1,Real>>>,

                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<1,4,Real>>>,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<1,3,Real>>>,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<1,2,Real>>>,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<1,1,Real>>>>& );

template ScalarMatrix_CRS<int>::ScalarMatrix_CRS(
    const BLA::MatrixBlock_4x4< DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<4,4,Real>>>,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<4,3,Real>>>,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<4,2,Real>>>,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<4,2,Real>>>,

                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<3,4,Real>>>,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<3,3,Real>>>,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<3,2,Real>>>,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<3,2,Real>>>,

                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<2,4,Real>>>,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<2,3,Real>>>,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<2,2,Real>>>,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<2,2,Real>>>,

                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<2,4,Real>>>,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<2,3,Real>>>,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<2,2,Real>>>,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<2,2,Real>>>>& );

template ScalarMatrix_CRS<int>::ScalarMatrix_CRS(
    const BLA::MatrixBlock_4x4< DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<8,8,Real> > >,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<8,1,Real> > >,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<8,2,Real> > >,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<8,1,Real> > >,

                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<1,8,Real> > >,
                                DLA::MatrixD<SparseMatrix_CRS<Real> >,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<1,2,Real> > >,
                                DLA::MatrixD<SparseMatrix_CRS<Real> >,

                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<2,8,Real> > >,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<2,1,Real> > >,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<2,2,Real> > >,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<2,1,Real> > >,

                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<1,8,Real> > >,
                                DLA::MatrixD<SparseMatrix_CRS<Real> >,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<1,2,Real> > >,
                                DLA::MatrixD<SparseMatrix_CRS<Real> > >& );

template ScalarMatrix_CRS<int>::ScalarMatrix_CRS(
    const BLA::MatrixBlock_4x4< DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<8,8,Real> > >,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<8,1,Real> > >,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<8,3,Real> > >,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<8,1,Real> > >,

                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<1,8,Real> > >,
                                DLA::MatrixD<SparseMatrix_CRS<Real> >,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<1,3,Real> > >,
                                DLA::MatrixD<SparseMatrix_CRS<Real> >,

                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<3,8,Real> > >,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<3,1,Real> > >,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<3,3,Real> > >,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<3,1,Real> > >,

                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<1,8,Real> > >,
                                DLA::MatrixD<SparseMatrix_CRS<Real> >,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<1,3,Real> > >,
                                DLA::MatrixD<SparseMatrix_CRS<Real> > >& );

template ScalarMatrix_CRS<int>::ScalarMatrix_CRS(
    const BLA::MatrixBlock_4x4< DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<8,8,Real> > >,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<8,1,Real> > >,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<8,4,Real> > >,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<8,1,Real> > >,

                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<1,8,Real> > >,
                                DLA::MatrixD<SparseMatrix_CRS<Real> >,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<1,4,Real> > >,
                                DLA::MatrixD<SparseMatrix_CRS<Real> >,

                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<4,8,Real> > >,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<4,1,Real> > >,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<4,4,Real> > >,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<4,1,Real> > >,

                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<1,8,Real> > >,
                                DLA::MatrixD<SparseMatrix_CRS<Real> >,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<1,4,Real> > >,
                                DLA::MatrixD<SparseMatrix_CRS<Real> > >& );

template ScalarMatrix_CRS<int>::ScalarMatrix_CRS(
    const BLA::MatrixBlock_4x4< DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<8,8,Real> > >,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<8,1,Real> > >,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<8,6,Real> > >,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<8,1,Real> > >,

                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<1,8,Real> > >,
                                DLA::MatrixD<SparseMatrix_CRS<Real> >,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<1,6,Real> > >,
                                DLA::MatrixD<SparseMatrix_CRS<Real> >,

                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<6,8,Real> > >,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<6,1,Real> > >,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<6,6,Real> > >,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<6,1,Real> > >,

                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<1,8,Real> > >,
                                DLA::MatrixD<SparseMatrix_CRS<Real> >,
                                DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<1,6,Real> > >,
                                DLA::MatrixD<SparseMatrix_CRS<Real> > >& );
}
}
