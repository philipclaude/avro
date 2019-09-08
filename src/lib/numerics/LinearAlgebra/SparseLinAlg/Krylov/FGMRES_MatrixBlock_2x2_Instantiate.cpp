// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define FGMRES_INSTANTIATE
#include "FGMRES_impl.h"

#include "FGMRES.h"

#include "LinearAlgebra/SparseLinAlg/tools/Array.h"
#include "LinearAlgebra/SparseLinAlg/tools/UpperTriangMatrix.h"
#include "LinearAlgebra/SparseLinAlg/SparseLinAlg_Add.h"
#include "LinearAlgebra/SparseLinAlg/SparseLinAlg_Mul.h"

#include "LinearAlgebra/DenseLinAlg/tools/norm.h"
#include "LinearAlgebra/DenseLinAlg/tools/dot.h"
#include "LinearAlgebra/SparseLinAlg/tools/norm.h"
#include "LinearAlgebra/SparseLinAlg/tools/dot.h"

#include "LinearAlgebra/BlockLinAlg/MatrixBlock_2x2.h"

#include "tools/minmax.h"


namespace SANS
{
namespace SLA
{

// define instantiation function using macros

#define INSTANTIATE(BLOCKMATRIX) \
template<> \
FGMRES<BLOCKMATRIX>::FGMRES( const PyDict& d, Solver_ptr M ) :\
  Base_type(M->systemSolve()), \
  params(FGMRESParam::params), \
  M_(M), \
  tol_(d.get(params.tol)), \
  tolFinal_(0), \
  nInner_(d.get(params.nInner)), \
  nOuter_(d.get(params.nOuter)), \
  printConv_(d.get(params.PrintCovergence)), \
  itTotal_(0) \
{ \
  SANS_DEVELOPER_EXCEPTION("FGMRES<BlockMatrix> not implemented"); \
} \
\
template<> \
FGMRES<BLOCKMATRIX>::~FGMRES() {} \
\
template<> \
void FGMRES<BLOCKMATRIX>::factorize() \
{ \
  SANS_DEVELOPER_EXCEPTION("FGMRES<BlockMatrix> not implemented"); \
} \
\
template<> \
LinearSolveStatus FGMRES<BLOCKMATRIX>::backsolve( const SparseVectorView_type& b, SparseVectorView_type& x ) const \
{ \
  SANS_DEVELOPER_EXCEPTION("FGMRES<BlockMatrix> not implemented"); return LinearSolveStatus(); \
}

// instantiations
typedef DLA::MatrixS<1,2,Real> MatrixS12;
typedef DLA::MatrixS<2,1,Real> MatrixS21;
typedef DLA::MatrixS<2,2,Real> MatrixS22;
typedef DLA::MatrixS<2,8,Real> MatrixS28;
typedef DLA::MatrixS<8,2,Real> MatrixS82;
typedef DLA::MatrixS<8,8,Real> MatrixS88;

typedef BLA::MatrixBlock_2x2< DLA::MatrixD<SLA::SparseMatrix_CRS<Real>>, DLA::MatrixD<SLA::SparseMatrix_CRS<Real>>,
                              DLA::MatrixD<SLA::SparseMatrix_CRS<Real>>, DLA::MatrixD<SLA::SparseMatrix_CRS<Real>> > BlockRealReal;

INSTANTIATE(BlockRealReal)

typedef BLA::MatrixBlock_2x2< DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS22>>, DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS21>>,
                              DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS12>>, DLA::MatrixD<SLA::SparseMatrix_CRS<Real>> > BlockMatrixS22Real;

INSTANTIATE(BlockMatrixS22Real)

typedef BLA::MatrixBlock_2x2< DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS22>>, DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS22>>,
                              DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS22>>, DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS22>> > BlockMatrixS22;

INSTANTIATE(BlockMatrixS22)

typedef BLA::MatrixBlock_2x2< DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS88>>, DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS88>>,
                              DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS88>>, DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS88>> > BlockMatrixS88;

INSTANTIATE(BlockMatrixS88)

typedef BLA::MatrixBlock_2x2< DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS88>>, DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS82>>,
                              DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS28>>, DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS22>> >
        BlockMatrixS88MatrixS22;

INSTANTIATE(BlockMatrixS88MatrixS22)

typedef DLA::MatrixS<3,3,Real> MatrixQ33;
typedef DLA::MatrixS<3,1,Real> MatrixQ31;
typedef DLA::MatrixS<1,3,Real> MatrixQ13;
typedef DLA::MatrixS<1,1,Real> MatrixQ11;

typedef BLA::MatrixBlock_2x2<DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixQ33> >, DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixQ33> >,
                             DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixQ33> >, DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixQ33> > > BlockMatrix33;

INSTANTIATE(BlockMatrix33)

typedef BLA::MatrixBlock_2x2<DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixQ11> >, DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixQ11> >,
                             DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixQ11> >, DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixQ11> > > BlockMatrix11;

INSTANTIATE(BlockMatrix11)

typedef BLA::MatrixBlock_2x2<DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixQ33> >, DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixQ31> >,
                             DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixQ13> >, DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixQ11> > > BlockMatrix3311;

INSTANTIATE(BlockMatrix3311)

typedef BLA::MatrixBlock_2x2<DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixQ33> >, DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixQ31> >,
                             DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixQ13> >, DLA::MatrixD<SLA::SparseMatrix_CRS<Real> > > BlockMatrix33RealReal;

INSTANTIATE(BlockMatrix33RealReal)

typedef DLA::MatrixS<4,4,Real> MatrixQ44;
typedef DLA::MatrixS<4,1,Real> MatrixQ41;
typedef DLA::MatrixS<1,4,Real> MatrixQ14;

typedef BLA::MatrixBlock_2x2<DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixQ44> >, DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixQ41> >,
                             DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixQ14> >, DLA::MatrixD<SLA::SparseMatrix_CRS<Real> > > BlockMatrix44RealReal;

INSTANTIATE(BlockMatrix44RealReal)

} //namespace SLA
} //namespace SANS
