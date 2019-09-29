// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define FGMRES_INSTANTIATE
#include "FGMRES_impl.h"

#include "FGMRES.h"

#include "numpack/sparse/tools/Array.h"
#include "numpack/sparse/tools/UpperTriangMatrix.h"
#include "numpack/sparse/sparse_Add.h"
#include "numpack/sparse/sparse_Mul.h"

#include "numpack/dense/tools/norm.h"
#include "numpack/dense/tools/dot.h"
#include "numpack/sparse/tools/norm.h"
#include "numpack/sparse/tools/dot.h"

#include "numpack/block/MatrixBlock_4x4.h"

#include "tools/minmax.h"


namespace numpack 
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
typedef DLA::MatrixS<1,8,Real> MatrixS18;
typedef DLA::MatrixS<2,1,Real> MatrixS21;
typedef DLA::MatrixS<2,2,Real> MatrixS22;
typedef DLA::MatrixS<2,8,Real> MatrixS28;
typedef DLA::MatrixS<8,1,Real> MatrixS81;
typedef DLA::MatrixS<8,2,Real> MatrixS82;
typedef DLA::MatrixS<8,8,Real> MatrixS88;

typedef DLA::MatrixS<1,3,Real> MatrixS13;
typedef DLA::MatrixS<3,1,Real> MatrixS31;
typedef DLA::MatrixS<3,3,Real> MatrixS33;
typedef DLA::MatrixS<3,8,Real> MatrixS38;
typedef DLA::MatrixS<8,3,Real> MatrixS83;

typedef DLA::MatrixS<1,4,Real> MatrixS14;
typedef DLA::MatrixS<4,1,Real> MatrixS41;
typedef DLA::MatrixS<4,4,Real> MatrixS44;
typedef DLA::MatrixS<4,8,Real> MatrixS48;
typedef DLA::MatrixS<8,4,Real> MatrixS84;

typedef DLA::MatrixS<1,6,Real> MatrixS16;
typedef DLA::MatrixS<6,1,Real> MatrixS61;
typedef DLA::MatrixS<6,6,Real> MatrixS66;
typedef DLA::MatrixS<6,8,Real> MatrixS68;
typedef DLA::MatrixS<8,6,Real> MatrixS86;

typedef BLA::MatrixBlock_4x4<
          DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS88> >, DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS81> >,
          DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS82> >, DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS81> >,
          //
          DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS18> >, DLA::MatrixD<SLA::SparseMatrix_CRS<Real> >,
          DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS12> >, DLA::MatrixD<SLA::SparseMatrix_CRS<Real> >,
          //
          DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS28> >, DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS21> >,
          DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS22> >, DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS21> >,
          //
          DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS18> >, DLA::MatrixD<SLA::SparseMatrix_CRS<Real> >,
          DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS12> >, DLA::MatrixD<SLA::SparseMatrix_CRS<Real> >
                            > BlockMatrix4x4_Coupling_Auxv_Auxi_IBL2D_Panel;

typedef BLA::MatrixBlock_4x4<
          DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS88> >, DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS81> >,
          DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS83> >, DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS81> >,
          //
          DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS18> >, DLA::MatrixD<SLA::SparseMatrix_CRS<Real> >,
          DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS13> >, DLA::MatrixD<SLA::SparseMatrix_CRS<Real> >,
          //
          DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS38> >, DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS31> >,
          DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS33> >, DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS31> >,
          //
          DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS18> >, DLA::MatrixD<SLA::SparseMatrix_CRS<Real> >,
          DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS13> >, DLA::MatrixD<SLA::SparseMatrix_CRS<Real> >
                            > BlockMatrix4x4_Coupling_Auxv_Auxi_IBL2D3Unknown_Panel;

typedef BLA::MatrixBlock_4x4<
          DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS88> >, DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS81> >,
          DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS84> >, DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS81> >,
          //
          DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS18> >, DLA::MatrixD<SLA::SparseMatrix_CRS<Real> >,
          DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS14> >, DLA::MatrixD<SLA::SparseMatrix_CRS<Real> >,
          //
          DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS48> >, DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS41> >,
          DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS44> >, DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS41> >,
          //
          DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS18> >, DLA::MatrixD<SLA::SparseMatrix_CRS<Real> >,
          DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS14> >, DLA::MatrixD<SLA::SparseMatrix_CRS<Real> >
                            > BlockMatrix4x4_Coupling_Auxv_Auxi_IBL2DtransitionLag_Panel;

typedef BLA::MatrixBlock_4x4<
          DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS88> >, DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS81> >,
          DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS86> >, DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS81> >,
          //
          DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS18> >, DLA::MatrixD<SLA::SparseMatrix_CRS<Real> >,
          DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS16> >, DLA::MatrixD<SLA::SparseMatrix_CRS<Real> >,
          //
          DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS68> >, DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS61> >,
          DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS66> >, DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS61> >,
          //
          DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS18> >, DLA::MatrixD<SLA::SparseMatrix_CRS<Real> >,
          DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixS16> >, DLA::MatrixD<SLA::SparseMatrix_CRS<Real> >
                            > BlockMatrix4x4_Coupling_Auxv_Auxi_IBL2DtransitionLag6unknownCutCell_Panel;

INSTANTIATE( BlockMatrix4x4_Coupling_Auxv_Auxi_IBL2D_Panel )
INSTANTIATE( BlockMatrix4x4_Coupling_Auxv_Auxi_IBL2D3Unknown_Panel )
INSTANTIATE( BlockMatrix4x4_Coupling_Auxv_Auxi_IBL2DtransitionLag_Panel )
INSTANTIATE( BlockMatrix4x4_Coupling_Auxv_Auxi_IBL2DtransitionLag6unknownCutCell_Panel )


} //namespace SLA
} //namespace numpack 
