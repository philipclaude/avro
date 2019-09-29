// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef dense_DOT_H
#define dense_DOT_H

#include "tools/SANSnumerics.h"
#include "tools/minmax.h"

#include "PromoteSurreal.h"

#include "numpack/dense/static/MatrixS_Type.h"
#include "numpack/dense/dynamic/MatrixD_Type.h"

#include "tools/SANSTraitsScalar.h"
#include "numpack/sparse/tools/sparse_Scalar.h"

namespace numpack 
{

namespace SLA
{
template< class Derived, bool useRF >
class sparseType;
}

template<class ExprL, class ExprR>
inline Real
dot(const SLA::sparseType<ExprL, false>& eL, const SLA::sparseType<ExprR, false>& eR);

//===========================================================================//
// scalar 'dot' products
//===========================================================================//

inline Real dot(const Real& a, const Real& b) { return a*b; }

template<int N>
inline SurrealS<N> dot(const SurrealS<N>& a, const SurrealS<N>& b) { return a*b; }

template<int N>
inline SurrealS<N> dot(const SurrealS<N>& a, const Real& b) { return a*b; }

template<int N>
inline SurrealS<N> dot(const Real& a, const SurrealS<N>& b) { return a*b; }


template<int M, class T2>
inline DLA::OpMulSScalar<DLA::MatrixS<M,1,T2>, Real, false, true >
dot(const Real& a, const DLA::MatrixS<M,1,T2>& b) { return a*b; }

template<int N, int M, class T2>
inline DLA::OpMulSScalar<DLA::MatrixS<M,1,T2>, SurrealS<N>, false, true >
dot(const SurrealS<N>& a, const DLA::MatrixS<M,1,T2>& b) { return a*b; }

template<int M, class T2>
inline DLA::OpMulSScalar<DLA::MatrixSymS<M,T2>, Real, true, false >
dot(const Real& a, const DLA::MatrixSymS<M,T2>& b) { return a*b; }

template<int M, class T2>
inline DLA::OpMulSScalar<DLA::MatrixSymS<M,T2>, Real, true, false >
dot(const DLA::MatrixSymS<M,T2>& a, const Real& b ) { return a*b; }

//===========================================================================//
// MatrixS dot products
//===========================================================================//

// Used to determine the appropriate scalar type returned from a dot product
// The 'scalar' might be a VectorS, e.g.
// dot( VectorS<Real>, VectorS<VectorS> )
// only performes the dot across the outermost VectorS
template<class T1, class T2>
struct dot_Scalar
{
  typedef typename Scalar< typename promote_Surreal< T1, T2 >::type >::type type;
};

template<int M2, class T2>
struct dot_Scalar<Real, DLA::VectorS<M2,T2>> { typedef DLA::VectorS<M2,T2> type; };

template<int M1, class T1>
struct dot_Scalar<DLA::VectorS<M1,T1>, Real> { typedef DLA::VectorS<M1,T1> type; };

template<int N, class T>
struct dot_Scalar<Real, DLA::MatrixSymS<N, T>> { typedef DLA::MatrixSymS<N, T> type; };

template<int N, class T>
struct dot_Scalar<DLA::MatrixSymS<N, T>, Real> { typedef DLA::MatrixSymS<N, T> type; };

template<int SN, class T1, int M2, class T2>
struct dot_Scalar<SurrealS<SN,T1>, DLA::VectorS<M2,T2>>
{
  typedef typename promote_Surreal<SurrealS<SN,T1>, T2>::type T;
  typedef DLA::VectorS<M2,T> type;
};

template<int M1, class T1, int SN, class T2>
struct dot_Scalar<DLA::VectorS<M1,T1>, SurrealS<SN,T2>>
{
  typedef typename promote_Surreal<SurrealS<SN,T2>, T1>::type T;
  typedef DLA::VectorS<M1,T> type;
};


//---------------------------------------------------------------------------//
template<int M, class T1, class T2>
inline typename dot_Scalar< T1, T2 >::type
dot(const DLA::MatrixS<M,1,T1>& a, const DLA::MatrixS<M,1,T2>& b)
{
  typename dot_Scalar< T1, T2 >::type val = 0;
  for (int i = 0; i < M; i++)
    val += dot(a(i,0),b(i,0));

  return val;
}

//---------------------------------------------------------------------------//
template<int M, class T>
inline typename Scalar< T >::type
dot(const DLA::MatrixS<M,1,T>& a, const DLA::MatrixS<M,1,T>& b)
{
  typename Scalar< T >::type val = 0;
  for (int i = 0; i < M; i++)
    val += dot(a(i,0),b(i,0));

  return val;
}

//---------------------------------------------------------------------------//
template<class ExprL, class ExprR>
inline typename dot_Scalar< typename ExprL::Ttype, typename ExprR::Ttype >::type
dot(const DLA::MatrixSType<ExprL, true, true>& eL, const DLA::MatrixSType<ExprR, true, true>& eR)
{
  typedef typename dot_Scalar< typename ExprL::Ttype, typename ExprR::Ttype >::type T;

  static_assert(ExprL::M == 1 || ExprL::N == 1, "Matrix expressions must be a vector");
  static_assert(ExprR::M == 1 || ExprR::N == 1, "Matrix expressions must be a vector");

  static_assert(ExprL::M == ExprR::M ||
                ExprL::N == ExprR::N ||
                ExprL::M == ExprR::N, "Matrix expressions must have same number of rows or columns");

  static const int M = MAX(ExprL::M, ExprL::N);

  DLA::MatrixS<ExprL::M,ExprL::N,typename ExprL::Ttype> a(eL);
  DLA::MatrixS<ExprR::M,ExprR::N,typename ExprR::Ttype> b(eR);

  T val = 0;
  for (int i = 0; i < M; i++)
    val += dot(a(i,0),b(i,0));

  return val;
}

//---------------------------------------------------------------------------//
template<class ExprL, class ExprR>
inline typename dot_Scalar< typename ExprL::Ttype, typename ExprR::Ttype >::type
dot(const DLA::MatrixSType<ExprL, true, true>& eL, const DLA::MatrixSType<ExprR, false, true>& eR)
{
  typedef typename dot_Scalar< typename ExprL::Ttype, typename ExprR::Ttype >::type T;

  static_assert(ExprL::M == 1 || ExprL::N == 1, "Matrix expressions must be a vector");
  static_assert(ExprR::M == 1 || ExprR::N == 1, "Matrix expressions must be a vector");

  static_assert(ExprL::M == ExprR::M ||
                ExprL::N == ExprR::N ||
                ExprL::M == ExprR::N, "Matrix expressions must have same number of rows or columns");

  static const int M = MAX(ExprL::M, ExprL::N);

  DLA::MatrixS<ExprL::M,ExprL::N,typename ExprL::Ttype> a(eL);
  const ExprR& b = eR.cast();

  T val = 0;
  for (int i = 0; i < M; i++)
    val += dot(a(i,0),b.value(i));

  return val;
}

//---------------------------------------------------------------------------//
template<class ExprL, class ExprR>
inline typename dot_Scalar< typename ExprL::Ttype, typename ExprR::Ttype >::type
dot(const DLA::MatrixSType<ExprL, false, true>& eL, const DLA::MatrixSType<ExprR, true, true>& eR)
{
  typedef typename dot_Scalar< typename ExprL::Ttype, typename ExprR::Ttype >::type T;

  static_assert(ExprL::M == 1 || ExprL::N == 1, "Matrix expressions must be a vector");
  static_assert(ExprR::M == 1 || ExprR::N == 1, "Matrix expressions must be a vector");

  static_assert(ExprL::M == ExprR::M ||
                ExprL::N == ExprR::N ||
                ExprL::M == ExprR::N, "Matrix expressions must have same number of rows or columns");

  static const int M = MAX(ExprL::M, ExprL::N);

  const ExprL& a = eL.cast();
  DLA::MatrixS<ExprR::M,ExprR::N,typename ExprR::Ttype> b(eR);

  T val = 0;
  for (int i = 0; i < M; i++)
    val += dot(a.value(i),b(i,0));

  return val;
}

//---------------------------------------------------------------------------//
template<class ExprL, class ExprR>
inline typename dot_Scalar< typename ExprL::Ttype, typename ExprR::Ttype >::type
dot(const DLA::MatrixSType<ExprL, false, true>& eL, const DLA::MatrixSType<ExprR, false, true>& eR)
{
  typedef typename dot_Scalar< typename ExprL::Ttype, typename ExprR::Ttype >::type T;

  static_assert(ExprL::M == 1 || ExprL::N == 1, "Matrix expressions must be a vector");
  static_assert(ExprR::M == 1 || ExprR::N == 1, "Matrix expressions must be a vector");

  static_assert(ExprL::M == ExprR::M ||
                ExprL::N == ExprR::N ||
                ExprL::M == ExprR::N, "Matrix expressions must have same number of rows or columns");

  static const int M = MAX(ExprL::M, ExprL::N);

  const ExprL& a = eL.cast();
  const ExprR& b = eR.cast();

  T val = 0;
  for (int i = 0; i < M; i++)
    val += dot(a.value(i),b.value(i));

  return val;
}


//===========================================================================//
// MatrixD dot products
//===========================================================================//
template<class T1, class T2>
inline typename dot_Scalar< T1, T2 >::type
dot(const DLA::VectorDView<T1>& a, const DLA::VectorDView<T2>& b)
{
  typedef typename dot_Scalar< T1, T2 >::type T;

  SANS_ASSERT(a.size() == b.size());

  const int size = a.size();

  T val = 0;
  for (int i = 0; i < size; i++)
    val += dot(a[i],b[i]);

  return val;
}

//---------------------------------------------------------------------------//
template<class T1, class T2>
inline typename dot_Scalar< T1, T2 >::type
dot(const DLA::MatrixDView<T1>& a, const DLA::MatrixDView<T2>& b)
{
  typedef typename dot_Scalar< T1, T2 >::type T;

  SANS_ASSERT(a.size() == b.size());
  SANS_ASSERT(a.n() == b.n());
  SANS_ASSERT(a.m() == b.m());
  SANS_ASSERT(a.m() == 1 || a.n() == 1);
  SANS_ASSERT(b.m() == 1 || b.n() == 1);

  const int size = a.size();

  T val = 0;
  for (int i = 0; i < size; i++)
    val += dot(a(i,0),b(i,0));

  return val;
}

//---------------------------------------------------------------------------//
template<class ExprL, bool useRFL, class ExprR, bool useRFR>
inline typename dot_Scalar< typename ExprL::node_type, typename ExprR::node_type >::type
dot(const DLA::MatrixDType<ExprL, useRFL>& eL, const DLA::MatrixDType<ExprR, useRFR>& eR)
{
  typedef typename dot_Scalar< typename ExprL::node_type, typename ExprR::node_type >::type T;
  const ExprL& a = eL.cast();
  const ExprR& b = eR.cast();

  SANS_ASSERT(a.size() == b.size());
  SANS_ASSERT(a.n() == b.n());
  SANS_ASSERT(a.m() == b.m());
  SANS_ASSERT(a.m() == 1 || a.n() == 1);
  SANS_ASSERT(b.m() == 1 || b.n() == 1);

  const int size = a.size();

  T val = 0;
  for (int i = 0; i < size; i++)
    val += dot(a.value(i),b.value(i));

  return val;
}

}

#endif //dense_DOT_H
